#include "WirechamberConstruction.hh"

#include <math.h>

#include <G4PVReplica.hh>
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4EqMagElectricField.hh>
#include <G4ClassicalRK4.hh>

WirechamberConstruction::WirechamberConstruction()
: dWindowThick(6*um), dMWPCEntranceR(7.0*cm), dMWPCExitR(7.5*cm),
mMWPCActiveRegionGas(WCPentane), dEntranceToCathodes(5.0*mm), dExitToCathodes(5.0*mm), dE0(0), fMyBField(NULL)
{
  // can again do stuff in this constructor. But class members already set.
}

void WirechamberConstruction::Build(int side)
{
  //----- construct active gas volume using supplementary class
  ActiveRegion.mMWPCGas = mMWPCActiveRegionGas;
  ActiveRegion.Build(side);

  dMWPCContainerHalf_Z = 0.5*(dEntranceToCathodes + dExitToCathodes + 2*cm);
  G4double gasVolumeWidth = 8.0*inch;	// MWPC gas box width

  dd = ActiveRegion.dWireSpacing;
  dL = ActiveRegion.dPlaneSpacing;
  dr = ActiveRegion.dAnodeRadius;

  //----- container volume for all MWPC
  G4Box* containerBox = new G4Box("mwpc_container_box", gasVolumeWidth/2., gasVolumeWidth/2., dMWPCContainerHalf_Z);
  container_log = new G4LogicalVolume(containerBox, mMWPCActiveRegionGas, Append(side, "mwpc_container_log_"));
  container_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // MWPC active gas volume placement with wireplane, relative to MWPC container volume
  vMyTranslation = G4ThreeVector(0, 0, (dEntranceToCathodes - dExitToCathodes)/2.);
  // Note: these lines place WirechamberActiveRegion inside mwpc.
  new G4PVPlacement(NULL, vMyTranslation, ActiveRegion.gas_log, Append(side, "mwpc_activeReg_phys_"), container_log, false, 0);

  //----- Kevlar strings and Mylar windows
  // rectangular cross section strings with equal volume to norminal 140um cylinders
  G4double kevRadius = 0.07*mm;
  G4double kevSpacing = 5.0*mm;
  G4int NbKevWires = 32;
  G4double kevLength = 15.0*cm;
  double kevAspectRatio = 16.0;    // aspect ratio, width:depth.
  G4double kevArea = M_PI*kevRadius*kevRadius;	// total cross section area
  G4double kevEffWidth = sqrt(kevArea*kevAspectRatio);	// effective width
  G4double kevEffThick = sqrt(kevArea/kevAspectRatio);	// effective thickness
  G4double kev_PosZ = -dMWPCContainerHalf_Z + kevEffThick/2.;

  // create shapes for kevlar strings. Will be replicated into an array
  G4Box* kevContainerBox = new G4Box("kevContainer_box", NbKevWires*kevSpacing/2., kevLength/2., kevEffThick/2.);
  G4Box* kevSegBox = new G4Box("kevSeg_box", kevSpacing/2., kevLength/2., kevEffThick/2.);
  G4Box* kevStripBox = new G4Box("kevStrip_box", kevEffWidth/2., kevLength/2., kevEffThick/2.);

  // Shapes for the Mylar windows in the MWPC
  G4Tubs* winInnerTube = new G4Tubs("winInnerTube", 0., dMWPCEntranceR, dWindowThick/2., 0., 2*M_PI);
  G4Tubs* winOuterTube = new G4Tubs("winOuterTube", 0., dMWPCExitR, dWindowThick/2., 0., 2*M_PI);

  // some visualization parameters
  G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
  G4VisAttributes* visKevlar = new G4VisAttributes(G4Colour(1,1.0,0,1.0));

  // logical volumes for the kevlar strips and Mylar windows
  kevContainer_log = new G4LogicalVolume(kevContainerBox, Vacuum, Append(side, "kevContainer_log_"));
  kevSeg_log = new G4LogicalVolume(kevSegBox, Vacuum, Append(side, "kevSeg_log_"));
  kevStrip_log = new G4LogicalVolume(kevStripBox, Kevlar, Append(side, "kevStrip_log_"));
  kevStrip_log -> SetVisAttributes(visKevlar);
  winIn_log = new G4LogicalVolume(winInnerTube, Mylar, Append(side, "winIn_log_"));
  winIn_log -> SetVisAttributes(visWindow);
  winOut_log = new G4LogicalVolume(winOuterTube, Mylar, Append(side, "winOut_log_"));
  winOut_log -> SetVisAttributes(visWindow);

  // placement of kevlar strings
  new G4PVPlacement(NULL, G4ThreeVector(0,0, kev_PosZ), kevContainer_log, Append(side, "kevContainer_phys_"),
                      container_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0,0), kevStrip_log, Append(side, "kevStrip_phys_"), kevSeg_log, false, 0);
  new G4PVReplica(Append(side, "kevlar_plane_"), kevSeg_log, kevContainer_log, kXAxis, NbKevWires, kevSpacing);

  // placement of Mylar windows
  new G4PVPlacement(NULL, G4ThreeVector(0,0, -dMWPCContainerHalf_Z + kevEffThick + dWindowThick/2.),
      winIn_log, Append(side, "winIn_phys_"), container_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, dMWPCContainerHalf_Z - dWindowThick/2.),
      winOut_log, Append(side, "winOut_phys_"), container_log, false, 0);

}

void WirechamberConstruction::GetFieldValue(const G4double Point[4], G4double* Bfield) const
{
  // set magnetic field
  if(fMyBField) fMyBField->GetFieldValue(Point,Bfield);
  else Bfield[0]=Bfield[1]=Bfield[2]=0;

  if(!dE0) { Bfield[3]=Bfield[4]=Bfield[5]=0; return; }

  // local position
  G4ThreeVector localPos = G4ThreeVector(Point[0],Point[1],Point[2])-vMyTranslation;
  if(fMyRotation) localPos = (*fMyRotation)(localPos);

  // electric field components
  G4ThreeVector E(0,0,0);
  double l = localPos[2];
  if(fabs(l) < dL)
  {
    double a = localPos[0]/dd;
    a = (a-floor(a)-0.5)*dd;
    if(a*a+l*l > dr*dr)
    {
      double denom = cosh(2*M_PI*l/dd)-cos(2*M_PI*a/dd);
      E[2] = dE0*sinh(2*M_PI*l/dd)/denom;
      E[0] = dE0*sin(2*M_PI*a/dd)/denom;
    }
  }

  // return to global coordinates
  if(fMyRotation) E = fMyRotation->inverse()(E);
  Bfield[3] = E[0];
  Bfield[4] = E[1];
  Bfield[5] = E[2];
}

void WirechamberConstruction::SetPotential(G4double Vanode)
{
  dE0 = M_PI*Vanode/dd/log(sinh(M_PI*dL/dd)/sinh(M_PI*dr/dd));
  G4cout << "Wirechamber voltage set to " << Vanode/volt <<" V => dE0 = " << dE0/(volt/cm) << " V/cm" << G4endl;

}

void WirechamberConstruction::ConstructField()
{
  G4cout << "Setting up wirechamber electromagnetic field...";
  // local field manager
  G4FieldManager* localFieldMgr = new G4FieldManager();
  localFieldMgr -> SetDetectorField(fMyBField);

  // equation of motion, stepper for field
  G4EqMagElectricField* pEquation = new G4EqMagElectricField(fMyBField);
  G4ClassicalRK4* pStepper = new G4ClassicalRK4(pEquation,8);
  G4MagInt_Driver* pIntgrDriver = new G4MagInt_Driver(0.01*um,pStepper,pStepper->GetNumberOfVariables());
  G4ChordFinder* pChordFinder = new G4ChordFinder(pIntgrDriver);
  localFieldMgr -> SetChordFinder(pChordFinder);

  // accuracy settings
  localFieldMgr -> GetChordFinder()->SetDeltaChord(10*um);
  localFieldMgr -> SetMinimumEpsilonStep(1e-6);
  localFieldMgr -> SetMaximumEpsilonStep(1e-5);
  localFieldMgr -> SetDeltaOneStep(0.1*um);

  // apply field manager to wirechamber and all daughter volumes
  container_log -> SetFieldManager(localFieldMgr,true);
  G4cout << " Done." << G4endl;
}
