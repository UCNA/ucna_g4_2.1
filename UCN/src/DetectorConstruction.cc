#include "DetectorConstruction.hh"
#include "GlobalField.hh"
#include "MWPCField.hh"

#include <G4UserLimits.hh>		// stole from Michael Mendenhall's code.

#include <Randomize.hh>			// Stolen from Analysis Manager
#include <G4ios.hh>			// Pretty sure needed for TrackerSD
#include <G4Run.hh>			// Leave them here since we use registerSD in DetectorConstruction
#include <G4Event.hh>			// And the registerSD is totally not working without it
#include <G4Track.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4IonTable.hh>
#include <G4SDManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4SDManager.hh>
#include <G4EventManager.hh>

#include <G4EqMagElectricField.hh>
#include <G4ClassicalRK4.hh>
#include <G4MagneticField.hh>		// Bottom half of detector construction
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4PropagatorInField.hh>
#include <G4TransportationManager.hh>
#include <G4UserLimits.hh>
#include <G4PVParameterised.hh>

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixMixedStepper.hh"


DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScintStepLimit(1.0*mm),	// note: fScintStepLimit initialized here
  fStorageIndex(0),	// this variable loops over our TrackerHit names storage index
  fSourceFoilThick(9.4*um),
  fCrinkleAngle(0*rad)
{

}


DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  SetVacuumPressure(0);	// this is the set vacuum pressure that was warned about in DefineMaterials()

  // user step limits
  G4UserLimits* UserCoarseLimits = new G4UserLimits();
  UserCoarseLimits->SetMaxAllowedStep(10*m);
  G4UserLimits* UserGasLimits = new G4UserLimits();
  UserGasLimits->SetMaxAllowedStep(1*cm);
  G4UserLimits* UserSolidLimits = new G4UserLimits();
  UserSolidLimits->SetMaxAllowedStep(fScintStepLimit);	// default value from Messenger class.

  // Experimental Hall. World volume.
  G4double expHall_x = 2.0*m;
  G4double expHall_y = 2.0*m;
  G4double expHall_z = 8.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box", expHall_x/2, expHall_y/2, expHall_z/2);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, Vacuum, "World_log");
  experimentalHall_log -> SetVisAttributes(G4VisAttributes::Invisible);
  experimentalHall_log -> SetUserLimits(UserCoarseLimits);
  experimentalHall_phys = new G4PVPlacement(NULL, G4ThreeVector(), "World_phys", experimentalHall_log, 0, false, 0);

  //----- Source holder object. Used if it is a calibration source.
  // fSourceFoilThick set to 9.4*um. In the geantgen_#.mac it is 9.5*um.
  // So ignore the Messenger default which is 7.2*um
  Source.dSourceWindowThickness = fSourceFoilThick/2.;	// should be same as source_windowThick = 4.7*um
  Source.Build();
  // And we don't place the source holder object until later.
//  source_phys = new G4PVPlacement(NULL, source_holderPos, source_container_log, "source_container_phys",
//                               experimentalHall_log, false, 0, true);


  //----- Decay Trap object (length 3m, main tube)
  Trap.dWindowThick = 0.50*um;	// M.Brown sets these before we enter geometry choice
  Trap.dCoatingThick = 0.150*um;

  // 2011/2012 geometry settings are:
  Trap.dWindowThick = 0.500*um;
  Trap.mDecayTrapWindowMat = Mylar;
  Trap.dInnerRadiusOfCollimator = 2.3*inch;
  Trap.dCollimatorThick = 0.7*inch;
  // Stuff pertaining to wirechamber volume cathode/anode radius
//  G4double wireVol_anodeRadius = 5*um;
//  G4double wireVol_cathodeRadius = 39.1*um;

  // make the DecayTrapConstruction object.
  Trap.Build(experimentalHall_log, fCrinkleAngle);

  //----- Attempt to place everything.
  G4RotationMatrix* EastSideRot = new G4RotationMatrix();
  EastSideRot -> rotateY(M_PI*rad);

/*    G4ThreeVector sideTransFrame = G4ThreeVector(0., 0., 2.2*m);
    G4ThreeVector sideTransScint = sideTransFrame + G4ThreeVector(0., 0., -Scint[0].GetScintFacePos());

    G4double frame_backWinFrameThick = 0.5*inch;  // originally placed further down but needed here
    G4double mwpc_containerHalf_Z = (Wirechamber[0].GetWidth())/2.;       // this should be same for 0 or 1
    G4double mwpc_PosZ = -mwpc_containerHalf_Z - frame_backWinFrameThick
                - (Scint[0].GetWidth()/2. + Scint[0].GetScintFacePos());
    G4ThreeVector sideTransMWPC = sideTransFrame + G4ThreeVector(0, 0, mwpc_PosZ);
*/
  for(int i = 0; i <= 1; i++)
  {
    G4RotationMatrix* sideRot = new G4RotationMatrix();
    if(i == 0) sideRot -> rotateY(M_PI*rad);

    Wirechamber[i].ActiveRegion.dAnodeRadius = 5*um;    // Michael Brown's pre-sets for 2011/2012
    Wirechamber[i].ActiveRegion.dCathodeRadius = 39.1*um;

    Scint[i].Build(i);
    Wirechamber[i].Build(i);
    Frame[i].Build(i, Wirechamber[i].mwpcOverall_shape, Scint[i].scintOverall_shape,
                Wirechamber[i].GetWidth(), Scint[i].GetWidth(), Scint[i].GetScintFacePos(),
                Wirechamber[i].dMWPCExitR, Wirechamber[i].dMWPCEntranceR);

    G4ThreeVector sideTransFrame = G4ThreeVector(0., 0., 2.2*m);
    G4ThreeVector sideTransScint = sideTransFrame + G4ThreeVector(0., 0., -Scint[0].GetScintFacePos());

    G4double frame_backWinFrameThick = 0.5*inch;  // originally placed further down but needed here
    G4double mwpc_containerHalf_Z = (Wirechamber[0].GetWidth())/2.;       // this should be same for 0 or 1
    G4double mwpc_PosZ = -mwpc_containerHalf_Z - frame_backWinFrameThick
                - (Scint[0].GetWidth()/2. + Scint[0].GetScintFacePos());

    G4ThreeVector sideTransMWPC = sideTransFrame + G4ThreeVector(0, 0, mwpc_PosZ);

    Wirechamber[i].vMyTranslation += Sign(i)*sideTransMWPC;

    scint_phys[i] = new G4PVPlacement(sideRot, Sign(i)*sideTransScint, Scint[i].container_log,
				Append(i, "scintContainer_"), experimentalHall_log, false, 0, true);
    mwpc_phys[i] = new G4PVPlacement(sideRot, Sign(i)*sideTransMWPC, Wirechamber[i].container_log,
				Append(i, "mwpcContainer_"), experimentalHall_log, false, 0, true);
    frame_phys[i] = new G4PVPlacement(sideRot, Sign(i)*sideTransFrame, Frame[i].container_log,
				Append(i, "framePackage_"), experimentalHall_log, false, 0, true);


  }
/*
    scint_phys[0] = new G4PVPlacement(EastSideRot, (-1)*sideTransScint, Scint[0].container_log,
                                Append(0, "scintContainer_"), experimentalHall_log, false, 0, true);
    mwpc_phys[0] = new G4PVPlacement(EastSideRot, (-1)*sideTransMWPC, Wirechamber[0].container_log,
                                Append(0, "mwpcContainer_"), experimentalHall_log, false, 0, true);
    frame_phys[0] = new G4PVPlacement(EastSideRot, (-1)*sideTransFrame, Frame[0].container_log,
                                Append(0, "framePackage_"), experimentalHall_log, false, 0, true);
    scint_phys[1] = new G4PVPlacement(NULL, sideTransScint, Scint[1].container_log,
                                Append(1, "scintContainer_"), experimentalHall_log, false, 0, true);
    mwpc_phys[1] = new G4PVPlacement(NULL, sideTransMWPC, Wirechamber[1].container_log,
                                Append(1, "mwpcContainer_"), experimentalHall_log, false, 0, true);
    frame_phys[1] = new G4PVPlacement(NULL, sideTransFrame, Frame[1].container_log,
                                Append(1, "framePackage_"), experimentalHall_log, false, 0, true);
*/

  //----- Scintillator construction. Used as Sensitive Volume
//  G4ThreeVector sideTransScintEast = G4ThreeVector(0., 0., (-1)*(2.2*m - scint_face_PosZ));
//  G4ThreeVector sideTransScintWest = G4ThreeVector(0., 0., (2.2*m - scint_face_PosZ));
/*  G4RotationMatrix* EastSideRot = new G4RotationMatrix();
  EastSideRot -> rotateY(M_PI*rad);

  for(int i = 0; i <= 1; i++)
  {
    Scint[i].Build(i);

  }
  scint_phys[0] = new G4PVPlacement(EastSideRot, G4ThreeVector(0,0, Sign(0)*(2.2*m - Scint[0].GetScintFacePos())),
                                Scint[0].container_log, "scintContainer_EAST", experimentalHall_log, false, 0, true);
  scint_phys[1] = new G4PVPlacement(NULL, G4ThreeVector(0,0, Sign(1)*(2.2*m - Scint[1].GetScintFacePos())),
                                Scint[1].container_log, "scintContainer_WEST", experimentalHall_log, false, 0, true);

  //----- Wirechamber construction. Interior ActiveRegion used as sensitive volume.
  for(int i = 0; i <= 1; i++)
  {
    Wirechamber[i].ActiveRegion.dAnodeRadius = 5*um;	// Michael Brown's pre-sets for 2011/2012
    Wirechamber[i].ActiveRegion.dCathodeRadius = 39.1*um;

    Wirechamber[i].Build(i);
  }

//  G4double mwpc_exitRadius = Wirechamber[0].dMWPCExitR;
//  G4double mwpc_entranceRadius = Wirechamber[0].dMWPCEntranceR;
  G4double frame_backWinFrameThick = 0.5*inch;  // originally placed further down but needed here
  G4double mwpc_containerHalf_Z = (Wirechamber[0].GetWidth())/2.;       // this should be same for 0 or 1
  G4double mwpc_PosZ = -mwpc_containerHalf_Z - frame_backWinFrameThick
                - (Scint[0].GetWidth()/2. + Scint[0].GetScintFacePos());

  G4ThreeVector sideTransMWPCEast = G4ThreeVector(0,0, (-1)*(2.2*m + mwpc_PosZ));
  G4ThreeVector sideTransMWPCWest = G4ThreeVector(0, 0, (2.2*m + mwpc_PosZ));

  mwpc_phys[0] = new G4PVPlacement(EastSideRot, sideTransMWPCEast, Wirechamber[0].container_log,
				"mwpcContainer_phys_EAST", experimentalHall_log, false, 0, true);
  mwpc_phys[1] = new G4PVPlacement(NULL, sideTransMWPCWest, Wirechamber[1].container_log,
				"mwpcContainer_phys_WEST", experimentalHall_log, false, 0, true);

  //----- Finish up detector construction. Need to place frame_container_log in experimentalHall
  G4ThreeVector frameTransEast = G4ThreeVector(0., 0., (-1)*(2.2*m));	// note: scint face position is 0 in local coord.
									// Also there's no offset. So it's just -2.2m
  G4ThreeVector frameTransWest = G4ThreeVector(0., 0., 2.2*m);

  for(int i = 0; i <= 1; i++)
  {
    Frame[i].Build(i, Wirechamber[i].mwpcOverall_shape, Scint[i].scintOverall_shape,
		Wirechamber[i].GetWidth(), Scint[i].GetWidth(), Scint[i].GetScintFacePos(),
		Wirechamber[i].dMWPCExitR, Wirechamber[i].dMWPCEntranceR);
  }
  frame_phys[0] = new G4PVPlacement(EastSideRot, frameTransEast, Frame[0].container_log, "Frame_Package_EAST",
					experimentalHall_log, false, 0, true);
  frame_phys[1] = new G4PVPlacement(NULL, frameTransWest, Frame[1].container_log, "Frame_Package_WEST",
					experimentalHall_log, false, 0, true);
*/

  for(int i = 0; i <= 1; i++)			// set user limits in specific volumes
  {
    Trap.decayTrapWin_log[i] -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].container_log -> SetUserLimits(UserGasLimits);
    Wirechamber[i].winIn_log -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].winOut_log -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].kevStrip_log -> SetUserLimits(UserSolidLimits);
  }

  // register logical volumes as sensitive detectors. Used for all info tracking during sim
  for(int i = 0; i <= 1; i++)
  {
    SD_scint_scintillator[i] = RegisterSD(Append(i, "SD_scint_"), Append(i, "HC_scint_"));
    Scint[i].scintillator_log -> SetSensitiveDetector(SD_scint_scintillator[i]);
//    scint_scintillator_log[i] -> SetSensitiveDetector(SD_scint_scintillator[i]);

/*    SD_scint_deadScint[i] = RegisterSD(Append(i, "SD_deadScint_"));
    scint_deadLayer_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);
    scint_container_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);
    frame_mwpcExitGasN2_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);	// include N2 vol here
    scint_lightGuide_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);	// and also light guides

    SD_scint_backing[i] = RegisterSD(Append(i, "SD_scint_backing_"));
    scint_backing_log[i] -> SetSensitiveDetector(SD_scint_backing[i]);

    SD_mwpc_winIn[i] = RegisterSD(Append(i, "SD_mwpc_winIn_"));
    mwpc_winIn_log[i] -> SetSensitiveDetector(SD_mwpc_winIn[i]);

    SD_mwpc_winOut[i] = RegisterSD(Append(i, "SD_mwpc_winOut_"));
    mwpc_winOut_log[i] -> SetSensitiveDetector(SD_mwpc_winOut[i]);

    SD_decayTrap_windows[i] = RegisterSD(Append(i, "SD_decayTrap_windows_"));
    decayTrap_mylarWindow_log[i] -> SetSensitiveDetector(SD_decayTrap_windows[i]);
    decayTrap_beWindow_log[i] -> SetSensitiveDetector(SD_decayTrap_windows[i]);
*/
    SD_wireVol[i] = RegisterSD(Append(i, "SD_wireVol_"), Append(i, "HC_wireVol_"));
    Wirechamber[i].ActiveRegion.gas_log -> SetSensitiveDetector(SD_wireVol[i]);
    Wirechamber[i].ActiveRegion.anodeSeg_log -> SetSensitiveDetector(SD_wireVol[i]);
    Wirechamber[i].ActiveRegion.cathSeg_log -> SetSensitiveDetector(SD_wireVol[i]);
//    wireVol_gas_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
//    wireVol_anodeSeg_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
//    wireVol_cathSeg_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
/*
    SD_wireVol_planes[i] = RegisterSD(Append(i, "SD_wireVol_planes_"));
    wireVol_cathodeWire_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);
    wireVol_cathPlate_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);
    wireVol_anodeWire_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);

    SD_mwpc_container[i] = RegisterSD(Append(i, "SD_mwpc_container_"));	// equivalently, dead region in mwpc
    mwpc_container_log[i] -> SetSensitiveDetector(SD_mwpc_container[i]);

    SD_mwpc_kevStrip[i] = RegisterSD(Append(i, "SD_mwpc_kevStrip_"));
    mwpc_kevStrip_log[i] -> SetSensitiveDetector(SD_mwpc_kevStrip[i]);
*/
  }

/*  SD_source = RegisterSD("SD_source");
  source_window_log -> SetSensitiveDetector(SD_source);
  for(int i = 0; i <= 1; i++)
    source_coating_log[i] -> SetSensitiveDetector(SD_source);

  for(int i = 0; i <= 1; i++)
  {
    SD_decayTrap_innerMonitors[i] = RegisterSD(Append(i, "SD_decayTrap_innerMonitors_"));
    decayTrap_innerMonitors_log[i] -> SetSensitiveDetector(SD_decayTrap_innerMonitors[i]);
  }

  // exp hall vacuum, decay tube, other inert parts. Stores all energy "lost".
  SD_world = RegisterSD("SD_world");
  experimentalHall_log -> SetSensitiveDetector(SD_world);
  decayTrap_tube_log -> SetSensitiveDetector(SD_world);
  for(int i = 0; i <= 1; i++)
  {
    frame_mwpcEntrance_log[i] -> SetSensitiveDetector(SD_world);
    frame_mwpcExit_log[i] -> SetSensitiveDetector(SD_world);	// confusing, since made of Al. But trust M.M.
    frame_container_log[i] -> SetSensitiveDetector(SD_world);
    decayTrap_collimator_log[i] -> SetSensitiveDetector(SD_world);	// and this is polyethylene?
    decayTrap_collimatorBack_log[i] -> SetSensitiveDetector(SD_world);
  }
*/

  G4double mwpc_fieldE0 = 2700*volt;
  G4ThreeVector mwpc_activeRegionTrans(0, 0, (Wirechamber[0].dEntranceToCathodes - Wirechamber[0].dExitToCathodes)/2.);

  // Create everything needed for global and local EM fields
//  G4ThreeVector East_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCEast;
//  G4ThreeVector West_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCWest;

  G4ThreeVector East_EMFieldLocation = mwpc_activeRegionTrans + Wirechamber[0].vMyTranslation;
  G4ThreeVector West_EMFieldLocation = mwpc_activeRegionTrans + Wirechamber[1].vMyTranslation;


  ConstructGlobalField();			// make magnetic and EM fields.

  ConstructEastMWPCField(Wirechamber[0].ActiveRegion.dWireSpacing,
			Wirechamber[0].ActiveRegion.dPlaneSpacing,
			Wirechamber[0].ActiveRegion.dAnodeRadius,
			mwpc_fieldE0, EastSideRot, East_EMFieldLocation);
  ConstructWestMWPCField(Wirechamber[1].ActiveRegion.dWireSpacing,
                        Wirechamber[1].ActiveRegion.dPlaneSpacing,
                        Wirechamber[1].ActiveRegion.dAnodeRadius,
                        mwpc_fieldE0, NULL, West_EMFieldLocation);
  return experimentalHall_phys;
}

TrackerSD* DetectorConstruction::RegisterSD(G4String sdName, G4String hcName)
{
  TrackerSD* sd = new TrackerSD(sdName, hcName);
  G4SDManager::GetSDMpointer() -> AddNewDetector(sd);

  fSDNamesArray[fStorageIndex] = sdName;
  fHCNamesArray[fStorageIndex] = hcName;
  fStorageIndex++;

  return sd;
}

void DetectorConstruction::ConstructGlobalField()
{
  G4cout << "Setting up global magnetic field. Call to global field object." << G4endl;

  GlobalField* magField = new GlobalField();
  G4FieldManager* globalFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  globalFieldManager -> SetDetectorField(magField);
  globalFieldManager -> CreateChordFinder(magField);

  G4MagIntegratorStepper* pStepper;
  G4Mag_UsualEqRhs* equationOfMotion = new G4Mag_UsualEqRhs(magField);
  //pStepper = new G4ClassicalRK4 (fEquation); // general case for "smooth" EM fields
  //pStepper = new G4SimpleHeum( fEquation ); // for slightly less smooth EM fields
  //pStepper = new G4HelixHeum( fEquation ); // for "smooth" pure-B fields
  //pStepper = new G4HelixImplicitEuler( fEquation ); // for less smooth pure-B fields; appears ~50% faster than above
  //pStepper = new G4HelixSimpleRunge( fEquation ); // similar speed to above
  //pStepper = new G4HelixExplicitEuler( fEquation ); // about twice as fast as above
  pStepper = new G4HelixMixedStepper(equationOfMotion,6); // avoids "Stepsize underflow in Stepper" errors
  globalFieldManager -> GetChordFinder() -> GetIntegrationDriver() -> RenewStepperAndAdjust(pStepper);

  globalFieldManager -> GetChordFinder() -> SetDeltaChord(100.0*um);
  globalFieldManager -> SetMinimumEpsilonStep(1e-6);
  globalFieldManager -> SetMaximumEpsilonStep(1e-5);
  globalFieldManager -> SetDeltaOneStep(0.1*um);
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);

  return;
}

void DetectorConstruction::ConstructEastMWPCField(G4double a, G4double b, G4double c, G4double d, G4RotationMatrix* e, G4ThreeVector f)
{
  G4cout << "Setting up East wirechamber electromagnetic field." << G4endl;
  MWPCField* eastLocalField = new MWPCField();
  eastLocalField -> SetActiveReg_d(a);
  eastLocalField -> SetActiveReg_L(b);
  eastLocalField -> SetActiveReg_r(c);
  eastLocalField -> SetSideRot(e);
  eastLocalField -> SetSideTrans(f);
  eastLocalField -> SetPotential(d);

  G4FieldManager* eastLocalFieldManager = new G4FieldManager();
  eastLocalFieldManager -> SetDetectorField(eastLocalField);

  G4EqMagElectricField* eastlocalEquation = new G4EqMagElectricField(eastLocalField);
  G4ClassicalRK4* eastlocalStepper = new G4ClassicalRK4(eastlocalEquation,8);
  G4MagInt_Driver* eastlocalIntgrDriver = new G4MagInt_Driver(0.01*um,eastlocalStepper,eastlocalStepper->GetNumberOfVariables());
  G4ChordFinder* eastlocalChordFinder = new G4ChordFinder(eastlocalIntgrDriver);
  eastLocalFieldManager -> SetChordFinder(eastlocalChordFinder);

  eastLocalFieldManager -> GetChordFinder() -> SetDeltaChord(10*um);
  eastLocalFieldManager -> SetMinimumEpsilonStep(1e-6);
  eastLocalFieldManager -> SetMaximumEpsilonStep(1e-5);
  eastLocalFieldManager -> SetDeltaOneStep(0.1*um);

  Wirechamber[0].container_log -> SetFieldManager(eastLocalFieldManager, true);
  return;
}

void DetectorConstruction::ConstructWestMWPCField(G4double a, G4double b, G4double c, G4double d, G4RotationMatrix* e, G4ThreeVector f)
{
  G4cout << "Setting up West wirechamber electromagnetic field." << G4endl;
  MWPCField* westLocalField = new MWPCField();
  westLocalField -> SetActiveReg_d(a);
  westLocalField -> SetActiveReg_L(b);
  westLocalField -> SetActiveReg_r(c);
  westLocalField -> SetSideRot(e);
  westLocalField -> SetSideTrans(f);
  westLocalField -> SetPotential(d);

  G4FieldManager* westLocalFieldManager = new G4FieldManager();
  westLocalFieldManager -> SetDetectorField(westLocalField);

  G4EqMagElectricField* westlocalEquation = new G4EqMagElectricField(westLocalField);
  G4ClassicalRK4* westlocalStepper = new G4ClassicalRK4(westlocalEquation,8);
  G4MagInt_Driver* westlocalIntgrDriver = new G4MagInt_Driver(0.01*um,westlocalStepper,westlocalStepper->GetNumberOfVariables());
  G4ChordFinder* westlocalChordFinder = new G4ChordFinder(westlocalIntgrDriver);
  westLocalFieldManager -> SetChordFinder(westlocalChordFinder);

  westLocalFieldManager -> GetChordFinder() -> SetDeltaChord(10*um);
  westLocalFieldManager -> SetMinimumEpsilonStep(1e-6);
  westLocalFieldManager -> SetMaximumEpsilonStep(1e-5);
  westLocalFieldManager -> SetDeltaOneStep(0.1*um);

  Wirechamber[1].container_log -> SetFieldManager(westLocalFieldManager, true);
  return;
}
