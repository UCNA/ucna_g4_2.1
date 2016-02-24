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

  //----- Scintillator construction. Used as Sensitive Volume
//  G4ThreeVector sideTransScintEast = G4ThreeVector(0., 0., (-1)*(2.2*m - scint_face_PosZ));
//  G4ThreeVector sideTransScintWest = G4ThreeVector(0., 0., (2.2*m - scint_face_PosZ));
  G4RotationMatrix* EastSideRot = new G4RotationMatrix();
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
  G4double mwpc_fieldE0 = 2700*volt;
  G4ThreeVector mwpc_activeRegionTrans(0, 0, (Wirechamber[0].dEntranceToCathodes - Wirechamber[0].dExitToCathodes)/2.);

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




  //----- Begin DetectorPackageConstruction. This is the frame that holds the scintillator and MWPC.
/*  G4double frame_packageRadius = 6.0*inch;
  G4double frame_mwpcEntranceThick = 0.375*inch;
  G4double frame_mwpcEntranceRadius = 3.0*inch;	// not the same value as mwpc_entranceRadius
  G4double frame_mwpcEntranceDepth = 5.0*inch;
  G4double frame_frontWinFrameThick = 1.0*inch;
//  G4double frame_backWinFrameThick = 0.5*inch;	// moved a few lines up since needed for mwpc placement

  // create the shapes that will be used for geometric objects below
  // aluminum entrance collimator to detector package
  G4double frame_entranceSectionLength = frame_mwpcEntranceDepth + frame_frontWinFrameThick;
  G4Tubs* frame_mwpcEntranceTube = new G4Tubs("mwpc_entrance_tube", 0., frame_packageRadius, 0.5*frame_entranceSectionLength, 0., 2*M_PI);
  G4Tubs* frame_entranceFrontTube = new G4Tubs("entrance_front_tube", frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					frame_packageRadius, 0.5*frame_mwpcEntranceThick, 0., 2*M_PI);
  G4Tubs* frame_entranceMidTube = new G4Tubs("entrance_mid_tube", frame_mwpcEntranceRadius, frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					0.5*frame_mwpcEntranceDepth, 0., 2*M_PI);
  G4Tubs* frame_entranceBackTube = new G4Tubs("entrance_back_tube", mwpc_entranceRadius, frame_packageRadius,
					0.5*frame_frontWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));

  // create overall detector package frame tube.
//  G4double frame_detFrameHalf_Z = frame_mwpcEntranceDepth + 2*mwpc_containerHalf_Z + 1.0*inch;
  G4double frame_detFrameHalf_Z = frame_mwpcEntranceDepth + Wirechamber[0].GetWidth() + 1.0*inch;

  G4Tubs* frame_framePackageTube = new G4Tubs("detPackage_tube_EAST", 0, frame_packageRadius, frame_detFrameHalf_Z, 0, 2*M_PI);

  // Overall container layer for the scintillator
  G4Tubs* scint_N2VolTube = new G4Tubs("N2_vol_tube", 0., Scint[0].dBackingRadius, Scint[0].GetWidth()/2., 0., 2*M_PI);

  // subtract off the scintillator container volume. Not rotated (since frame will be made and then rotated).
  // But displaced by -scint_face_PosZ relative to the local coordinates of the detector package frame.
  G4SubtractionSolid* frame_frameMinusScint = new G4SubtractionSolid("DPC_frame_minus_scint_container_log", frame_framePackageTube,
					scint_N2VolTube, NULL, G4ThreeVector(0., 0., -Scint[0].GetScintFacePos()));
  // subtract off the mwpc container volume. Not rotated (since frame will be made then rotated).
  G4SubtractionSolid* frame_containerShape = new G4SubtractionSolid("frame_container_minus_Scint_MWPC", frame_frameMinusScint,
					mwpc_containerBox, NULL, G4ThreeVector(0., 0., mwpc_PosZ));

  // "place componenets relative to scintillator face at 0" is what old code says.
  // Already successfully placed the scint and the mwpc in experimentalHall_log. This code is remaining geom that wasn't done
  G4double frame_entrance_PosZ = mwpc_PosZ - (2*mwpc_containerHalf_Z + frame_entranceSectionLength)/2.;

  // aluminum exit window and N2 volume at back of gas box
  G4Tubs* frame_mwpcExitTube = new G4Tubs("mwpc_exit_tube", mwpc_exitRadius, frame_packageRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));

  G4double frame_exitWin_PosZ = mwpc_PosZ + (2*mwpc_containerHalf_Z + frame_backWinFrameThick)/2.;
  G4Tubs* frame_mwpcExitGasN2Tube = new G4Tubs("mwpc_exit_N2_tube", 0, mwpc_exitRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);

  // material behind the detector. Misc stuff that can cause back scattering events.
  G4double frame_backStuffThick = 1.0*inch;
  G4Tubs* frame_backStuffTube = new G4Tubs("backstuff_tube_EAST", 0, 0.5*frame_packageRadius, frame_backStuffThick, 0., 2*M_PI);

  // create logical volumes from all the shapes defined above. As always, array of 2 for East/West (0/1)
  for(int i = 0; i <= 1; i++)
  {
    frame_mwpcEntrance_log[i] = new G4LogicalVolume(frame_mwpcEntranceTube, Vacuum, Append(i, "mwpc_entrance_log_"));
    frame_entranceFront_log[i] = new G4LogicalVolume(frame_entranceFrontTube, Al, Append(i, "entrance_front_log_"));
    frame_entranceMid_log[i] = new G4LogicalVolume(frame_entranceMidTube, Al, Append(i, "entrance_mid_log_"));
    frame_entranceBack_log[i] = new G4LogicalVolume(frame_entranceBackTube, Al, Append(i, "entrance_back_log_"));
    frame_mwpcEntrance_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);
    frame_entranceFront_log[i] -> SetVisAttributes(visMWPCEntrance);
    frame_entranceMid_log[i] -> SetVisAttributes(visMWPCEntrance);
    frame_entranceBack_log[i] -> SetVisAttributes(visMWPCEntrance);

    frame_container_log[i] = new G4LogicalVolume(frame_containerShape, Vacuum, Append(i, "frame_container_log_"));
    frame_container_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);
//    frame_container_log[i] -> SetVisAttributes(new G4VisAttributes(G4Color(1, 0, 1, 1)));

    frame_mwpcExit_log[i] = new G4LogicalVolume(frame_mwpcExitTube, Al, Append(i, "mwpc_exit_log_"));
    frame_mwpcExit_log[i] -> SetVisAttributes(visMWPCExit);
    frame_mwpcExitGasN2_log[i] = new G4LogicalVolume(frame_mwpcExitGasN2Tube, WCNitrogen, Append(i, "mwpc_exit_N2_log_"));
    frame_mwpcExitGasN2_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);

    frame_backStuff_log[i] = new G4LogicalVolume(frame_backStuffTube, SS304, Append(i, "backStuff_log_"));
  }

  // place all the physical volumes (class members) in their respective mother volumes.
  for(int i = 0; i <= 1; i++)
  {  // These three are placed in the logical vol. mwpcEntrance.
    frame_entranceFront_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*(frame_entranceSectionLength - frame_mwpcEntranceThick)),
                                  frame_entranceFront_log[i], Append(i, "entrance_front_phys_"), frame_mwpcEntrance_log[i], false, 0);
    frame_entranceMid_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*frame_frontWinFrameThick),
                                  frame_entranceMid_log[i], Append(i, "entrance_mid_phys_"), frame_mwpcEntrance_log[i], false, 0);
    frame_entranceBack_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, 0.5*(frame_entranceSectionLength - frame_frontWinFrameThick)),
                                  frame_entranceBack_log[i], Append(i, "entrance_back_phys_"), frame_mwpcEntrance_log[i], false, 0);

    // All this stuff makes up the frame "container_log" i.e. the overall frame logical.
    frame_mwpcEntrance_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0., 0., frame_entrance_PosZ),
                                  frame_mwpcEntrance_log[i], Append(i, "frame_mwpc_entrance_"), frame_container_log[i], false, 0);
    frame_mwpcExit_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExit_log[i],
                                  Append(i, "mwpc_exit_"), frame_container_log[i], false, 0);
    frame_mwpcExitGasN2_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExitGasN2_log[i],
                                  Append(i, "mwpc_exit_N2_phys_"), frame_container_log[i], false, 0);
    frame_backStuff_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0, 0, frame_detFrameHalf_Z - 0.5*frame_backStuffThick),
                                  frame_backStuff_log[i], Append(i, "backStuff_phys_"), frame_container_log[i], false, 0);
  }
*/



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



//  frame_container_phys[0] = new G4PVPlacement(EastSideRot, frameTransEast, frame_container_log[0],
//				"Detector_Package_Frame_EAST", experimentalHall_log, false, 0, true);
//  frame_container_phys[1] = new G4PVPlacement(NULL, frameTransWest, frame_container_log[1],
//				"Detector_Package_Frame_WEST", experimentalHall_log, false, 0, true);

  for(int i = 0; i <= 1; i++)			// set user limits in specific volumes
  {
    Trap.decayTrapWin_log[i] -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].container_log -> SetUserLimits(UserGasLimits);
    Wirechamber[i].winIn_log -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].winOut_log -> SetUserLimits(UserSolidLimits);
    Wirechamber[i].kevStrip_log -> SetUserLimits(UserSolidLimits);
/*    mwpc_container_log[i] -> SetUserLimits(UserGasLimits);
    mwpc_winIn_log[i] -> SetUserLimits(UserSolidLimits);
    mwpc_winOut_log[i] -> SetUserLimits(UserSolidLimits);
    mwpc_kevStrip_log[i] -> SetUserLimits(UserSolidLimits); */
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

  // Create everything needed for global and local EM fields
  G4ThreeVector East_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCEast;
  G4ThreeVector West_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCWest;

  ConstructGlobalField();			// make magnetic and EM fields.

  ConstructEastMWPCField(Wirechamber[0].ActiveRegion.dWireSpacing,
			Wirechamber[0].ActiveRegion.dPlaneSpacing,
			Wirechamber[0].ActiveRegion.dAnodeRadius,
			mwpc_fieldE0, EastSideRot, East_EMFieldLocation);
  ConstructWestMWPCField(Wirechamber[1].ActiveRegion.dWireSpacing,
                        Wirechamber[1].ActiveRegion.dPlaneSpacing,
                        Wirechamber[1].ActiveRegion.dAnodeRadius,
                        mwpc_fieldE0, NULL, West_EMFieldLocation);


/*  ConstructEastMWPCField(wireVol_wireSpacing, wireVol_planeSpacing, wireVol_anodeRadius,
			mwpc_fieldE0, EastSideRot, East_EMFieldLocation);
  ConstructWestMWPCField(wireVol_wireSpacing, wireVol_planeSpacing, wireVol_anodeRadius,
			mwpc_fieldE0, NULL, West_EMFieldLocation); */
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
//  mwpc_container_log[0] -> SetFieldManager(eastLocalFieldManager, true);
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
//  mwpc_container_log[1] -> SetFieldManager(westLocalFieldManager, true);
  return;
}
