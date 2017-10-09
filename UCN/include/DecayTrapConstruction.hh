#ifndef DecayTrapConstruction_HH
#define DecayTrapConstruction_HH

#include "DetectorConstructionUtils.hh"

#include <G4LogicalVolume.hh>

class DecayTrapConstruction: public MaterialUser
{
public:
  DecayTrapConstruction();	// constructor

//  G4double dWindowThick;
//  G4double dCoatingThick;

  G4double dWindowThick[2];	// these are M.Brown's additional class properties
  G4double dCoatingThick[2];	// needed like this for 2012-2013 geometry.
  G4double dThicknessTrapWindow[2];
  G4double dBeWinPosZ[2];
  G4double dMylarWinPosZ[2];
  G4double dTrapWinPosZ[2];


  G4double dInnerRadiusOfTrap;	// the class properties below here are from the old iteration of Xuan's sim.
  G4double dTubeWallThickness;
  G4double dInnerRadiusOfCollimator;

  G4double dCollimatorThick;	// Michael Brown's addition. Has been defaulted to M.M's 0.8*inch.
				// But in the 2011/2012 geometry gets set to 0.7*inch
  G4double dCollimatorLength;

  G4Material* mTubeMat;
  G4Material* mCollimatorMat;
  G4Material* mDecayTrapWindowMat;
  G4Material* mDecayTrapCoatingMat;

  G4Tubs* tTrapWindowTube[2];	// these G4Tubs* are added by M.Brown.
  G4Tubs* tMylarTube[2];
  G4Tubs* tBeTube[2];

  G4LogicalVolume* decayTrapTube_log; 	///< decay trap tube
  G4LogicalVolume* decayTrapWin_log[2]; 	///< trap window volume
  G4LogicalVolume* mylarWin_log[2]; 	///< mylar layer of window
  G4LogicalVolume* beWin_log[2]; 	///< berillium layer of window
  G4LogicalVolume* trapMonitor_log[2];	///< extra event monitoring region
  G4LogicalVolume* collimator_log[2]; 	///< collimator
  G4LogicalVolume* collimatorBack_log[2]; ///< bracket behind collimator

  //WiggleSheet wigglefoils[2]; ///< optional replacement crinkly endcap foils

  // construct is passed-to world volume, since we are not using an overall trap "container" vol.
  void Build(G4LogicalVolume* world, float crinkleAngle = 0);

};

#endif
