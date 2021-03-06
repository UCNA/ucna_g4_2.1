#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"	// original example used these 3
#include "G4ParticleGun.hh"
#include "globals.hh"

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleGun.hh>
#include <G4Event.hh>
#include <G4VUserEventInformation.hh>

#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWithAString.hh>


class G4ParticleGun;
class G4Event;
class TTree;
class TFile;

// user event information for recording primary event weighting
class PrimEvtWeighting : public G4VUserEventInformation
{
  public:
    PrimEvtWeighting(double W): w(W) {}	// constructor
    void Print() const { G4cout << "Primary weighting: " << w << G4endl; }

    double w;				// event primary weight
};

using namespace std;

struct event
{
  G4int event_gen_id;
  G4double event_energy;	// turns into keV
  G4int event_speciesFlag;      // 11 means electron, 22 is gamma
  G4double event_xMo;		// momentum is unitless vector
  G4double event_yMo;
  G4double event_zMo;
  G4double event_xPos;		// turns into m
  G4double event_yPos;
  G4double event_zPos;
  G4double event_time;		// turns into s but always 0. I think supposed to be ns
  G4double event_weight;	// unitless
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, G4UImessenger
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // UI interface
    virtual void SetNewValue(G4UIcommand* command, G4String newValue);

  private:
    G4ParticleGun*  fParticleGun; 	// pointer a to G4 gun class
    DetectorConstruction* fMyDetector;	// pointer to the detector geometry class

    bool bIsLoaded;			// check if we've already loaded the input ptcl file

    vector <G4double> dCoincidenceEnergies;
    vector <G4int> iCoincidencePtcl;		// flag to see if we want to produce a coincidence ptcl
					//	-1 = default value in code
					//	+1 = Auger electron, ...
    G4int iNbCoincidence;		// number of coincidence particles we'll produce
    G4int iMaxNbCoin;
    G4bool bCoincidenceWasFired;	// keeps track of the ordering of Set_113Sn vs SavePrimPtclInfo

    G4double fSourceRadius; // spread radius for source droplets
    G4ThreeVector fPosOffset;		// base positioning offset
    event fEvtsArray[2000000];		// size has to be number of lines in input file

    void UseExternalKinFile(int nID);
    void LoadFile(G4String fileName);
    void ConvertTreeToTxt(G4String treeName, G4String txtFileName);

    void DiskRandom(G4double radius, G4double& x, G4double& y);
    void DisplayGunStatus();
    void SavePrimPtclInfo(int index);
    void SavePrimPtclInfo2(int eventID);	// used for particles generated in the GEANT4 sim

    void Set_113SnSource();
    void Set_139CeSource();
    void Set_207BiSource();

    // UI command variables for messenger class input/output files
    G4UIdirectory* uiGenDir;	// UI directory for primaryGeneratorAction related commands

    G4UIcmdWithAString* uiInputFileCmd;	// which input file name to take
    G4String sInputFileName;

    G4UIcmdWithAString* uiOutputFileCmd;	// which output file name to take
    G4String sOutputFileName;

    G4UIcmdWithAString* uiParticleGenTypeCmd;
    G4String sPtclType;

};

#endif


