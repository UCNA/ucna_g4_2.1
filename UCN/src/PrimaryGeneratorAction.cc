#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
using   namespace       std;

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC)
: G4VUserPrimaryGeneratorAction(), G4UImessenger(),
  fParticleGun(0),
  fMyDetector(myDC),
  fSourceRadius(3.*mm),	// this is default set. If you aren't using DiskRandom, don't care.
  fPosOffset(G4ThreeVector(0,0,0)),	// base positioning offset. Is non-zero only if main Decay Trap is offset (it is not)
  bIsLoaded(false), bUseExternal(false),
  iCoincidencePtcl(-1), iNbCoincidence(0), bCoincidenceWasFired(false)
{
  //----- Below is messenger class

  uiGenDir = new G4UIdirectory("/particleGun/");
  uiGenDir -> SetGuidance("/primaryGeneratorAction control");

  uiInputFileCmd = new G4UIcmdWithAString("/particleGun/inputName", this);
  uiInputFileCmd -> SetGuidance("Set the input file name of the primaries.");
  uiInputFileCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
  sInputFileName = "none.txt";

  uiOutputFileCmd = new G4UIcmdWithAString("/particleGun/outputName", this);
  uiOutputFileCmd -> SetGuidance("Set the output file name of the primaries.");
  uiOutputFileCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
  sOutputFileName = "none.txt";

  //----- Above is messenger class

  G4int nPtcls = 1;
  fParticleGun = new G4ParticleGun(nPtcls);

  // At the end of constructor, GEANT4 default calls GeneratePrimaries method
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == uiInputFileCmd)
  {
    sInputFileName = G4String(newValue);
    G4cout << "Setting input primary particles file to " << sInputFileName << G4endl;
  }
  else if(command == uiOutputFileCmd)
  {
    sOutputFileName = G4String(newValue);
    G4cout << "Setting output primary particles information file to " << sOutputFileName << G4endl;
  }
  else
    G4cout << "COMMAND DOES NOT MATCH PRIMARY GENERATOR ACTION OPTIONS." << G4endl;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // use GEANT4 event id to track which part of fEvtsArray we will use for generated event
  int evtID = anEvent -> GetEventID();

  if(bUseExternal == true)
  {
    if(evtID == 0)
    {
      G4cout << "Using external kinematics file to generate initial particles... " << G4endl;
    }
    UseExternalKinFile(evtID);
    SavePrimPtclInfo(evtID);
  }
  else if(bUseExternal == false)
  {
    // note this will look different than loading the external kinematics file
    // because we need to be able to account for Auger and other coincidence particles
    if(evtID == 0)
    {
      G4cout << "Using GEANT4 particle generation to create initial particles..." << G4endl;
    }

    Set_113SnSource();	// processes all decays of 113Sn that we are interested in

    if(bCoincidenceWasFired == false)
    {
      SavePrimPtclInfo(evtID);
    }
    else if(bCoincidenceWasFired == true)
    {
      // save primary ptcl info as same event number so we can add up later
      SavePrimPtclInfo(evtID - iNbCoincidence);
      // since we've now accounted for one coincidence, lower the number of coincidences in the counter
      iNbCoincidence = iNbCoincidence - 1;
      if(iNbCoincidence < 0)
	G4cout << "RUN-TIME ERROR IN CODE. iNbCoincidence went below 0. Makes no sense." << G4endl;

      if(iNbCoincidence == 0)
        iCoincidencePtcl = -1;		// if we have no more coincidences, reset the ptcl flag to -1

    }

  }


  // Call to method to save primary particle initial info.
  // Need to be super careful here. A Priori, there's no reason that the momentum vector
  // is normalized when I'm reading it from final.
  // And momentum direction is a normalized vector once GEANT4 gets a hold of it.
  // So when we print out what was read in vs. what we get from particle gun may not be the same.
  // Will need to check this. For now it is the same within float to double rounding.
//  SavePrimPtclInfo(evtID);

  fParticleGun -> GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::UseExternalKinFile(int nID)
{
  if(bIsLoaded == false)
  {
    G4cout << "------> Fetching initial particles info from: " << sInputFileName << G4endl;
    LoadFile(sInputFileName);

    if(fMyDetector -> GetUseSourceHolder() == false)
    {
      G4cout << "\n Loading neutron beta decay events..." << G4endl;
    }
    else if(fMyDetector -> GetUseSourceHolder() == true)
    {
      G4cout << "\n Adjusting initial ptcl position for source holder use..." << G4endl;
      G4cout << " Randomly throwing events in source holder radius: " << fSourceRadius/mm << " mm." << G4endl;
    }

  }

  // use the particle species flag to set particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
  if(fEvtsArray[nID].event_speciesFlag == 11)
  {
    particle = particleTable->FindParticle("e-");
  }
  else if(fEvtsArray[nID].event_speciesFlag == 22)
  {
    particle = particleTable->FindParticle("gamma");
  }
  else
  {
    G4cout << "No matching particle species definition." << G4endl;
  }

  fParticleGun -> SetParticleDefinition(particle);

  fParticleGun -> SetParticleEnergy(fEvtsArray[nID].event_energy*keV);

  G4ThreeVector pos = fPosOffset;       // so far in my simulation this is 0
  if(fMyDetector -> GetUseSourceHolder() == false)
  {
    pos[0] += fEvtsArray[nID].event_xPos*m;
    pos[1] += fEvtsArray[nID].event_yPos*m;
    pos[2] += fEvtsArray[nID].event_zPos*m;
    fParticleGun -> SetParticlePosition(pos);
  }
  else if(fMyDetector -> GetUseSourceHolder() == true)
  {
    pos += fMyDetector ->  GetSourceHolderPos();
    if(fSourceRadius != 0)
    {
      G4double x,y;
      DiskRandom(fSourceRadius, x, y);
      pos[0] += x;
      pos[1] += y;

      pos[0] += fEvtsArray[nID].event_xPos*m;
      pos[1] += fEvtsArray[nID].event_yPos*m;
      pos[2] += fEvtsArray[nID].event_zPos*m;
    }
    fParticleGun -> SetParticlePosition(pos);
  }

  fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fEvtsArray[nID].event_xMo,
                                                fEvtsArray[nID].event_yMo,
                                                fEvtsArray[nID].event_zMo));

  fParticleGun -> SetParticleTime(fEvtsArray[nID].event_time*ns);

}

void PrimaryGeneratorAction::LoadFile(G4String fileName)
{
  event eRead;
  int i = 0;

  string buf;
  ifstream infile;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> eRead.event_gen_id
		>> eRead.event_speciesFlag
		>> eRead.event_energy
		>> eRead.event_xPos >> eRead.event_yPos >> eRead.event_zPos
		>> eRead.event_xMo >> eRead.event_yMo >> eRead.event_zMo
		>> eRead.event_time
		>> eRead.event_weight;

      fEvtsArray[i] = eRead;
      i++;
    }
  }

  // Set this class member to true so that we don't reload the input file every time
  bIsLoaded = true;
}

void PrimaryGeneratorAction::DiskRandom(G4double radius, G4double& x, G4double& y)
{
  while(true)
  {
    x = (2.0*G4UniformRand()-1.)*radius;	// random number between 0 and 2, subtract 1
    y = (2.0*G4UniformRand()-1.)*radius;	// goes randomly between -1 to +1
    if(x*x+y*y<=radius*radius) break;
  }
}

void PrimaryGeneratorAction::DisplayGunStatus()
{
  G4cout
  << fParticleGun->GetParticleDefinition()->GetParticleName()
  << " gun from " << fParticleGun->GetParticlePosition()/m
  << "m towards " << fParticleGun->GetParticleMomentumDirection()
  << " at " << fParticleGun->GetParticleTime()/ns
  << "ns : gun firing at energy " << fParticleGun->GetParticleEnergy()/keV
  << "keV" <<
  G4endl;
}

void PrimaryGeneratorAction::SavePrimPtclInfo(int index)
{
  ofstream outfile;     // output initial particle information into same file as final sim output
  outfile.open(sOutputFileName, ios::app);
  outfile << fEvtsArray[index].event_gen_id << "\t"
	<< fEvtsArray[index].event_speciesFlag << "\t"
	<< fEvtsArray[index].event_energy << "\t"
	<< fEvtsArray[index].event_xPos << "\t"
	<< fEvtsArray[index].event_yPos << "\t"
	<< fEvtsArray[index].event_zPos << "\t"
	<< fEvtsArray[index].event_xMo << "\t"	// be careful here with momentum
	<< fEvtsArray[index].event_yMo << "\t"	// see comment back in GeneratePrimaries
	<< fEvtsArray[index].event_zMo << "\t"
	<< fEvtsArray[index].event_time << "\t"
	<< fEvtsArray[index].event_weight << "\t";	// has to be a \t since getting appended in EventAction
  outfile.close();
}


void PrimaryGeneratorAction::Set_113SnSource()	// don't need additional arguments since we set the particle gun.
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  int r1;	// total percentage of stuff we're interested in for 113Sn
  r1 = rand() % 1029306 + 1;
  double percentage = r1/10000.;

  if(iNbCoincidence > 0)
  {
    if(iCoincidencePtcl == 1)		// create a K shell Auger electron for 113Sn
    {
      fParticleGun -> SetParticleEnergy(20.1*keV);
      particle = particleTable -> FindParticle(particleName="e-");
    }
    // technically, here you would put other options for other coincidences.
    // for 113Sn, there is no other coincidence that we care about.

    bCoincidenceWasFired = true;        // sets the flag stating that we created a ptcl
                                        // that is intended to be recorded in coincidence
  }
  else if(iNbCoincidence == 0)
  {
    bCoincidenceWasFired = false;	// sets the flag stating we are using a coincidence for this event to false
					// meaning this is a ptcl generated that is not supposed to be a coincidence.
    if((percentage >= 0) && (percentage <= 64.97))
    {	// fires a gamma in the decay from 1st excited to ground state
      fParticleGun -> SetParticleEnergy(391.698*keV);
      particle = particleTable->FindParticle(particleName="gamma");
    }
    else if((percentage > 64.97) && (percentage <= 100.7246))
    {	// produces a Conversion Electron from 1st excited state to GS

      int r2 = rand() % 1000000 + 1;
      double AugerPercent = r2/10000.;	// check whether the CE also comes with an Auger
      if((AugerPercent >= 0) && (AugerPercent <= 11.8056))
      {
        iCoincidencePtcl = 1;	// sets the flag to produce an Auger
        iNbCoincidence++;	// increments our coincidence counter so next event, we produce an Auger
      }

      // whether or not the Auger flag checks out, we produce all the CE down here.
      // These are, in order, the K, L, M, N, O shell Conversion Electron energies.
      if((percentage > 64.97) && (percentage <= 93.77))
      {
        fParticleGun -> SetParticleEnergy(363.758*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 93.77) && (percentage <= 99.37))
      {
        fParticleGun -> SetParticleEnergy(387.461*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 99.37) && (percentage <= 100.507))
      {
        fParticleGun -> SetParticleEnergy(390.872*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 100.507) && (percentage <= 100.712))
      {
        fParticleGun -> SetParticleEnergy(391.576*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 100.712) && (percentage <= 100.7246))
      {
        fParticleGun -> SetParticleEnergy(391.697*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
    }
    else if((percentage > 100.7246) && (percentage <= 102.8346))
    {	// fires a gamma ray in the decay of the second excited state to 1st excited
      fParticleGun -> SetParticleEnergy(255.134*keV);
      particle = particleTable -> FindParticle(particleName="gamma");
    }
    else if((percentage > 102.8346) && (percentage <= 102.9306))
    {
      // here is where we would typically throw another Auger check.
      // Since there is a non-zero probability that we'll get another, second, coincidence Auger.
      // But since it is such small prob (11% of 0.06%), we'll ignore it altogether.

	// fires CE from 2nd excited to 1st excited
	// in order, the CE are K, L, M, N shell CE being released. We assume neglible # of Auger's accompany them.
      if((percentage > 102.8346)  && (percentage <= 102.9166))
      {
        fParticleGun -> SetParticleEnergy(227.194*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 102.9166) && (percentage <= 102.928))
      {
        fParticleGun -> SetParticleEnergy(250.896*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 102.9280) && (percentage <= 102.9302))
      {
        fParticleGun -> SetParticleEnergy(254.308*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 102.9302) && (percentage <= 102.9306))
      {
        fParticleGun -> SetParticleEnergy(255.012*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
    }

  }


  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);        // Michael's has this line. Idk why.

  //----- Setting isotropic particle momentum direction
  G4double alphaMin = 0*deg;    // alpha is apex angle
  G4double alphaMax = 180*deg;  // 180* ensures cone -> full sphere
  G4double cosAlphaMin = cos(alphaMin);
  G4double cosAlphaMax = cos(alphaMax);
  G4double phiMin = 0*deg;      // phi in 0 to 2pi
  G4double phiMax = 360.*deg;   // Presumably the rotation angle. 2pi makes it a cone.

  G4double cosAlpha = cosAlphaMin-G4UniformRand()*(cosAlphaMin-cosAlphaMax);
  G4double sinAlpha = sqrt(1. - cosAlpha*cosAlpha);
  G4double phi = phiMin + G4UniformRand()*(phiMax - phiMin);

  G4double ux = sinAlpha*cos(phi);
  G4double uy = sinAlpha*sin(phi);
  G4double uz = cosAlpha;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));


  //----- Setting the particle generation position
  G4double x0 = 0;              // Says it is negligibly thin.
  G4double y0 = 0;              // Brad told me the source radius
  G4double z0 = 0;

  x0 += fPosOffset[0];		// this is (0,0,0) in my sim and Mendenhall's
  y0 += fPosOffset[1];
  z0 += fPosOffset[2];

  if(fSourceRadius != 0)
  {
    DiskRandom(fSourceRadius, x0, y0);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

}

/*  G4ThreeVector newUz;        // fires isotropic cone where cone axis can be rotated.
  G4double theta, phi, apex;

  G4double xCentre = 0*cm;
  G4double yCentre = 0*cm;
  G4double zCentre = -1*cm;

  phi = (180 + atan(yCentre/xCentre)*(180.0/3.1416))*deg; // basic geometry and then converted to degrees.
  theta = acos(zCentre/sqrt(xCentre*xCentre + yCentre*yCentre + zCentre*zCentre))*(180.0/3.1416)*deg;

  apex = 180*deg;

  newUz = G4ThreeVector(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

  G4double cosAlpha = 1. - G4UniformRand()*(1.- std::cos(apex));
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi      = 2*3.1416*G4UniformRand();  //psi uniform in [0, 2*pi]
  G4ThreeVector dir(sinAlpha*std::cos(psi), sinAlpha*std::sin(psi), cosAlpha);
  dir.rotateUz(newUz);

  fParticleGun->SetParticleMomentumDirection(dir);
*/

