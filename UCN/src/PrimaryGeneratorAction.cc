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
  iNbCoincidence(0), iMaxNbCoin(0),
  bCoincidenceWasFired(false)
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
    // Call to method to save primary particle initial info.
    // Need to be super careful here. A Priori, there's no reason that the momentum vector
    // is normalized when I'm reading it from final.
    // And momentum direction is a normalized vector once GEANT4 gets a hold of it.
    // So when we print out what was read in vs. what we get from particle gun may not be the same.
    // Will need to check this. For now it is the same within float to double rounding.
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

//    Set_113SnSource();	// processes all decays of 113Sn that we are interested in
//    Set_139CeSource();
    Set_207BiSource();

    if(bCoincidenceWasFired == false)
    {
      SavePrimPtclInfo2(evtID);
    }
/*    // 113Sn source print-out
    else if(bCoincidenceWasFired == true)
    {
      // save primary ptcl info as same event number so we can add up later
      SavePrimPtclInfo2(evtID - 1);
      // since we've now accounted for one coincidence, lower the number of coincidences in the counter
      iNbCoincidence = iNbCoincidence - 1;
      if(iNbCoincidence < 0)
	G4cout << "RUN-TIME ERROR IN CODE. iNbCoincidence went below 0. Makes no sense." << G4endl;

      if(iNbCoincidence == 0)
      {
        iCoincidencePtcl = -1;		// if we have no more coincidences, reset the ptcl flag to -1
      }
    } */

    // 139Ce print out
    // This is fine. At each event, iNbCoincidence and iMaxNbCoin increment to same number.
    // Hence it can't overlap with iNbCoincidence == 0 unless you have already reset iMaxNbCoin to 0
    // which means as you go through again, it'll switch to to generating a new ptcl, not firing another coincidence
/*    else if(bCoincidenceWasFired == true)
    {
      // check if we've fired any coincidences (and if we need to)
      if(iNbCoincidence == iMaxNbCoin)
      {
	SavePrimPtclInfo2(evtID - 1);	// if we need to, and we haven't done so, do it and lower event # by 1
      }
      else if(iNbCoincidence == iMaxNbCoin - 1)
      {
	SavePrimPtclInfo2(evtID - 2);	// if we have fired 1 already, lower the event number by 2
      }

	// THIS IS GENERALIZABLE. IF THERE ARE MORE THAN 2 COINCIDENCES, KEEP ADDING IF-ELSE STATEMENTS.
	// But do not add more than number of coincidences since iNbCoincidence == 0 is a check later

      // since we've now accounted for one coincidence, lower the number of coincidences in the counter
      iNbCoincidence = iNbCoincidence - 1;
      if(iNbCoincidence < 0)
        G4cout << "RUN-TIME ERROR IN CODE. iNbCoincidence went below 0. Makes no sense." << G4endl;

      if(iNbCoincidence == 0)
      {
        iCoincidencePtcl.clear();
        iMaxNbCoin = 0;                 // this is reset but only used for incidences where more than...
      }                                 // ...1 coincidence ptcl i.e. 139Ce
    }
*/
    // 207Bi
    else if(bCoincidenceWasFired == true)
    {
      // check if we've fired any coincidences (and if we need to)
      if(iNbCoincidence == iMaxNbCoin)
      {
        SavePrimPtclInfo2(evtID - 1);   // if we need to, and we haven't done so, do it and lower event # by 1
      }
      else if(iNbCoincidence == iMaxNbCoin - 1)
      {
        SavePrimPtclInfo2(evtID - 2);   // if we have fired 1 already, lower the event number by 2
      }
      else if(iNbCoincidence == iMaxNbCoin - 2)
      {
        SavePrimPtclInfo2(evtID - 3);
      }
      else if(iNbCoincidence == iMaxNbCoin - 3)
      {
        SavePrimPtclInfo2(evtID - 4);
      }
      else if(iNbCoincidence == iMaxNbCoin - 4)
      {
	SavePrimPtclInfo2(evtID - 5);
      }
      // since we've now accounted for one coincidence, lower the number of coincidences in the counter
      iNbCoincidence = iNbCoincidence - 1;
      if(iNbCoincidence < 0)
        G4cout << "RUN-TIME ERROR IN CODE. iNbCoincidence went below 0. Makes no sense." << G4endl;

      if(iNbCoincidence == 0)
      {
        iCoincidencePtcl.clear();
        iMaxNbCoin = 0;                 // this is reset but only used for incidences where more than...
      }                                 // ...1 coincidence ptcl i.e. 139Ce
    }



  }
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

void PrimaryGeneratorAction::SavePrimPtclInfo2(int eventID)
{
  ofstream outfile;
  outfile.open(sOutputFileName, ios::app);
  outfile << eventID << "\t"
	  << fParticleGun->GetParticleDefinition()->GetPDGEncoding() << "\t"
	  << fParticleGun->GetParticleEnergy()/keV << "\t"
	  << fParticleGun->GetParticlePosition()[0]/m << "\t"
	  << fParticleGun->GetParticlePosition()[1]/m << "\t"
	  << fParticleGun->GetParticlePosition()[2]/m << "\t"
	  << fParticleGun->GetParticleMomentumDirection()[0] << "\t"
	  << fParticleGun->GetParticleMomentumDirection()[1] << "\t"
	  << fParticleGun->GetParticleMomentumDirection()[2] << "\t"
	  << fParticleGun->GetParticleTime()/ns << "\t"
	  << "1 \t";	// this is done since all weights are 1
  outfile.close();

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


void PrimaryGeneratorAction::Set_207BiSource()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  if(iNbCoincidence > 0)
  {
    if(iCoincidencePtcl[iCoincidencePtcl.size() - 1] == 1)           // create a K shell Auger electron for 139Ce
    {
      fParticleGun -> SetParticleEnergy(dCoincidenceEnergies[dCoincidenceEnergies.size() - 1]);
      particle = particleTable -> FindParticle(particleName="e-");
      dCoincidenceEnergies.pop_back();
      iCoincidencePtcl.pop_back();
    }
    else if(iCoincidencePtcl[iCoincidencePtcl.size() - 1] == 2)
    {
      fParticleGun -> SetParticleEnergy(dCoincidenceEnergies[dCoincidenceEnergies.size() - 1]);
      particle = particleTable -> FindParticle(particleName="gamma");
      dCoincidenceEnergies.pop_back();
      iCoincidencePtcl.pop_back();
    }
    bCoincidenceWasFired = true;        // sets the flag stating that we created a ptcl
                                        // that is intended to be recorded in coincidence
  }
  else if(iNbCoincidence == 0)
  {
    bCoincidenceWasFired = false;       // sets the flag stating we are using a coincidence for this event to false
                                        // meaning this is a ptcl generated that is not supposed to be a coincidence.

    G4double dr1 = G4UniformRand()*99.7886307;	// gives us a range from (0, 99.7886307)
						// total prob of headed to 3 branches for 207Bi
						// The upper 2 are calculated from nndc data, the 1st excited
						// is taken from the less precise decay scheme on nndc
    // headed to the stable branch (this entire thing is 1 cascade of coincidences)
    if((dr1 >= 0) && (dr1 <= 83.86))
    {
      G4double dr5 = G4UniformRand()*100;
      if((dr5 >= 0) && (dr5 <= 2.50733))	// fire an Auger without any other coincidences (since stable level)
      {
        fParticleGun -> SetParticleEnergy(56.7*keV);
        particle = particleTable->FindParticle(particleName="e-");
      }
      else	// singular Auger was not fired
      {	// so now we're in the 3rd excited level which is stable

	// if gamma, then the remainder of the code is the prob. split of decay from 1st -> GS
        if((dr1 >= 0) && (dr1 <= 74.5))	// if gamma, then we do the entire decay from 1st excited
	{
	  fParticleGun -> SetParticleEnergy(1063.656*keV);
	  particle = particleTable->FindParticle(particleName="gamma");

	  G4double dr6 = G4UniformRand()*99.84;	// percentage of decays from 1st excited to GS
	  if((dr6 >= 0) && (dr6 <= 97.75))
	  {
	    iCoincidencePtcl.push_back(2);
	    dCoincidenceEnergies.push_back(569.698*keV);
     	    iNbCoincidence++;
     	    iMaxNbCoin++;
	  }
	  else if((dr6 > 97.75) && (dr6 <= 99.84))
	  {
	    // i.e. if we are firing a CE, check if we get an Auger
	    G4double dr7 = G4UniformRand()*100;
	    if((dr7 >= 0) && (dr7 <= 2.50733)) // sets an Auger in coincidence with init gamma ray + second CE
	    {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(56.7*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
	    }

	    // since we are firing a CE, check the spread on energies
	    if((dr6 > 97.75) && (dr6 <= 99.287))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(481.6935*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
            else if((dr6 > 99.287) && (dr6 <= 99.729))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(553.8372*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
            else if((dr6 > 99.729) && (dr6 <= 99.84))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(565.8473*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
	  }
	}
	// if CE, we do the exact same thing. Except we also need to check for an Auger coincidence beforehand.
	else if((dr1 > 74.5) && (dr1 <= 83.86))
	{
	  // CE is taking us from 3rd excited to 1st excited. Check for Auger
          G4double dr17 = G4UniformRand()*100;
          if((dr17 >= 0) && (dr17 <= 2.50733))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(56.7*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }

	  // regardless of whether an Auger was fired or not, now sample the CE that do 3->1
	  if((dr1 > 74.5) && (dr1 <= 74.94))
	  {	// this is the first particle that decays from this branch.
		// So it needs to be set explicitly and not as a coincidence ptcl
	    fParticleGun -> SetParticleEnergy(1059.805*keV);
	    particle = particleTable -> FindParticle(particleName="e-");
	  }
	  else if((dr1 > 74.94) && (dr1 <= 76.78))
	  {
	    fParticleGun -> SetParticleEnergy(1047.795*keV);
	    particle = particleTable -> FindParticle(particleName="e-");
	  }
	  else if((dr1 > 76.78) && (dr1 <= 83.86))
	  {
	    fParticleGun -> SetParticleEnergy(975.651*keV);
	    particle = particleTable -> FindParticle(particleName="e-");
	  }
	  // Now we are in 1st excited via CE decay from 3rd excited, with possibly an additional Auger released.
	  // What remains is to put the entire 1st -> GS as a series of coincidence ptcls i.e. a cascade
	  G4double dr18 = G4UniformRand()*99.84;	// percentage of decays from 1st excited to GS
	  if((dr18 >= 0) && (dr18 <= 97.75))
	  {
	    iCoincidencePtcl.push_back(2);
	    dCoincidenceEnergies.push_back(569.698*keV);
     	    iNbCoincidence++;
     	    iMaxNbCoin++;
	  }
	  else if((dr18 > 97.75) && (dr18 <= 99.84))
	  {
	    // i.e. if we are firing a CE, from 1st to GS, check if we get an Auger
	    G4double dr19 = G4UniformRand()*100;
	    if((dr19 >= 0) && (dr19 <= 2.50733)) // sets an Auger in coincidence with init gamma ray + second CE
	    {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(56.7*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
	    }

	    // since we are firing a CE, from 1st to GS, check the spread on energies
	    if((dr18 > 97.75) && (dr18 <= 99.287))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(481.6935*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
            else if((dr18 > 99.287) && (dr18 <= 99.729))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(553.8372*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
            else if((dr18 > 99.729) && (dr18 <= 99.84))
            {
              iCoincidencePtcl.push_back(1);
              dCoincidenceEnergies.push_back(565.8473*keV);
              iNbCoincidence++;
              iMaxNbCoin++;
            }
	  }
	}
      }
    }
    // headed to the upper branch
    else if((dr1 > 83.86) && (dr1 <= 90.8886307))
    {
      // on the way to upper branch, check if we get an Auger electron
      G4double dr8 = G4UniformRand()*100;
      if((dr8 >= 0) && (dr8 <= 2.50733)) // sets an Auger in coincidence with neutron capture decay to 4th excited state
      {	// note this is added as a coincidence since we are guaranteed a ptcl in the decay down from 4th excited
        iCoincidencePtcl.push_back(1);
        dCoincidenceEnergies.push_back(56.7*keV);
        iNbCoincidence++;
        iMaxNbCoin++;
      }

      if((dr1 > 83.86) && (dr1 <= 90.73))
      {	// gamma takes us from 4th excited to 1st excited
        fParticleGun -> SetParticleEnergy(1770.228*keV);
        particle = particleTable->FindParticle(particleName="gamma");

	// so now we're in 1st excited, need to set a cascade of coincidences from 1st excited down to GS
	G4double dr13 = G4UniformRand()*99.84;
	if((dr13 >= 0) && (dr13 <= 97.75))
	{	// if true, a gamma takes us to GS. So no Auger coincidence.
	  iCoincidencePtcl.push_back(2);
	  dCoincidenceEnergies.push_back(569.698*keV);
	  iNbCoincidence++;
	  iMaxNbCoin++;
	}
	else if((dr13 > 97.75) && (dr13 <= 99.84))
	{
	  G4double dr12 = G4UniformRand()*100;
	  // since a CE drops us from 1st excited to GS, check for an Auger in coincidence
	  if((dr12 >= 0) && (dr12 <= 2.50733))
	  {
	    iCoincidencePtcl.push_back(1);
	    dCoincidenceEnergies.push_back(56.7*keV);
	    iNbCoincidence++;
	    iMaxNbCoin++;
	  }

	  // now, evaluate whether that CE was M, L or K shell (in that order)
	  if((dr13 > 97.75) && (dr13 <= 97.861))
	  {
	    iCoincidencePtcl.push_back(1);
	    dCoincidenceEnergies.push_back(565.8473*keV);
	    iNbCoincidence++;
	    iMaxNbCoin++;
	  }
	  else if((dr13 > 97.861) && (dr13 <= 98.303))
	  {
	    iCoincidencePtcl.push_back(1);
	    dCoincidenceEnergies.push_back(553.8372*keV);
	    iNbCoincidence++;
	    iMaxNbCoin++;
	  }
	  else if((dr13 > 98.303) && (dr13 <= 99.84))
	  {
	    iCoincidencePtcl.push_back(1);
	    dCoincidenceEnergies.push_back(481.6935);
	    iNbCoincidence++;
	    iMaxNbCoin++;
	  }
	}

      }
      else if((dr1 > 90.73) && (dr1 <= 90.7572))
      {	// if we fire a CE in the decay from 4th excited to 1st excited, check if we get an Auger in coincidence
        G4double dr9 = G4UniformRand()*100;
        if((dr9 >= 0) && (dr9 <= 2.50733)) // sets an Auger in coincidence with neutron capture decay to 4th excited state
        {
          iCoincidencePtcl.push_back(1);
          dCoincidenceEnergies.push_back(56.7*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
        }
	// these are the CE we could fire to drop from 4th to 1st. L and K shells
        if((dr1 > 90.73) && (dr1 <= 90.7334))
        {
          fParticleGun -> SetParticleEnergy(1754.367*keV);
	  particle = particleTable -> FindParticle(particleName="e-");
        }
	else if((dr1 > 90.7334) && (dr1 <= 90.7572))
	{
	  fParticleGun -> SetParticleEnergy(1682.224*keV);
	  particle = particleTable -> FindParticle(particleName="e-");
	}

        // so now we're in 1st excited, need to set a cascade of coincidences from 1st excited down to GS
        G4double dr11 = G4UniformRand()*99.84;
        if((dr11 >= 0) && (dr11 <= 97.75))
        {       // if true, a gamma takes us to GS. So no Auger coincidence.
          iCoincidencePtcl.push_back(2);
          dCoincidenceEnergies.push_back(569.698*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
        }
        else if((dr11 > 97.75) && (dr11 <= 99.84))
        {
          G4double dr12 = G4UniformRand()*100;
          // since a CE drops us from 1st excited to GS, check for an Auger in coincidence
          if((dr12 >= 0) && (dr12 <= 2.50733))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(56.7*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }

          // now, evaluate whether that CE was M, L or K shell (in that order)
          if((dr11 > 97.75) && (dr11 <= 97.861))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(565.8473*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          else if((dr11 > 97.861) && (dr11 <= 98.303))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(553.8372*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          else if((dr11 > 98.303) && (dr11 <= 99.84))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(481.6935);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
        }
      }
      // this is 4th excited to 2nd excited. We'll ignore the 2->1 excited decay since the probability is so small.
      // So we'll only code the 2nd to GS decay (about 0.13% of decays).
      else if((dr1 > 90.7572) && (dr1 <= 90.8882))
      {	// gamma is taking us from 4th excited to 2nd excited state
	fParticleGun -> SetParticleEnergy(1442.2*keV);
	particle = particleTable -> FindParticle(particleName="e-");

	// now we're in 2nd excited. A gamma or a CE drops us to ground state (GS)
	G4double dr14 = G4UniformRand()*0.130962;	// total intensities sum for 2nd->GS decays
	if((dr14 >= 0) && (dr14 <= 0.128))
	{	// fire a gamma in coincidence with w/e ptcl got us to 2nd excited state
          iCoincidencePtcl.push_back(2);
          dCoincidenceEnergies.push_back(897.77*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
	}
	else if((dr14 > 0.128) && (dr14 <= 0.130962))
	{	// this means a CE is taking us to GS. So check for Augers
	  G4double dr15 = G4UniformRand()*100;
	  if((dr15 >= 0) && (dr15 <= 2.50733))
	  {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(56.7*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
	  }
	  // regardless of whether an Auger was produced or not (in coincidence), check which CE
	  // will drop us to GS, which needs to be fired in coincidence as well
	  if((dr14 > 0.128) && (dr14 <= 0.128095))
	  {
	    iCoincidencePtcl.push_back(1);
	    dCoincidenceEnergies.push_back(893.92*keV);
	    iNbCoincidence++;
	    iMaxNbCoin++;
	  }
	  else if((dr14 > 0.128095) && (dr14 <= 0.128502))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(881.91*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          else if((dr14 > 0.128502) && (dr14 <= 0.130962))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(809.77*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
	}
      }
      else if((dr1 > 90.8882) && (dr1 <= 90.8886307))
      {	// if we fire a CE in the decay from 4th excited to 2nd excited, check if we get an Auger in coincidence
        G4double dr10 = G4UniformRand()*100;
        if((dr10 >= 0) && (dr10 <= 2.50733)) // sets an Auger in coincidence with neutron capture decay to 4th excited state
        {
          iCoincidencePtcl.push_back(1);
          dCoincidenceEnergies.push_back(56.7*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
        }
	// these are the CE we could fire to drop from 4th to 2nd. M, L, K shells
	if((dr1 > 90.8882) && (dr1 <= 90.8882144))
	{
	  fParticleGun -> SetParticleEnergy(1438.35*keV);
	  particle = particleTable -> FindParticle(particleName="e-");
	}
	else if((dr1 > 90.8882144) && (dr1 <= 90.8882757))
	{
	  fParticleGun -> SetParticleEnergy(1426.34*keV);
	  particle = particleTable -> FindParticle(particleName="e-");
	}
	else if((dr1 > 90.8882757) && (dr1 <= 90.8886307))
	{
	  fParticleGun -> SetParticleEnergy(1354.20*keV);
	  particle = particleTable -> FindParticle(particleName="e-");
	}
	// At this point, in this branch, a CE has taken us from 4th to 2nd and (possibly) and Auger fired in coincidence.
        // now we're in 2nd excited. A gamma or a CE drops us to ground state (GS)
        G4double dr14 = G4UniformRand()*0.130962;       // total intensities sum for 2nd->GS decays
        if((dr14 >= 0) && (dr14 <= 0.128))
        {       // fire a gamma in coincidence with w/e ptcl got us to 2nd excited state
          iCoincidencePtcl.push_back(2);
          dCoincidenceEnergies.push_back(897.77*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
        }
        else if((dr14 > 0.128) && (dr14 <= 0.130962))
        {       // this means a CE is taking us to GS. So check for Augers
          G4double dr15 = G4UniformRand()*100;
          if((dr15 >= 0) && (dr15 <= 2.50733))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(56.7*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          // regardless of whether an Auger was produced or not (in coincidence), check which CE
          // will drop us to GS, which needs to be fired in coincidence as well
          if((dr14 > 0.128) && (dr14 <= 0.128095))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(893.92*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          else if((dr14 > 0.128095) && (dr14 <= 0.128502))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(881.91*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
          else if((dr14 > 0.128502) && (dr14 <= 0.130962))
          {
            iCoincidencePtcl.push_back(1);
            dCoincidenceEnergies.push_back(809.77*keV);
            iNbCoincidence++;
            iMaxNbCoin++;
          }
        }
      }
    }
    // going to lowest branch (1st excited state).
    // this branching ratio is taken from the less precise decay scheme since I didn't want to calculate it.
    // The percentages will be the fraction of decays given that you are in 1st excited (since all decays end up there)
    else if((dr1 > 90.8886307) && (dr1 <= 99.7886307))
    {
      // since this neutron capture is decaying to 1st excited, check for an Auger
      G4double dr20 = G4UniformRand()*100;
      if((dr20 >= 0) && (dr20 <= 2.50733))
      {
        iCoincidencePtcl.push_back(1);
        dCoincidenceEnergies.push_back(56.7*keV);
        iNbCoincidence++;
        iMaxNbCoin++;
      }

      // if we get an Auger, add it in coincidence to the coincidence stack
      // but regardless, then process the decay down to GS
      G4double percent1stToGS = G4UniformRand()*99.84;	// percentages taken from nndc for 1st -> GS
      if((percent1stToGS >= 0) && (percent1stToGS <= 97.75))
      {
	fParticleGun -> SetParticleEnergy(569.698*keV);
	particle = particleTable->FindParticle(particleName = "gamma");
      }
      else if((percent1stToGS > 97.75) && (percent1stToGS < 99.84))
      {	// decay 1st to GS via CE. Check for an Auger.
        G4double dr21 = G4UniformRand()*100;
        if((dr21 >= 0) && (dr21 <= 2.50733))
        {
          iCoincidencePtcl.push_back(1);
          dCoincidenceEnergies.push_back(56.7*keV);
          iNbCoincidence++;
          iMaxNbCoin++;
        }


	// then, process the CE from 1st to GS
        if((percent1stToGS > 97.75) && (percent1stToGS <= 97.861))
        {
	  fParticleGun -> SetParticleEnergy(565.8473*keV);
	  particle = particleTable->FindParticle(particleName = "e-");
        }
        else if((percent1stToGS > 97.861) && (percent1stToGS <= 98.303))
        {
	  fParticleGun -> SetParticleEnergy(553.8372*keV);
	  particle = particleTable->FindParticle(particleName = "e-");
        }
        else if((percent1stToGS > 98.303) && (percent1stToGS <= 99.84))
        {
	  fParticleGun -> SetParticleEnergy(481.6935*keV);
	  particle = particleTable->FindParticle(particleName = "e-");
        }
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

  x0 += fPosOffset[0];          // this is (0,0,0) in my sim and Mendenhall's
  y0 += fPosOffset[1];
  z0 += fPosOffset[2];

  if(fSourceRadius != 0)
  {
    DiskRandom(fSourceRadius, x0, y0);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
}


void PrimaryGeneratorAction::Set_139CeSource()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  G4double percentage = G4UniformRand()*99.9231;	// total percentage of stuff we're interested in for 139Ce

  if(iNbCoincidence > 0)
  {
    if(iCoincidencePtcl[iCoincidencePtcl.size() - 1] == 1)           // create a K shell Auger electron for 139Ce
    {
      fParticleGun -> SetParticleEnergy(dCoincidenceEnergies[dCoincidenceEnergies.size() - 1]);
      particle = particleTable -> FindParticle(particleName="e-");
      dCoincidenceEnergies.pop_back();
      iCoincidencePtcl.pop_back();
    }
    // technically, here you would put other options for other coincidences.
    // for 139Ce, there is no other coincidence that we care about.

    bCoincidenceWasFired = true;        // sets the flag stating that we created a ptcl
                                        // that is intended to be recorded in coincidence
  }
  else if(iNbCoincidence == 0)
  {
    bCoincidenceWasFired = false;       // sets the flag stating we are using a coincidence for this event to false
                                        // meaning this is a ptcl generated that is not supposed to be a coincidence.
    					// since the 1st excited state isn't very long lived
    G4double AugerPercentFromInitDecay = G4UniformRand()*100;	// need to check if we get an Auger in first decay.
    if((AugerPercentFromInitDecay >= 0) && (AugerPercentFromInitDecay <= 6.92))
    {
      iCoincidencePtcl.push_back(1);		// if we do, then maximum # of ptcl's fired in coincidence is 3
      dCoincidenceEnergies.push_back(27.4*keV);
      iNbCoincidence++;
      iMaxNbCoin++;
    }


    if((percentage >= 0) && (percentage <= 80))
    {   // fires a gamma in the decay from 1st excited to ground state
      fParticleGun -> SetParticleEnergy(391.698*keV);
      particle = particleTable->FindParticle(particleName="gamma");
    }
    else if((percentage > 80) && (percentage <= 99.9231))
    {   // produces a Conversion Electron from 1st excited state to GS

      double AugerPercent = G4UniformRand()*100;  // check whether the CE also comes with an Auger
      if((AugerPercent >= 0) && (AugerPercent <= 6.92))
      {
        iCoincidencePtcl.push_back(1);   // sets the flag to produce an Auger
	dCoincidenceEnergies.push_back(27.4*keV);
        iNbCoincidence++;       // increments our coincidence counter so next event, we produce an Auger
	iMaxNbCoin++;
      }

      // whether or not the Auger flag checks out, we produce all the CE down here.
      // These are, in order, the K, L, M shell Conversion Electron energies.
      if((percentage > 80) && (percentage <= 97.15))
      {
        fParticleGun -> SetParticleEnergy(126.9329*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 97.15) && (percentage <= 99.448))
      {
        fParticleGun -> SetParticleEnergy(159.5912*keV);
        particle = particleTable -> FindParticle(particleName="e-");
      }
      else if((percentage > 99.448) && (percentage <= 99.9231))
      {
        fParticleGun -> SetParticleEnergy(164.4962*keV);
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

  x0 += fPosOffset[0];          // this is (0,0,0) in my sim and Mendenhall's
  y0 += fPosOffset[1];
  z0 += fPosOffset[2];

  if(fSourceRadius != 0)
  {
    DiskRandom(fSourceRadius, x0, y0);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
}


void PrimaryGeneratorAction::Set_113SnSource()	// don't need additional arguments since we set the particle gun.
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  double percentage = G4UniformRand()*102.9306;	// total percentage of things we're interested in for 113Sn

  if(iNbCoincidence > 0)
  {
    if(iCoincidencePtcl[iCoincidencePtcl.size() - 1] == 1)		// create a K shell Auger electron for 113Sn
    {
      fParticleGun -> SetParticleEnergy(dCoincidenceEnergies[dCoincidenceEnergies.size() - 1]);
      particle = particleTable -> FindParticle(particleName="e-");
      dCoincidenceEnergies.pop_back();
      iCoincidencePtcl.pop_back();
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

      double AugerPercent = G4UniformRand()*100;	// check whether the CE also comes with an Auger
      if((AugerPercent >= 0) && (AugerPercent <= 11.8056))
      {
        iCoincidencePtcl.push_back(1);	// sets the flag to produce an Auger
        dCoincidenceEnergies.push_back(20.1*keV);
        iNbCoincidence++;	// increments our coincidence counter so next event, we produce an Auger
        iMaxNbCoin++;
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

