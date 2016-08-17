#include "ttree.hh"

#define		N_SD			4
#define		NB_MAX_COINCIDENCES	7
#define		NB_INPUT_FILES		100

// note to self: this code is meant to get hacked up in order to correctly get the GEANT4 sim output into a TTree

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void ReadAndPrint(TString inputFileName, TTree* outTree);

int main()
{
  // loop, sum and store the .txt version
  // You only need this for multiple files.
/*  for(int j = 0; j < NB_INPUT_FILES; j++)
  {
    FillEventInfo(Form("/home/xuansun/Documents/g4_data/10mill_BetaSimData/base/UCNASimOutput_external_Base_C-Geom_%i.txt", j));
    cout << "Event info array is size " << InfoArray.size() << endl;
    PrintToTxtFile(Form("CoincidenceSummed_betaBase_C-Geom_%i.txt", j));
    InfoArray.clear();	//  this resets all the entries so we can count again
  }
*/
  // if you have betas, just read in a file (do not use FillEventInfo since you can have memory leaks)
  // Also, if you are using sources, the PrintToTxtFile argument and ReadAndPrint argument need to be the same

  // Store the TTree version
  for(int t = 0; t < NB_INPUT_FILES; t++)
  {
    TFile file(Form("xuan_analyzed_%i.root", t), "RECREATE");
    // note: for arrays you don't need to use & deallocator.
    // Probably cause the array variable is already a pointer.
    // For primitive types (including wrapper classes _t), you need &.
    TTree* anaTree = new TTree("anaTree", "tree for analysis");
//    ReadAndPrint(Form("CoincidenceSummed_betaBase_C-Geom_%i.txt", t), anaTree);
    ReadAndPrint(Form("/home/xuansun/Documents/g4_data/100mill_BetaSimData/raw_base/UCNASimOutput_external_Base_C-Geom_%i.txt", t), anaTree);
    anaTree -> Write();
    anaTree -> Delete();
  }

  cout << "-------------- End of Program ---------------" << endl;

  return 0;
}




void PrintToTxtFile(TString outName)
{
  cout << "Printing summed data to " << outName << endl;

  // print out all the info line-by-line into output file which is .txt
  ofstream outfile;
  outfile.open(outName, ios::app);
  for(int j = 0; j < InfoArray.size(); j++)
  {
    outfile << InfoArray[j].eid << "\t" << InfoArray[j].sid << "\t" << InfoArray[j].ke << "\t"
            << InfoArray[j].x << "\t" << InfoArray[j].y << "\t" << InfoArray[j].z << "\t"
            << InfoArray[j].px << "\t" << InfoArray[j].py << "\t" << InfoArray[j].pz << "\t"
            << InfoArray[j].t << "\t" << InfoArray[j].w << "\t"
            << InfoArray[j].tFlag << "\t" << InfoArray[j].compT << "\t"
            << InfoArray[j].se_name << "\t" << InfoArray[j].se_e << "\t" << InfoArray[j].se_eq << "\t"
            << InfoArray[j].se_hitT << "\t" << InfoArray[j].se_nTracks << "\t"
            << InfoArray[j].se_ex << "\t" << InfoArray[j].se_ey << "\t" << InfoArray[j].se_ez << "\t"
            << InfoArray[j].se_e2x << "\t" << InfoArray[j].se_e2y << "\t" << InfoArray[j].se_e2z << "\t"
            << InfoArray[j].se_thetaIn << "\t" << InfoArray[j].se_thetaOut << "\t"
            << InfoArray[j].se_keIn << "\t" << InfoArray[j].se_keOut << "\t"
            << InfoArray[j].we_name << "\t" << InfoArray[j].we_e << "\t" << InfoArray[j].we_eq << "\t"
            << InfoArray[j].we_hitT << "\t" << InfoArray[j].we_nTracks << "\t"
            << InfoArray[j].we_ex << "\t" << InfoArray[j].we_ey << "\t" << InfoArray[j].we_ez << "\t"
            << InfoArray[j].we_e2x << "\t" << InfoArray[j].we_e2y << "\t" << InfoArray[j].we_e2z << "\t"
            << InfoArray[j].we_thetaIn << "\t" << InfoArray[j].we_thetaOut << "\t"
            << InfoArray[j].we_keIn << "\t" << InfoArray[j].we_keOut << "\t"
            << InfoArray[j].sw_name << "\t" << InfoArray[j].sw_e << "\t" << InfoArray[j].sw_eq << "\t"
            << InfoArray[j].sw_hitT << "\t" << InfoArray[j].sw_nTracks << "\t"
            << InfoArray[j].sw_ex << "\t" << InfoArray[j].sw_ey << "\t" << InfoArray[j].sw_ez << "\t"
            << InfoArray[j].sw_e2x << "\t" << InfoArray[j].sw_e2y << "\t" << InfoArray[j].sw_e2z << "\t"
            << InfoArray[j].sw_thetaIn << "\t" << InfoArray[j].sw_thetaOut << "\t"
            << InfoArray[j].sw_keIn << "\t" << InfoArray[j].sw_keOut << "\t"
            << InfoArray[j].ww_name << "\t" << InfoArray[j].ww_e << "\t" << InfoArray[j].ww_eq << "\t"
            << InfoArray[j].ww_hitT << "\t" << InfoArray[j].ww_nTracks << "\t"
            << InfoArray[j].ww_ex << "\t" << InfoArray[j].ww_ey << "\t" << InfoArray[j].ww_ez << "\t"
            << InfoArray[j].ww_e2x << "\t" << InfoArray[j].ww_e2y << "\t" << InfoArray[j].ww_e2z << "\t"
            << InfoArray[j].ww_thetaIn << "\t" << InfoArray[j].ww_thetaOut << "\t"
            << InfoArray[j].ww_keIn << "\t" << InfoArray[j].ww_keOut << "\n";

    if(j%100000 == 0)
    {
      cout << "Printing out event " << j << endl;
    }
  }
  outfile.close();

}

void FillEventInfo(TString fileName)
{
  EventInfo evt;        // create an EventInfo object to store everything.

  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  int i = 0;
  bool helper = false;
  vector <int> indexArray;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);

    if(!bufstream.eof())
    {
      bufstream >> evt.eid >> evt.sid >> evt.ke >> evt.x >> evt.y >> evt.z >> evt.px >> evt.py >> evt.pz >> evt.t >> evt.w
                >> evt.tFlag >> evt.compT
                >> evt.se_name >> evt.se_e >> evt.se_eq >> evt.se_hitT >> evt.se_nTracks
                >> evt.se_ex >> evt.se_ey >> evt.se_ez >> evt.se_e2x >> evt.se_e2y >> evt.se_e2z
                >> evt.se_thetaIn >> evt.se_thetaOut >> evt.se_keIn >> evt.se_keOut
                >> evt.we_name >> evt.we_e >> evt.we_eq >> evt.we_hitT >> evt.we_nTracks
                >> evt.we_ex >> evt.we_ey >> evt.we_ez >> evt.we_e2x >> evt.we_e2y >> evt.we_e2z
                >> evt.we_thetaIn >> evt.we_thetaOut >> evt.we_keIn >> evt.we_keOut
                >> evt.sw_name >> evt.sw_e >> evt.sw_eq >> evt.sw_hitT >> evt.sw_nTracks
                >> evt.sw_ex >> evt.sw_ey >> evt.sw_ez >> evt.sw_e2x >> evt.sw_e2y >> evt.sw_e2z
                >> evt.sw_thetaIn >> evt.sw_thetaOut >> evt.sw_keIn >> evt.sw_keOut
                >> evt.ww_name >> evt.ww_e >> evt.ww_eq >> evt.ww_hitT >> evt.ww_nTracks
                >> evt.ww_ex >> evt.ww_ey >> evt.ww_ez >> evt.ww_e2x >> evt.ww_e2y >> evt.ww_e2z
                >> evt.ww_thetaIn >> evt.ww_thetaOut >> evt.ww_keIn >> evt.ww_keOut;

      if(evt.eid == 0)
      {
        InfoArray.push_back(evt);
      }
      else
      {
        for(int t = (InfoArray.size() - NB_MAX_COINCIDENCES); t < InfoArray.size(); t++)
        {
          if(evt.eid == InfoArray[t].eid)
          {
            helper = true;      // it means we found a match
            indexArray.push_back(t);            // save where the match was
          }
        }
        if(helper == true)
        {
          if(indexArray.size() > 0)
          {
            int index = indexArray[0];
	    if(evt.sid == 11)	// check if it's an electron
	    {
		InfoArray[index].sid = 11;	// if in our cascade we got an electron
	    }					// record the whole ptcl ID as an electron

            InfoArray[index].ke += evt.ke;
            InfoArray[index].se_e += evt.se_e;
            InfoArray[index].se_eq += evt.se_eq;
            InfoArray[index].se_ex += evt.se_ex;
            InfoArray[index].se_ey += evt.se_ey;
            InfoArray[index].se_ez += evt.se_ez;
            InfoArray[index].se_e2x += evt.se_e2x;
            InfoArray[index].se_e2y += evt.se_e2y;
            InfoArray[index].se_e2z += evt.se_e2z;
            InfoArray[index].we_e += evt.we_e;
            InfoArray[index].we_eq += evt.we_eq;
            InfoArray[index].we_ex += evt.we_ex;
            InfoArray[index].we_ey += evt.we_ey;
            InfoArray[index].we_ez += evt.we_ez;
            InfoArray[index].we_e2x += evt.we_e2x;
            InfoArray[index].we_e2y += evt.we_e2y;
            InfoArray[index].we_e2z += evt.we_e2z;
            InfoArray[index].sw_e += evt.sw_e;
            InfoArray[index].sw_eq += evt.sw_eq;
            InfoArray[index].sw_ex += evt.sw_ex;
            InfoArray[index].sw_ey += evt.sw_ey;
            InfoArray[index].sw_ez += evt.sw_ez;
            InfoArray[index].sw_e2x += evt.sw_e2x;
            InfoArray[index].sw_e2y += evt.sw_e2y;
            InfoArray[index].sw_e2z += evt.sw_e2z;
            InfoArray[index].ww_e += evt.ww_e;
            InfoArray[index].ww_eq += evt.ww_eq;
            InfoArray[index].ww_ex += evt.ww_ex;
            InfoArray[index].ww_ey += evt.ww_ey;
            InfoArray[index].ww_ez += evt.ww_ez;
            InfoArray[index].ww_e2x += evt.ww_e2x;
            InfoArray[index].ww_e2y += evt.ww_e2y;
            InfoArray[index].ww_e2z += evt.ww_e2z;

            indexArray.pop_back();
          }
          if(indexArray.size() == 0)
          {
            helper = false;
          }
        }

        else
        {
          InfoArray.push_back(evt);
        }

      }

    }
    if(i%100000 == 0)
        cout << "Reading in event " << i << " from file " << fileName << endl;

    if(i == 1500000)
	break;

    i++;
  }
}

void ReadAndPrint(TString inputFileName, TTree* outTree)
{
  Int_t eventID;
  Int_t ptclSpecies;

  Double_t primKE;
  Double_t primPos[3];
  Double_t primMo[3];
  Double_t primTime;
  Double_t primWeight;

  Int_t trapped;
  Double_t compTime;

  // variables related to scintillator crystal and mwpc active region
  Double_t edepScint[2];
  Double_t edepQScint[2];
  Double_t mwpcPos[2][3];       // 2 is for EAST/WEST, 3 is for 3 entries of G4ThreeVector
  Double_t mwpcPos2[2][3];
  Double_t scintPos[2][3];
  Double_t scintPos2[2][3];
  Double_t edepMWPC[2];
  Double_t scintHitTime[2];

  // variables related to all SD's. ATM it's just the above 4. But eventually 24.
  Double_t edepSD[N_SD];
  Double_t thetaInSD[N_SD];
  Double_t thetaOutSD[N_SD];
  Double_t keInSD[N_SD];
  Double_t keOutSD[N_SD];
//  Double_t edepAll;   // can't currently read all on event-by-event basis.
//  Int_t hitCountSD[N_SD];     // counts number of hits per SD. No ability to do this yet.
//  Double_t hitTimeSD[N_SD];   // earliest hit time in each SD

  // initial variables read out by PrimaryGeneratorAction
  outTree -> Branch("primaryEventID", &eventID, "eventID/I");
  outTree -> Branch("primaryParticleSpecies", &ptclSpecies, "ptclSpecies/I");
  outTree -> Branch("primaryKE", &primKE, "primKE/D");
  outTree -> Branch("primaryPosition", primPos, "primPos[3]/D");
  outTree -> Branch("primaryMomentum", primMo, "primMo[3]/D");
  outTree -> Branch("primaryTime", &primTime, "primTime/D");
  outTree -> Branch("primaryWeight", &primWeight, "primWeight/D");

  // global event print out from EventAction
  outTree -> Branch("trapped", &trapped, "trapped/I");
  outTree -> Branch("myCompTime", &compTime, "compTime/D");

  // variables stored by TrackerHit that are printed out by EventAction
  outTree -> Branch("scintillatorEdep", edepScint, "EdepE/D:EdepW/D");
  outTree -> Branch("scintillatorEdepQuenched", edepQScint, "EdepQE/D:EdepQW/D");
//  anaTree -> Branch("all SD edep", &edepAll, "edepAll/D");    // can't currently read
  outTree -> Branch("mwpcEnergy", edepMWPC, "EdepMWPC_E/D:EdepMWPC_W/D");
  outTree -> Branch("scintTimeToHit", scintHitTime, "scintTimeE/D:scintTimeW/D");

  char tmp[1024];
  sprintf(tmp, "EdepSD[%i]/D", N_SD);
  outTree -> Branch("EdepSD", edepSD, tmp);
  sprintf(tmp, "thetaInSD[%i]/D", N_SD);
  outTree -> Branch("thetaInSD", thetaInSD, tmp);
  sprintf(tmp, "thetaOutSD[%i]/D", N_SD);
  outTree -> Branch("thetaOutSD", thetaOutSD, tmp);
  sprintf(tmp, "keInSD[%i]/D", N_SD);
  outTree -> Branch("keInSD", keInSD, tmp);
  sprintf(tmp, "keOutSD[%i]/D", N_SD);
  outTree -> Branch("keOutSD", keOutSD, tmp);

  outTree -> Branch("MWPCPos", mwpcPos, "mwpcPosE[3]/D:mwpcPosW[3]/D");
  outTree -> Branch("ScintPos", scintPos, "scintPosE[3]/D:scintPosW[3]/D");
  outTree -> Branch("MWPCPosSigma", mwpcPos2, "mwpcPos2E[3]/D:mwpcPos2W[3]/D");
  outTree -> Branch("ScintPosSigma", scintPos2, "scintPos2E[3]/D:scintPos2W[3]/D");



  Int_t eid, sid;
  Double_t ke, x, y, z, px, py, pz, t, w;	// units keV, m, m, m, none, none, none, ns, none

  Int_t tFlag;
  Double_t compT;	// units s

  string se_name;
  Double_t se_e, se_eq, se_hitT;	// units keV, keV, s
  Int_t se_nTracks;
  Double_t se_ex, se_ey, se_ez, se_e2x, se_e2y, se_e2z, se_thetaIn, se_thetaOut, se_keIn, se_keOut;
        // units m*keV, m*keV, m*keV, m*m*keV... deg, deg, keV, keV

  string we_name;
  Double_t we_e, we_eq, we_hitT;
  Int_t we_nTracks;
  Double_t we_ex, we_ey, we_ez, we_e2x, we_e2y, we_e2z, we_thetaIn, we_thetaOut, we_keIn, we_keOut;

  string sw_name;
  Double_t sw_e, sw_eq, sw_hitT;        // units keV, keV, s
  Int_t sw_nTracks;
  Double_t sw_ex, sw_ey, sw_ez, sw_e2x, sw_e2y, sw_e2z, sw_thetaIn, sw_thetaOut, sw_keIn, sw_keOut;
        // units m*keV, m*keV, m*keV, m*m*keV... deg, deg, keV, keV

  string ww_name;
  Double_t ww_e, ww_eq, ww_hitT;
  Int_t ww_nTracks;
  Double_t ww_ex, ww_ey, ww_ez, ww_e2x, ww_e2y, ww_e2z, ww_thetaIn, ww_thetaOut, ww_keIn, ww_keOut;

  TString fileName = inputFileName;
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  int counter = 0;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);

    if(!bufstream.eof())
    {

      bufstream >> eid >> sid >> ke >> x >> y >> z >> px >> py >> pz >> t >> w
                >> tFlag >> compT
                >> se_name >> se_e >> se_eq >> se_hitT >> se_nTracks
                >> se_ex >> se_ey >> se_ez >> se_e2x >> se_e2y >> se_e2z >> se_thetaIn >> se_thetaOut >> se_keIn >> se_keOut
                >> we_name >> we_e >> we_eq >> we_hitT >> we_nTracks
                >> we_ex >> we_ey >> we_ez >> we_e2x >> we_e2y >> we_e2z >> we_thetaIn >> we_thetaOut >> we_keIn >> we_keOut
                >> sw_name >> sw_e >> sw_eq >> sw_hitT >> sw_nTracks
                >> sw_ex >> sw_ey >> sw_ez >> sw_e2x >> sw_e2y >> sw_e2z >> sw_thetaIn >> sw_thetaOut >> sw_keIn >> sw_keOut
                >> ww_name >> ww_e >> ww_eq >> ww_hitT >> ww_nTracks
                >> ww_ex >> ww_ey >> ww_ez >> ww_e2x >> ww_e2y >> ww_e2z >> ww_thetaIn >> ww_thetaOut >> ww_keIn >> ww_keOut;


    }

    // related to primaries and event info
    eventID = eid;
    ptclSpecies = sid;
    primKE = ke;
    primPos[0] = x;
    primPos[1] = y;
    primPos[2] = z;
    primMo[0] = px;
    primMo[1] = py;
    primMo[2] = pz;
    primTime = t;
    primWeight = w;
    trapped = tFlag;
    compTime = compT;

    // now we begin recording SD info
    edepScint[0] = se_e;
    edepScint[1] = sw_e;
    edepQScint[0] = se_eq;
    edepQScint[1] = sw_eq;
    scintHitTime[0] = se_hitT;
    scintHitTime[1] = sw_hitT;
    scintPos[0][0] = se_ex; scintPos[0][1] = se_ey; scintPos[0][2] = se_ez;
    scintPos[1][0] = sw_ex; scintPos[1][1] = sw_ey; scintPos[1][2] = sw_ez;
    scintPos2[0][0] = se_e2x; scintPos2[0][1] = se_e2y; scintPos2[0][2] = se_e2z;
    scintPos2[1][0] = sw_e2x; scintPos2[1][1] = sw_e2y; scintPos2[1][2] = sw_e2z;
    edepMWPC[0] = we_e;
    edepMWPC[1] = ww_e;
    mwpcPos[0][0] = we_ex; mwpcPos[0][1] = we_ey; mwpcPos[0][2] = we_ez;
    mwpcPos[1][0] = ww_ex; mwpcPos[1][1] = ww_ey; mwpcPos[1][2] = ww_ez;
    mwpcPos2[0][0] = we_e2x; mwpcPos2[0][1] = we_e2y; mwpcPos2[0][2] = we_e2z;
    mwpcPos2[1][0] = ww_e2x; mwpcPos2[1][1] = ww_e2y; mwpcPos2[1][2] = ww_e2z;


    outTree -> Fill();

    if(counter%100000 == 0)
	cout << "Completed filling TTree with event " << counter << endl;

    counter++;
  }

}
