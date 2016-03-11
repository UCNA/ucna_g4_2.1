#ifndef TTree_HH
#define TTree_HH

#include         <iostream>
#include         <fstream>
#include         <TGaxis.h>
#include         <sstream>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <stdlib.h>
#include         <TF1.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <math.h>
#include         <TStyle.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <vector>
#include         <string.h>
#include         <fstream>
#include         <TROOT.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <time.h>
#include         <TH2F.h>
#include         <assert.h>
#include         <string>
#include         <TRandom.h>
#include         <TTree.h>
using            namespace std;


void PrintToTxtFile(TString outName);
void FillEventInfo(TString fileName);
void ReadAndPrint(TString inputFileName, TTree* outTree);

struct EventInfo
{
  // set up variables to read in .txt data file
  Int_t eid, sid;
  Double_t ke, x, y, z, px, py, pz, t, w;       // units keV, m, m, m, none, none, none, ns, none

  Int_t tFlag;
  Double_t compT;       // units s

  string se_name;
  Double_t se_e, se_eq, se_hitT;        // units keV, keV, s
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

};


vector <EventInfo> InfoArray;   // Store all the events to track coincidences

#endif
