#include	 "BetaSpectrum.hh"

#include	 <complex.h>
#include	 <stdio.h>
#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include	 <TRandom3.h>
#include 	 <TTree.h>
#include	 <TMath.h>
using		 namespace std;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);

int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    cout << "Incorrect format. Execute with: \n";
    cout << "(executable) (number of events) (output file name) (0 for East, 1 for West)" << endl;
    return 0;
  }

  Int_t n_events = atoi(argv[1]);	// converts argument to int. C lib command.
  Int_t polFlag = atoi(argv[3]);	// polarization value
  TString outputName(argv[2]);
  cout << "Saving initial events kinematics in file " << outputName << endl;

  Int_t pol = 10000;
  if(polFlag == 0) { pol = -1; }
  else if(polFlag == 1) { pol = +1; }
  else { cout << "polarization flag is incorrect." << endl; }

  TFile fTree(outputName, "RECREATE");
  TTree* Evts = new TTree("Evts", "initial events kinematics");

  Int_t event_id = -1;	// event ID will be incremented
  Int_t event_ptclID = 11;	// PDG flag 11 means electron

  // randomly throw the kinetic energy
  Double_t test_prob, pdf_value, Te_test, /*theta_test,*/ phi_test, cosTheta_test;
  Double_t event_KE = -1;	// keV
  Double_t event_theta = 0;
  Double_t event_phi = 0;
  TRandom3 factor(0);

  Double_t event_pos[3];	// position in m
  Double_t event_dir[3];        // momentum direction, unit vector
  Double_t event_time;		// time in ns or s, it's all zero anyway
  Double_t event_weight;

  Evts -> Branch("num", &event_id, "event_id/I");
  Evts -> Branch("PID", &event_ptclID, "event_ptclID/I");
  Evts -> Branch("KE", &event_KE, "event_KE/D");
  Evts -> Branch("vertex", event_pos, "event_pos[3]/D");
  Evts -> Branch("direction", event_dir, "event_dir[3]/D");
  Evts -> Branch("time", &event_time, "event_time/D");
  Evts -> Branch("weight", &event_weight, "event_weight/D");

  Double_t normalizer = -1;	// find max value of prob distribution, to normalize beta PDF
  for(int i = 0; i < 9999; i++)
  {
//    Double_t value = neutronCorrectedBetaSpectrum((neutronBetaEp*i)/10000)*(1 - 1*A0_PDG*(beta((neutronBetaEp*i)/10000)));
    Double_t value = neutronCorrectedBetaSpectrum((neutronBetaEp*i)/10000)*(1 + correctedAsymmetry((neutronBetaEp*i)/10000, -1));
    if(value > normalizer)
    {
      normalizer = value;
    }
  }

  Double_t maxPDF = 0;

  // fetch number of events equal to max t in for loop
  for(int t = 0; t < n_events; t++)
  {
    test_prob = 1;
    pdf_value = 0;
    // rejection-acceptance sampling using BetaSpectrum.cc neutronCorrectedBetaSpectrum(...), correctedAsymmetry(...)
    while(test_prob > pdf_value)
    {
      cosTheta_test = 2*factor.Rndm() - 1;	// uniformly sample cos(theta) from -1 to 1

//      theta_test = TMath::ACos(cosTheta_test);
      phi_test = 2*M_PI*factor.Rndm();
      Te_test = neutronBetaEp*factor.Rndm();	// sample random kinetic energy

// Xuan's old change trying to add the Asymmetry term by hand before using MPM's
//      pdf_value = (neutronCorrectedBetaSpectrum(Te_test)*(1 + pol*A0_PDG*beta(Te_test)*TMath::Cos(theta_test))) / normalizer;

      // normalizer is calculated to within 1/10000 precision and used to set the max value of the energy distribution.
      pdf_value = (neutronCorrectedBetaSpectrum(Te_test)*(1 + pol*correctedAsymmetry(Te_test, cosTheta_test))) / normalizer;

      test_prob = factor.Rndm();
    }

    if(pdf_value > maxPDF)
    {
      maxPDF = pdf_value;
    }

//    event_theta = theta_test;
    event_theta = TMath::ACos(cosTheta_test);
    event_phi = phi_test;
    // set momentum direction. Sufficient to use spherical angles since it is unit vector
    event_dir[0] = TMath::Sin(event_theta)*TMath::Cos(event_phi);
    event_dir[1] = TMath::Sin(event_theta)*TMath::Sin(event_phi);
    event_dir[2] = TMath::Cos(event_theta);

    // set the values of each event for the TTree storage.
    event_pos[0] = 1;
    event_pos[1] = 1;
    while(event_pos[0]*event_pos[0] + event_pos[1]*event_pos[1] >= 0.0034129)
    {
      event_pos[0] = 0.05842*(2*factor.Rndm()-1);	// randomly choose x, y between -0.058 and 0.058m
      event_pos[1] = 0.05842*(2*factor.Rndm()-1);	// once it passes the radial cut, exit loop
    }
    event_pos[2] = 1.5*(2*factor.Rndm()-1);	// randomly choose z between -1.5 to 1.5m

    event_KE = Te_test;
    event_ptclID = 11;	// stays as electron
    event_id = t;	// one unique event ID # for each ptcl
    event_time = 0;	// ns or s, doesn't matter since it's zero
    event_weight = 1;	// all equally weighted

    Evts -> Fill();
  }

  if(maxPDF > 1)
  {
    cout << "Badness. Max value of pdf is " << maxPDF << endl;
  }


  fTree.Write();
  fTree.Close();


  cout << "-------------- End of Program ---------------" << endl;
  return 0;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Hits");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
  }

  hPlot -> Draw(command);
}

