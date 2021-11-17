//root
#include "TVector2.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TString.h"

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

Int_t generate_2circles_curve(){
    
  Double_t xMin = -7.0;
  Double_t xMax =  7.0;
  Double_t yMin = -7.0;
  Double_t yMax =  7.0;

  Double_t r_circle1 = 2.0;
  Double_t r_circle2 = 4.0;
  Double_t x0_circle1 = 0.0;
  Double_t y0_circle1 = 2.0;
  Double_t x0_circle2 = 0.0;
  Double_t y0_circle2 = 4.0;

  TVector2 v_center_circle1(x0_circle1,y0_circle1);
  TVector2 v_center_circle2(x0_circle2,y0_circle2);

  Double_t phi_start_circle1 = -170.0/180.0*TMath::Pi();
  Double_t phi_stop_circle1  =  -90.0/180.0*TMath::Pi();
  Double_t phi_start_circle2 =  -90.0/180.0*TMath::Pi();
  Double_t phi_stop_circle2  =  -10.0/180.0*TMath::Pi();
  const Int_t nn_circle1 = 1000;
  const Int_t nn_circle2 = 1000;

  ////////////////// INITIAL DATA /////////////////////////
  TVector2 v;
  Double_t phi;
  Double_t x[nn_circle1+nn_circle2];
  Double_t y[nn_circle1+nn_circle2];

  v.SetMagPhi(r_circle1,phi_start_circle1);  
  for(Int_t i = 0;i<nn_circle1;i++){
    phi = (phi_stop_circle1 - phi_start_circle1)/(nn_circle1 - 1)*i;
    TVector2 v_tmp = (v.Rotate(phi) + v_center_circle1);
    x[i] = v_tmp.X();
    y[i] = v_tmp.Y();
  }

  v.SetMagPhi(r_circle2,phi_start_circle2);
  for(Int_t i = 0;i<nn_circle2;i++){
    phi = (phi_stop_circle2 - phi_start_circle2)/(nn_circle2 - 1)*i;
    TVector2 v_tmp = (v.Rotate(phi) + v_center_circle2);
    x[i+nn_circle1] = v_tmp.X();
    y[i+nn_circle1] = v_tmp.Y();
  }

  TGraph *gr = new TGraph((nn_circle1+nn_circle2),x,y);
  gr->SetTitle("");
  gr->SetLineWidth(3.0);  

  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();

  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.05);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.05);
  //
  gr->GetXaxis()->SetLimits(xMin,xMax);
  gr->Draw("AP");
  gr->SetMinimum(yMin);
  gr->SetMaximum(yMax);
  ///////////////// plotting //////////////////////////////

  /////////////////////////////////////////////////////////
  FILE * fp;
  fp = fopen ("curve_inv.dat", "w+");
  for(int i = 0;i<(nn_circle1+nn_circle2);i++)
    fprintf(fp, "%20.10f %20.10f \n",x[i],y[i]);
  fclose(fp);
  /////////////////////////////////////////////////////////
  
  return 0;
}
