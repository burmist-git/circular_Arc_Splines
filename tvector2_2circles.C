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

Int_t tvector2_2circles(Double_t xA = -1.1, TString outGifname = "tvector2_2circles_fit_-1.4.gif"){
    
  Double_t xMin = -11.0;
  Double_t xMax =  11.0;
  Double_t yMin = -11.0;
  Double_t yMax =  11.0;

  Double_t r_circle1 = 2.0;
  Double_t r_circle2 = 4.0;
  Double_t x0_circle1 = 0.0;
  Double_t y0_circle1 = 0.0;
  Double_t x0_circle2 = 0.0;
  Double_t y0_circle2 = -2.0;

  TVector2 v_center_circle1(x0_circle1,y0_circle1);
  TVector2 v_center_circle2(x0_circle2,y0_circle2);

  Double_t phi_start_circle1 = 150.0/180.0*TMath::Pi();
  Double_t phi_stop_circle1  =  90.0/180.0*TMath::Pi();
  Double_t phi_start_circle2 =  90.0/180.0*TMath::Pi();
  Double_t phi_stop_circle2  =  30.0/180.0*TMath::Pi();
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
  gr->SetLineWidth(3.0);
  
  ////////////////// Circular Arc Splines /////////////////////////
  //Approximation of Smooth Planar Curves by Circular Arc Splines
  //http://kaj.uniwersytetradom.pl/docs/Biarcs.pdf
  
  Double_t x0; 
  Double_t x2;
  Double_t y0;
  Double_t y2;
  gr->GetPoint(0,x0,y0);
  gr->GetPoint((gr->GetN()-1),x2,y2);
  
  Double_t delta_x_fit = 0.1;
  gr->Fit("pol8","","",x0-delta_x_fit,x2+delta_x_fit);
  TF1 *f_fit = (TF1*)gr->GetListOfFunctions()->FindObject("pol8");
  //Double_t y0=f_fit->Eval(x0);
  //Double_t y2=f_fit->Eval(x2);
  Double_t yDer0=f_fit->Derivative(x0);
  Double_t yDer2=f_fit->Derivative(x2);

  //Double_t xC = (yDer2*x2 - yDer0*x0 - y2 + y0)/(yDer2 - yDer0);
  //Double_t xC = x0 + (yDer2*(x2 - x0) - y2 + y0)/(yDer2 - yDer0);
  Double_t xC = x2 + (yDer0*(x2 - x0) - y2 + y0)/(yDer2 - yDer0);
  Double_t yC = yDer0*(xC - x0) + y0;

  //Double_t xA = 0.4;
  Double_t yA = yDer0*(xA-x0)+y0;
  Double_t D = (xA - x0)*TMath::Sqrt(1.0 + yDer0*yDer0) + (x2 - xA)*TMath::Sqrt(1.0 + yDer2*yDer2);
  Double_t E = y2 - y0 - yDer0*(xA - x0) - yDer2*(x2-xA);
  Double_t F = yDer2*D+TMath::Sqrt(1.0+yDer2*yDer2);
  Double_t m = yDer2 + 2*E*(yDer2*E + TMath::Sqrt(1.0 + yDer2*yDer2)*D)/(D*D - E*E);
  
  Double_t xB = xA + (y2 - yA - yDer2*(x2 - xA))/(m - yDer2);
  Double_t yB = yA + m*(xB - xA);
  //Double_t xB = x2 + (y2 - yA - m*(x2 - xA))/(m - yDer2);
  //Double_t yB = y2 - yDer2*(x2 - xB);

  Double_t x1 = xA + (xA - x0)*TMath::Sqrt((1.0+yDer0*yDer0)/(1.0+m*m));
  Double_t y1 = yA + m*(x1 - xA);
  //
  Double_t q0 = y0 + (x1-x0 + m*(y1-y0))/(m-yDer0);
  Double_t p0 = x0 - yDer0*(q0 - y0);
  Double_t r0 = (q0 - y0)*TMath::Sqrt(1 + yDer0*yDer0);
  //
  Double_t q1 = y2 - (x2-x1 + m*(y2-y1))/(m-yDer2);
  Double_t p1 = x2 - yDer2*(q1 - y2);
  Double_t r1 = (q1 - y2)*TMath::Sqrt(1 + yDer2*yDer2);
  ////////////////// Circular Arc Splines /////////////////////////


  ////////////////// FITED DATA ///////////////////////////////////
  const Int_t nn_fit_circle1 = 10000;
  const Int_t nn_fit_circle2 = 10000;
  TVector2 v_fit;
  TVector2 v_r0_fit(p0,q0);
  TVector2 v_r1_fit(p1,q1);
  Double_t phi_fit;
  Double_t x_fit[nn_fit_circle1+nn_fit_circle2];
  Double_t y_fit[nn_fit_circle1+nn_fit_circle2];
  //
  v_fit.SetMagPhi(TMath::Abs(r0),0.0);
  for(Int_t i = 0;i<nn_fit_circle1;i++){
    phi_fit = (2.0*TMath::Pi() - 0.0)/(nn_fit_circle1 - 1)*i;
    TVector2 v_tmp = (v_fit.Rotate(phi_fit) + v_r0_fit);
    if(v_tmp.X()>=x0 && v_tmp.X()<=x1 && v_tmp.Y()>y2){
      x_fit[i] = v_tmp.X();
      y_fit[i] = v_tmp.Y();
    }
  }
  //
  v_fit.SetMagPhi(TMath::Abs(r1),0.0);
  for(Int_t i = 0;i<nn_fit_circle2;i++){
    phi_fit = (2.0*TMath::Pi() - 0.0)/(nn_fit_circle1 - 1)*i;
    TVector2 v_tmp = (v_fit.Rotate(phi_fit) + v_r1_fit);
    if(v_tmp.X()>x1 && v_tmp.X()<=x2 && v_tmp.Y()>y2){
      x_fit[i+nn_fit_circle1] = v_tmp.X();
      y_fit[i+nn_fit_circle2] = v_tmp.Y();
    }
  }
  //
  TGraph *gr_fit = new TGraph((nn_fit_circle1+nn_fit_circle2),x_fit,y_fit);
  gr_fit->SetLineWidth(3.0);
  gr_fit->SetMarkerStyle(7);
  gr_fit->SetMarkerColor(kBlue+2);

  ////////////////// FITED DATA ///////////////////////////////////  

  ///////////////// support graphs ////////////////////////
  TGraph *gr_support = new TGraph();
  gr_support->SetPoint(0,0,0);
  gr_support->SetPoint(1,x0_circle1,y0_circle1);
  gr_support->SetPoint(2,x0_circle2,y0_circle2);
  gr_support->SetPoint(3,x0,y0);
  gr_support->SetPoint(4,x2,y2);
  gr_support->SetPoint(5,xC,yC);
  gr_support->SetPoint(6,xA,yA);
  gr_support->SetPoint(7,xB,yB);
  gr_support->SetPoint(8,x1,y1);
  //
  gr_support->SetMarkerStyle(20);  
  gr_support->SetMarkerColor(kMagenta);  
  //
  TGraph *gr_pq = new TGraph();  
  gr_pq->SetPoint(0,p0,q0);
  gr_pq->SetPoint(1,p1,q1);
  gr_pq->SetPoint(2,x1,y1);
  //
  gr_pq->SetMarkerStyle(20);  
  gr_pq->SetMarkerColor(kBlue);  
  //
  
  // Tangent lines
  const Int_t nn_tl_P0C = 1000;
  const Int_t nn_tl_P2C = 1000;
  Double_t xMin_tl_P0C = -3.0;
  Double_t xMax_tl_P0C =  2.0; 
  Double_t xMin_tl_P2C =  0.0;
  Double_t xMax_tl_P2C =  5.0; 
  Double_t x_tl_P0C[nn_tl_P0C];
  Double_t y_tl_P0C[nn_tl_P0C]; 
  Double_t x_tl_P2C[nn_tl_P2C];
  Double_t y_tl_P2C[nn_tl_P2C]; 
  //
  for(Int_t i = 0;i<nn_tl_P0C;i++){
    x_tl_P0C[i] = xMin_tl_P0C + (xMax_tl_P0C - xMin_tl_P0C)/(nn_tl_P0C-1)*i;
    y_tl_P0C[i] = yDer0*(x_tl_P0C[i] - x0) + y0; 
  }
  for(Int_t i = 0;i<nn_tl_P2C;i++){
    x_tl_P2C[i] = xMin_tl_P2C + (xMax_tl_P2C - xMin_tl_P2C)/(nn_tl_P2C-1)*i;
    y_tl_P2C[i] = yDer2*(x_tl_P2C[i] - x2) + y2; 
  }
  TGraph *gr_tl_P0C = new TGraph(nn_tl_P0C,x_tl_P0C,y_tl_P0C);
  TGraph *gr_tl_P2C = new TGraph(nn_tl_P2C,x_tl_P2C,y_tl_P2C);

  // The m line
  const Int_t nn_ml = 1000;
  Double_t xMin_ml = -4.0;
  Double_t xMax_ml =  4.0; 
  Double_t x_ml[nn_ml];
  Double_t y_ml[nn_ml]; 
  //
  for(Int_t i = 0;i<nn_tl_P0C;i++){
    x_ml[i] = xMin_ml + (xMax_ml - xMin_ml)/(nn_ml-1)*i;
    y_ml[i] = m*(x_ml[i] - xA) + yA; 
  }
  TGraph *gr_ml = new TGraph(nn_ml,x_ml,y_ml);
  ///////////////// support graphs ////////////////////////

  

  ///////////////// plotting //////////////////////////////

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

  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_support);
  mg->Add(gr_tl_P0C);
  mg->Add(gr_tl_P2C);
  mg->Add(gr_ml);
  mg->Add(gr_pq);
  mg->Add(gr_fit);
  mg->GetXaxis()->SetLimits(xMin,xMax);
  mg->Draw("AP");
  mg->SetMinimum(yMin);
  mg->SetMaximum(yMax);
  
  c1->SaveAs(outGifname.Data());

  ///////////////// plotting //////////////////////////////
    return 0;
}
