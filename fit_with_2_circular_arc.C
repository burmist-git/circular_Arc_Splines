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
#include <iomanip>
#include <vector>

#include <time.h>

using namespace std;


////////////////// Circular Arc Splines //////////////////////////
//Approximation of Smooth Planar Curves by Circular Arc Splines //
//http://kaj.uniwersytetradom.pl/docs/Biarcs.pdf                //
////////////////// Circular Arc Splines //////////////////////////

//Double_t xMin = -11.0;
//Double_t xMax =  11.0;
//Double_t yMin = -11.0;
//Double_t yMax =  11.0;

Double_t xMin = -310.0;
Double_t xMax = 310.0;
Double_t yMin = 0.0;
Double_t yMax = 25.0;

void plot_gr( TGraph *gr, TCanvas *c1);
Double_t getRMSERR(TF1 *f_data,TF1 *f_fit, Double_t xmin, Double_t xmax, Int_t nn);

struct fit_state_str {
  Double_t x0;
  Double_t x2;
  Double_t y0;
  Double_t y2;
  Double_t xC;
  Double_t yC;
  Double_t xA;
  Double_t yA;
  Double_t xB;
  Double_t yB;
  Double_t x1;
  Double_t y1;
  //
  Double_t yDer0;
  Double_t yDer2;
  Double_t m;
  //
  Double_t q0;
  Double_t p0;
  Double_t r0;
  Double_t q1;
  Double_t p1;
  Double_t r1;
  //
  Double_t rmserr;
  //
  std::vector<Double_t> vec_x_fit;
  std::vector<Double_t> vec_y_fit;
  fit_state_str(){
  }
};

void fit(TF1 *f_fit, Double_t x0, Double_t x2, const Int_t nSteps, Int_t nnRMSERR, Double_t *rmserr_arr, Double_t *xA_arr, fit_state_str *fit_state);
void print_str(fit_state_str *fit_state);

Int_t fit_with_2_circular_arc( TString dataFileIn = "curve.dat", TString outGifname = "curve.dat.gif"){

  //dataFileIn = "curve_inv.dat";
  //outGifname = "curve_inv.dat.gif";
  
  //dataFileIn = "lens_profile_n1.4.dat";
  //outGifname = "lens_profile_n1.4.dat.gif";
  dataFileIn = "lens_profile_n1.4_long.dat";
  outGifname = "lens_profile_n1.4_long.dat.gif";
  
  //////////////////////// LOAD DATA ////////////////////////////
  TGraph *gr = new TGraph();
  TGraph *gr_nofit = new TGraph();
  gr->SetTitle("");
  gr->SetLineWidth(3.0);
  gr_nofit->SetTitle("");
  gr_nofit->SetLineWidth(3.0);
  Double_t x,y;
  ifstream fileIn (dataFileIn.Data());
  if (fileIn.is_open()) {
    while ( fileIn>>x>>y ) {
      gr->SetPoint(gr->GetN(),x,y);
      gr_nofit->SetPoint(gr->GetN(),x,y);
    }
    fileIn.close();
  }
  else cout<<"Unable to open file"<<endl; 
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  plot_gr(gr,c1);
  
  //
  Double_t pL_x;
  Double_t pL_y;
  Double_t pR_x;
  Double_t pR_y;
  //
  gr->GetPoint(0,pL_x,pL_y);
  gr->GetPoint((gr->GetN()-1),pR_x,pR_y);
  gr->Fit("pol8","Q","",pL_x,pR_x);
  TF1 *f_fit = (TF1*)gr->GetListOfFunctions()->FindObject("pol8");
  //Double_t x0 = -1.6;
  //Double_t x2 =  3.0;
  Double_t x0 = -295.0;
  Double_t x2 =  -75.0;

  //////////////////////// LOAD DATA ////////////////////////////

  const Int_t nSteps = 400;
  Int_t nnRMSERR = 100;
  //  
  Double_t rmserr_arr[nSteps];
  Double_t xA_arr[nSteps];
  //
  fit_state_str *fit_state = new fit_state_str();
  //
  fit( f_fit, x0, x2, nSteps, nnRMSERR, rmserr_arr, xA_arr, fit_state);
  print_str(fit_state);
 
  ///////////////// plotting //////////////////////////////
  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TGraph *gr_rmse_xA = new TGraph(nSteps,xA_arr,rmserr_arr);
  gr_rmse_xA->SetTitle("");
  gr_rmse_xA->Draw("APL");

  //
  TCanvas *c3 = new TCanvas("c3","c3",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.05);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.05);


  ///////////////// support graphs ////////////////////////
  TGraph *gr_support = new TGraph();
  gr_support->SetPoint(0,0,0);
  //gr_support->SetPoint(1,x0_circle1,y0_circle1);
  //gr_support->SetPoint(2,x0_circle2,y0_circle2);
  gr_support->SetPoint(1,fit_state->x0,fit_state->y0);
  gr_support->SetPoint(2,fit_state->x2,fit_state->y2);
  gr_support->SetPoint(3,fit_state->xC,fit_state->yC);
  gr_support->SetPoint(4,fit_state->xA,fit_state->yA);
  gr_support->SetPoint(5,fit_state->xB,fit_state->yB);
  gr_support->SetPoint(6,fit_state->x1,fit_state->y1);
  //
  gr_support->SetMarkerStyle(20);  
  gr_support->SetMarkerColor(kMagenta);  

  TGraph *gr_pq = new TGraph();  
  gr_pq->SetPoint(0,fit_state->p0,fit_state->q0);
  gr_pq->SetPoint(1,fit_state->p1,fit_state->q1);
  gr_pq->SetPoint(2,fit_state->x1,fit_state->y1);
  //
  gr_pq->SetMarkerStyle(20);  
  gr_pq->SetMarkerColor(kBlue);  
  //
  
  // Tangent lines
  const Int_t nn_tl_P0C = 10000;
  const Int_t nn_tl_P2C = 10000;
  Double_t xMin_tl_P0C = xMin;
  Double_t xMax_tl_P0C = xMax; 
  Double_t xMin_tl_P2C = xMin;
  Double_t xMax_tl_P2C = xMax; 
  Double_t x_tl_P0C[nn_tl_P0C];
  Double_t y_tl_P0C[nn_tl_P0C]; 
  Double_t x_tl_P2C[nn_tl_P2C];
  Double_t y_tl_P2C[nn_tl_P2C]; 
  //
  for(Int_t j = 0;j<nn_tl_P0C;j++){
    x_tl_P0C[j] = xMin_tl_P0C + (xMax_tl_P0C - xMin_tl_P0C)/(nn_tl_P0C-1)*j;
    y_tl_P0C[j] = fit_state->yDer0*(x_tl_P0C[j] - fit_state->x0) + fit_state->y0; 
  }
  for(Int_t j = 0;j<nn_tl_P2C;j++){
    x_tl_P2C[j] = xMin_tl_P2C + (xMax_tl_P2C - xMin_tl_P2C)/(nn_tl_P2C-1)*j;
    y_tl_P2C[j] = fit_state->yDer2*(x_tl_P2C[j] - fit_state->x2) + fit_state->y2; 
  }
  TGraph *gr_tl_P0C = new TGraph(nn_tl_P0C,x_tl_P0C,y_tl_P0C);
  TGraph *gr_tl_P2C = new TGraph(nn_tl_P2C,x_tl_P2C,y_tl_P2C);
  
  // The m line
  const Int_t nn_ml = 10000;
  Double_t xMin_ml = xMin;
  Double_t xMax_ml = xMax; 
  Double_t x_ml[nn_ml];
  Double_t y_ml[nn_ml]; 
  //
  for(Int_t j = 0;j<nn_ml;j++){
    x_ml[j] = xMin_ml + (xMax_ml - xMin_ml)/(nn_ml-1)*j;
    y_ml[j] = fit_state->m*(x_ml[j] - fit_state->xA) + fit_state->yA; 
  }
  TGraph *gr_ml = new TGraph(nn_ml,x_ml,y_ml);
  ///////////////// support graphs ////////////////////////

  ///////////////// the fit ///////////////////////////////
  TGraph *gr_fit = new TGraph();
  gr_fit->SetTitle("");
  gr_fit->SetMarkerColor(kBlue);
  for(unsigned int ii = 0;ii<fit_state->vec_x_fit.size();ii++){
    gr_fit->SetPoint(gr_fit->GetN(),fit_state->vec_x_fit.at(ii),fit_state->vec_y_fit.at(ii));
  }
  ///////////////// plotting //////////////////////////////
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_nofit);
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
  
  //c1->SaveAs(outGifname.Data());

  ///////////////// plotting //////////////////////////////

  return 0;  
  
}

void fit(TF1 *f_fit, Double_t x0, Double_t x2, const Int_t nSteps,Int_t nnRMSERR,Double_t *rmserr_arr, Double_t *xA_arr, fit_state_str *fit_state){
  Double_t y0=f_fit->Eval(x0);
  Double_t y2=f_fit->Eval(x2);
  Double_t yDer0=f_fit->Derivative(x0);
  Double_t yDer2=f_fit->Derivative(x2);
  //
  Double_t xC = x2 + (yDer0*(x2 - x0) - y2 + y0)/(yDer2 - yDer0);
  Double_t yC = yDer0*(xC - x0) + y0;
  //
  Double_t xA;
  Double_t xA_min = x0;
  Double_t xA_max = xC;
  //
  for(Int_t ii = 0;ii<nSteps;ii++){
    //for(Int_t ii = 50;ii<51;ii++){  
    xA = xA_min + (xA_max - xA_min)/(nSteps - 1)*ii;
    //xA = -1.1;
    Double_t yA = yDer0*(xA-x0)+y0;
    Double_t D = (xA - x0)*TMath::Sqrt(1.0 + yDer0*yDer0) + (x2 - xA)*TMath::Sqrt(1.0 + yDer2*yDer2);
    Double_t E = y2 - y0 - yDer0*(xA - x0) - yDer2*(x2-xA);
    Double_t F = yDer2*D+TMath::Sqrt(1.0+yDer2*yDer2);
    Double_t m = yDer2 + 2*E*(yDer2*E + TMath::Sqrt(1.0 + yDer2*yDer2)*D)/(D*D - E*E);
    //
    Double_t xB = xA + (y2 - yA - yDer2*(x2 - xA))/(m - yDer2);
    Double_t yB = yA + m*(xB - xA);
    //    
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

    ////////////////// FITED DATA ///////////////////////////////////
    const Int_t nn_fit_circle1 = 1000;
    const Int_t nn_fit_circle2 = 1000;
    TVector2 v_fit;
    TVector2 v_r0_fit(p0,q0);
    TVector2 v_r1_fit(p1,q1);
    Double_t phi_fit;
    //
    TVector2 v_r0_fit_left;
    TVector2 v_r0_fit_right;
    TVector2 v_r1_fit_left;
    TVector2 v_r1_fit_right;
    v_r0_fit_left.Set(x0,y0);
    v_r0_fit_right.Set(x1,y1);
    v_r1_fit_left.Set(x1,y1);
    v_r1_fit_right.Set(x2,y2);
    TVector2 v_tmp_tmp;
    v_tmp_tmp = v_r0_fit_left - v_r0_fit;
    Double_t r0_fit_phi_min = v_tmp_tmp.Phi();
    v_tmp_tmp = v_r0_fit_right - v_r0_fit;
    Double_t r0_fit_phi_max = v_tmp_tmp.Phi();
    v_tmp_tmp = v_r1_fit_left - v_r1_fit;
    Double_t r1_fit_phi_min = v_tmp_tmp.Phi();
    v_tmp_tmp = v_r1_fit_right - v_r1_fit;
    Double_t r1_fit_phi_max = v_tmp_tmp.Phi();
    //
    Double_t x_fit[nn_fit_circle1+nn_fit_circle2];
    Double_t y_fit[nn_fit_circle1+nn_fit_circle2];
    //
    v_fit.SetMagPhi(TMath::Abs(r0),0.0);
    //phi_fit_min = ;
    //phi_fit_max = ;
    for(Int_t j = 0;j<nn_fit_circle1;j++){
      //phi_fit = (2.0*TMath::Pi() - 0.0)/(nn_fit_circle1 - 1)*i;
      phi_fit = r0_fit_phi_min + (r0_fit_phi_max - r0_fit_phi_min)/(nn_fit_circle1 - 1)*j;
      TVector2 v_tmp = (v_fit.Rotate(phi_fit) + v_r0_fit);
      //if(v_tmp.X()>=x0 && v_tmp.X()<=x1 && v_tmp.Y()>y2){
      x_fit[j] = v_tmp.X();
      y_fit[j] = v_tmp.Y();
      //}
    }
    //
    v_fit.SetMagPhi(TMath::Abs(r1),0.0);
    for(Int_t j = 0;j<nn_fit_circle2;j++){
      phi_fit = (2.0*TMath::Pi() - 0.0)/(nn_fit_circle1 - 1)*j;
      phi_fit = r1_fit_phi_min + (r1_fit_phi_max - r1_fit_phi_min)/(nn_fit_circle2 - 1)*j;
      TVector2 v_tmp = (v_fit.Rotate(phi_fit) + v_r1_fit);
      //if(v_tmp.X()>x1 && v_tmp.X()<=x2 && v_tmp.Y()>y2){
      x_fit[j+nn_fit_circle1] = v_tmp.X();
      y_fit[j+nn_fit_circle2] = v_tmp.Y();
      //}
    }
    //
    TGraph *gr_fit = new TGraph((nn_fit_circle1+nn_fit_circle2),x_fit,y_fit);
    gr_fit->SetLineWidth(3.0);
    gr_fit->SetMarkerStyle(7);
    gr_fit->SetMarkerColor(kBlue+2);
    gr_fit->Fit("pol8","Q","",x0,x2);
    TF1 *f_fit_fit = (TF1*)gr_fit->GetListOfFunctions()->FindObject("pol8");    
    Double_t rmserr = getRMSERR(f_fit,f_fit_fit,x0,x2,nnRMSERR);
    rmserr_arr[ii] = rmserr;
    xA_arr[ii] = xA;
    if(ii>=1){
      if(rmserr<rmserr_arr[ii-1]){
	fit_state->x0 = x0;
	fit_state->x2 = x2;
	fit_state->y0 = y0;
	fit_state->y2 = y2;
	fit_state->xC = xC;
	fit_state->yC = yC;
	fit_state->xA = xA;
	fit_state->yA = yA;
	fit_state->xB = xB;
	fit_state->yB = yB;
	fit_state->yDer0 = yDer0;
	fit_state->yDer2 = yDer2;
	fit_state->m = m;
	fit_state->x1 = x1;
	fit_state->y1 = y1;
	fit_state->q0 = q0;
	fit_state->p0 = p0;
	fit_state->r0 = r0;
	fit_state->q1 = q1;
	fit_state->p1 = p1;
	fit_state->r1 = r1;
	fit_state->rmserr = rmserr;
	fit_state->vec_x_fit.clear();
	fit_state->vec_y_fit.clear();
	for(Int_t jj = 0;jj<(nn_fit_circle1+nn_fit_circle2);jj++){
	  fit_state->vec_x_fit.push_back(x_fit[jj]);
	  fit_state->vec_y_fit.push_back(y_fit[jj]);
	}
      }
      //cout<<"rmserr "<<rmserr<<endl;
      ////////////////// FITED DATA ///////////////////////////////////      
    }
  }
}
  
void plot_gr(TGraph *gr, TCanvas *c1){
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.05);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.05);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->GetXaxis()->SetLimits(xMin,xMax);
  mg->Draw("APL");
  mg->SetMinimum(yMin);
  mg->SetMaximum(yMax);
  
  //c1->SaveAs(outGifname.Data());
}

Double_t getRMSERR(TF1 *f_data,TF1 *f_fit, Double_t xmin, Double_t xmax, Int_t nn){
  Double_t x;
  Double_t ydata;
  Double_t yfit;
  Double_t rmse = 0.0;
  for(Int_t i = 0;i<nn;i++){
    x = xmin + (xmax - xmin)/(nn-1)*i;
    ydata = f_data->Eval(x);
    yfit = f_fit->Eval(x);
    rmse = rmse + (ydata-yfit)*(ydata-yfit);
  }
  return TMath::Sqrt(rmse)/nn;
}

void print_str(fit_state_str *fit_state){
  cout<<"x0     "<<fit_state->x0<<endl
      <<"x2     "<<fit_state->x2<<endl
      <<"y0     "<<fit_state->y0<<endl
      <<"y2     "<<fit_state->y2<<endl;
  cout<<"xC     "<<fit_state->xC<<endl
      <<"yC     "<<fit_state->yC<<endl
      <<"xA     "<<fit_state->xA<<endl
      <<"yA     "<<fit_state->yA<<endl;
  cout<<"xB     "<<fit_state->xB<<endl
      <<"yB     "<<fit_state->yB<<endl
      <<"x1     "<<fit_state->x1<<endl
      <<"y1     "<<fit_state->y1<<endl;
  cout<<"yDer0  "<<fit_state->yDer0<<endl
      <<"yDer2  "<<fit_state->yDer2<<endl
      <<"m      "<<fit_state->m<<endl
      <<"q0     "<<fit_state->q0<<endl;
  cout<<"p0     "<<fit_state->p0<<endl
      <<"r0     "<<fit_state->r0<<endl
      <<"q1     "<<fit_state->q1<<endl
      <<"p1     "<<fit_state->p1<<endl;
  cout<<"r1     "<<fit_state->r1<<endl
      <<"rmserr "<<fit_state->rmserr<<endl;
}
