#include "TemplateParametrisation.hh"

#include "AnalysisTools.hh"
#include "BinningTools.hh"
#include "ValueWithError.hh"
#include "Utilities.hh"
#include "LithiumBinning.hh"

#include <TH1D.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TAttLine.h>

#define INFO_OUT_TAG "TemplateParametrisation"
#include "debugging.hh"

double function_MassTemplateWidth(double* xEkinN, double* par){
  double mass = par[0];
  double charge = 3.0;
  double massPerA = 0.9338; //GeV/n                                                                                                                                                            
  double massInGeV = mass * massPerA;
  double R = massInGeV/charge * TMath::Sqrt(TMath::Power(xEkinN[0]/massPerA + 1, 2) -1);
  double gamma = xEkinN[0]/massPerA +1;
  double pdf = TMath::Sqrt(TMath::Power(par[1] + par[2]*R + par[3]/R, 2) + par[4] * par[4] * TMath::Power(gamma, 4));
  return pdf;
}


double MyGaus(double* x, double* par){
  double g1  = 0.0;
  double arg1 = (par[2] != 0.0) ? (x[0] - par[1])/(par[2]) : 0.0;
  g1 = par[0] * exp(-0.5*arg1*arg1);

  return g1;
}


double AsyGaus(double* x, double* par){
  double PDF = 0.0;
  double g1  = 0.0;
  double g2  = 0.0;

  double arg1 = (par[2] != 0.0) ? (x[0] - par[1])/(par[2]/par[3]) : 0.0;
  double arg2 = (par[2] != 0.0) ? (x[0] - par[1])/(par[2]) : 0.0;

  g1 = par[0] * exp(-0.5*arg1*arg1);
  g2 = par[0] * exp(-0.5*arg2*arg2);

  if (x[0] <= par[1])
    PDF = g1;
  else
    PDF = g2;

  return PDF;
}


double GausPlusAsyGausFunction(double *x, double* par){

  double PDF = 0.0;
  double gaus_One  = 0.0;
  double AsyGaus_left  = 0.0;
  double AsyGaus_right  = 0.0;

  double norm_One = par[0];
  double norm_Two =  par[0] * par[1];
  double mean = par[2];
  double sigma_One = par[3];
  double sigma_Two = par[4];
  double AsyFactor = par[5];

  double arg1 = (sigma_One != 0.0) ? (x[0] - mean)/(sigma_One) : 0.0;
  gaus_One = exp(-0.5*arg1*arg1);

  double arg21 = (sigma_Two != 0.0) ? (x[0] - mean)/(sigma_Two/AsyFactor) : 0.0;
  double arg22 = (sigma_Two != 0.0) ? (x[0] - mean)/(sigma_Two) : 0.0;
  AsyGaus_left = exp(-0.5* arg21 * arg21);
  AsyGaus_right = exp(-0.5* arg22 * arg22);

  if (x[0] <= mean){
    PDF = norm_One * gaus_One + norm_Two * AsyGaus_left;
  }
  else {
    PDF = norm_One * gaus_One + norm_Two * AsyGaus_right;
  }
 return PDF;
}

double DoubleGausFunction(double *x, double* par){
  double PDF = 0.0;
  double g1  = 0.0;
  double g2  = 0.0;

  double mean = par[2];

  double arg1 = (par[3] != 0.0) ? (x[0] - mean)/(par[3]) : 0.0;
  g1 = exp(-0.5*arg1*arg1);

  double arg2 = (par[4] != 0.0) ? (x[0] - mean)/(par[4]) : 0.0;
  g2 = exp(-0.5* arg2 * arg2);
  PDF = par[0] * g1 + par[0] *par[1] * g2;

  return PDF;
}

double TwoIsotopesMassFunction_FixFactor(double *x, double* par){

  double norm_P_Li7 = par[0];
  double norm_S_Li7 = par[0] * par[1];
  double mean_Li7 = par[2];
  double sigma_Li7 = par[3];
  double sigmaSecond_Li7 = par[4];
  double AsyFactor_Li7 = par[5];
  double ratio = par[6];

  double arg1_Li7 = (sigma_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigma_Li7) : 0.0;
  double gausOne_Li7 = exp(-0.5*arg1_Li7*arg1_Li7);
  double arg21_Li7 = (sigmaSecond_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigmaSecond_Li7/AsyFactor_Li7) : 0.0;
  double arg22_Li7 = (sigmaSecond_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigmaSecond_Li7) : 0.0;
  double AsyGaus_left_Li7 = exp(-0.5* arg21_Li7 * arg21_Li7);
  double AsyGaus_right_Li7 = exp(-0.5* arg22_Li7 * arg22_Li7);

  double factor = 6.0/7.0;   //A small factor need to be added for R6 != R7 for each EkinN bin
  double norm_P_Li6 = par[0] * factor;
  double norm_S_Li6 = par[0] * par[1] * factor;
  double mean_Li6 = mean_Li7/factor;
  double sigma_Li6 = sigma_Li7/factor;
  double sigmaSecondary_Li6 = sigmaSecond_Li7/factor;
  double AsyFactor_Li6 = par[5];
  double arg1_Li6 = (sigma_Li6 != 0.0) ? (x[0] - mean_Li6)/(sigma_Li6) : 0.0;
  double gausOne_Li6 = exp(-0.5*arg1_Li6*arg1_Li6);
  double arg21_Li6 = (sigmaSecondary_Li6!= 0.0) ? (x[0] - mean_Li6)/(sigmaSecondary_Li6/AsyFactor_Li6) : 0.0;
  double arg22_Li6 = (sigmaSecondary_Li6!= 0.0) ? (x[0] - mean_Li6)/(sigmaSecondary_Li6) : 0.0;
  double AsyGaus_left_Li6 = exp(-0.5* arg21_Li6 * arg21_Li6);
  double AsyGaus_right_Li6 = exp(-0.5* arg22_Li6 * arg22_Li6);

  double pdf1_Li6 = norm_P_Li6 * gausOne_Li6 + norm_S_Li6 * AsyGaus_left_Li6;
  double pdf2_Li6 = norm_P_Li6 * gausOne_Li6 + norm_S_Li6 * AsyGaus_right_Li6;

  double pdf1_Li7 = norm_P_Li7 * gausOne_Li7 + norm_S_Li7 * AsyGaus_left_Li7;
  double pdf2_Li7 = norm_P_Li7 * gausOne_Li7 + norm_S_Li7 * AsyGaus_right_Li7;

  double PDF =0.0;

  if (x[0] <= mean_Li7){
    PDF = ratio * pdf1_Li6 + pdf1_Li7;
  }
  else if ( mean_Li7 < x[0] < mean_Li6) {
    PDF = ratio * pdf1_Li6 + pdf2_Li7;
  }
  else {
    PDF = ratio * pdf2_Li6 + pdf2_Li7;
  }
  return PDF;
}


//compared with the last function with a fixed factor of 6.0/7.0, this function has par[7] to par[12],
//describe the difference between the parameteres of Li6 and Li7, these information are from MC
double TwoIsotopesMassFunction(double *x, double* par){

  double norm_P_Li7 = par[0];
  double norm_S_Li7 = par[0] * par[1];
  double mean_Li7 = par[2];
  double sigma_Li7 = par[3];
  double sigmaSecond_Li7 = par[4];
  double AsyFactor_Li7 = par[5];
  double ratio = par[6];
  
  double norm_P_Li6 = par[0] * par[7];
  double norm_S_Li6 = par[0] * par[1] * par[8];
  double mean_Li6 = mean_Li7 * par[9];
  double sigma_Li6 = sigma_Li7 * par[10];
  double sigmaSecondary_Li6 = sigmaSecond_Li7 * par[11];
  double AsyFactor_Li6 = par[5] * par[12];

  double arg1_Li7 = (sigma_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigma_Li7) : 0.0;
  double gausOne_Li7 = exp(-0.5*arg1_Li7*arg1_Li7);
  double arg21_Li7 = (sigmaSecond_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigmaSecond_Li7/AsyFactor_Li7) : 0.0;
  double arg22_Li7 = (sigmaSecond_Li7 != 0.0) ? (x[0] - mean_Li7)/(sigmaSecond_Li7) : 0.0;
  double AsyGaus_left_Li7 = exp(-0.5* arg21_Li7 * arg21_Li7);
  double AsyGaus_right_Li7 = exp(-0.5* arg22_Li7 * arg22_Li7);

  double arg1_Li6 = (sigma_Li6 != 0.0) ? (x[0] - mean_Li6)/(sigma_Li6) : 0.0;
  double gausOne_Li6 = exp(-0.5*arg1_Li6*arg1_Li6);
  double arg21_Li6 = (sigmaSecondary_Li6!= 0.0) ? (x[0] - mean_Li6)/(sigmaSecondary_Li6/AsyFactor_Li6) : 0.0;
  double arg22_Li6 = (sigmaSecondary_Li6!= 0.0) ? (x[0] - mean_Li6)/(sigmaSecondary_Li6) : 0.0;
  double AsyGaus_left_Li6 = exp(-0.5* arg21_Li6 * arg21_Li6);
  double AsyGaus_right_Li6 = exp(-0.5* arg22_Li6 * arg22_Li6);

  double pdf1_Li6 = norm_P_Li6 * gausOne_Li6 + norm_S_Li6 * AsyGaus_left_Li6;
  double pdf2_Li6 = norm_P_Li6 * gausOne_Li6 + norm_S_Li6 * AsyGaus_right_Li6;

  double pdf1_Li7 = norm_P_Li7 * gausOne_Li7 + norm_S_Li7 * AsyGaus_left_Li7;
  double pdf2_Li7 = norm_P_Li7 * gausOne_Li7 + norm_S_Li7 * AsyGaus_right_Li7;

  double PDF =0.0;

  if (x[0] <= mean_Li7){
    PDF = ratio * pdf1_Li6 + pdf1_Li7;
  }
  else if ( mean_Li7 < x[0] < mean_Li6) {
    PDF = ratio * pdf1_Li6 + pdf2_Li7;
  }
  else {
    PDF = ratio * pdf2_Li6 + pdf2_Li7;
  }
  return PDF;
}



void SetGraphStyle(TGraph* graphs, Style_t markerStyle, Color_t markerColor, Size_t markerSize){
  graphs->SetMarkerColor(markerColor);
  graphs->SetMarkerSize(markerSize);
  graphs->SetMarkerStyle(markerStyle);
  graphs->SetLineColor(markerColor);
}

void SetFunctionStyle(TF1* Fuction, Style_t lineStyle, Color_t lineColor, Size_t lineWidth){
  Fuction->SetLineColor(lineColor);
  Fuction->SetLineStyle(lineStyle);
  Fuction->SetLineWidth(lineWidth);
}

void SetHistStyle(TH1* Hist, Style_t MarkerStyle, Color_t MarkerColor, Color_t lineColor, Size_t MarkerSize){
  Hist->SetLineColor(lineColor);
  Hist->SetMarkerColor(MarkerColor);
  Hist->SetMarkerStyle(MarkerStyle);
  Hist->SetMarkerSize(MarkerSize);
  Hist->SetLineColor(lineColor);

  Hist->GetXaxis()->SetRange(25,126);
  Hist->GetXaxis()->SetLabelFont(42);
  Hist->GetXaxis()->SetTitleFont(42);
  Hist->GetYaxis()->SetLabelFont(42);
  Hist->GetYaxis()->SetTitleSize(0.06);
  Hist->GetYaxis()->SetTitleOffset(0.67);
  Hist->GetYaxis()->SetTitleFont(42);
}

void SetGraphAxis(TMultiGraph* multigraph_SigmaPrimary){
  multigraph_SigmaPrimary->GetXaxis()->SetLabelFont(22);
  multigraph_SigmaPrimary->GetXaxis()->SetLabelSize(0.08);
  multigraph_SigmaPrimary->GetXaxis()->SetTitleOffset(1);
  multigraph_SigmaPrimary->GetXaxis()->SetTitleFont(42);
  multigraph_SigmaPrimary->GetYaxis()->SetLabelFont(42);
  multigraph_SigmaPrimary->GetYaxis()->SetLabelSize(0.08);
  multigraph_SigmaPrimary->GetYaxis()->SetTitleSize(0.1);
  multigraph_SigmaPrimary->GetYaxis()->SetTitleOffset(0.4);
  multigraph_SigmaPrimary->GetYaxis()->SetTitleFont(42);
}

void SetResidualPlot(TH1* hist, TH1* pullHistogram, TF1* Function){
  pullHistogram->Reset();

  for (int iBin = 0; iBin < hist->GetNbinsX(); iBin++) {
    pullHistogram->SetBinContent(iBin, (Function->Eval(hist->GetBinCenter(iBin)) - hist->GetBinContent(iBin)) / (hist->GetBinError(iBin) != 0 ? hist->GetBinError(iBin) : 1.0));
  }
  pullHistogram->GetXaxis()->SetLabelSize(0.09);
  pullHistogram->GetXaxis()->SetTitleSize(0.10);
  pullHistogram->GetYaxis()->SetRangeUser(-4, 4);
  pullHistogram->GetYaxis()->SetLabelSize(0.09);
  pullHistogram->GetYaxis()->SetTitleSize(0.30);

  pullHistogram->SetBit(TH1::kNoTitle);
   pullHistogram->GetXaxis()->SetTitle("1/mass [amu]");
   pullHistogram->GetXaxis()->SetRange(25,126);
   pullHistogram->GetXaxis()->SetLabelFont(42);
   pullHistogram->GetXaxis()->SetLabelSize(0.13);
   pullHistogram->GetXaxis()->SetTitleSize(0.12);
   pullHistogram->GetXaxis()->SetTitleOffset(1);
   pullHistogram->GetXaxis()->SetTitleFont(42);
   pullHistogram->GetYaxis()->SetTitle("(model - data) / #sigma");
   pullHistogram->GetYaxis()->SetLabelFont(42);
   pullHistogram->GetYaxis()->SetLabelSize(0.09);
   pullHistogram->GetYaxis()->SetTitleSize(0.11);
   pullHistogram->GetYaxis()->SetTitleOffset(0.31);
   pullHistogram->GetYaxis()->SetTitleFont(42);

   pullHistogram->SetMarkerStyle(kFullCircle);
   pullHistogram->SetMarkerSize(0.6);
   pullHistogram->SetLineWidth(2);
   pullHistogram->SetLineColor(kBlack);
   pullHistogram->SetMarkerColor(kBlack);
}


 void format_line(TAttLine* line,int col,int sty){
    line->SetLineWidth(5);
    line->SetLineColor(col);
    line->SetLineStyle(sty);
 }
