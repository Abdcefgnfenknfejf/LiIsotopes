#include "PlotTools.hh"
#include "Utilities.hh"

#include <iostream>
#include <map>

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>

void SetGraphStyle(TGraphAsymmErrors* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName,
                   std::string graphName, float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

  assert(graph);

  graph->SetTitle(graphName.c_str());
  graph->GetXaxis()->SetTitle(xAxisName.c_str());
  graph->GetXaxis()->SetTitleFont(22);
  graph->GetXaxis()->SetLabelFont(22);
  graph->GetXaxis()->SetLabelSize(labelSizeX);
  graph->GetXaxis()->SetTitleSize(titleSizeX);
  graph->GetXaxis()->SetTitleOffset(titleOffsetX);
  graph->GetYaxis()->SetTitle(yAxisName.c_str());
  graph->GetYaxis()->SetTitleFont(22);
  graph->GetYaxis()->SetLabelFont(22);
  graph->GetYaxis()->SetLabelSize(labelSizeY);
  graph->GetYaxis()->SetTitleSize(titleSizeY);
  graph->GetYaxis()->SetTitleOffset(titleOffsetY);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(markerColor);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
}



void SetGraphStyle(TGraphErrors* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName,
		   std::string graphName, float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

   assert(graph);

  graph->SetTitle(graphName.c_str());
  graph->GetXaxis()->SetTitle(xAxisName.c_str());
  graph->GetXaxis()->SetTitleFont(22);
  graph->GetXaxis()->SetLabelFont(22);
  graph->GetXaxis()->SetLabelSize(labelSizeX);
  graph->GetXaxis()->SetTitleSize(titleSizeX);
  graph->GetXaxis()->SetTitleOffset(titleOffsetX);
  graph->GetYaxis()->SetTitle(yAxisName.c_str());
  graph->GetYaxis()->SetTitleFont(22);
  graph->GetYaxis()->SetLabelFont(22);
  graph->GetYaxis()->SetLabelSize(labelSizeY);
  graph->GetYaxis()->SetTitleSize(titleSizeY);
  graph->GetYaxis()->SetTitleOffset(titleOffsetY);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(markerColor);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
}


void SetGraphStyle(TGraph* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName, std::string graphName,
		   float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

   assert(graph);

  graph->SetTitle(graphName.c_str());
  graph->GetXaxis()->SetTitle(xAxisName.c_str());
  graph->GetXaxis()->SetTitleFont(22);
  graph->GetXaxis()->SetLabelFont(22);
  graph->GetXaxis()->SetLabelSize(labelSizeX);
  graph->GetXaxis()->SetTitleSize(titleSizeX);
  graph->GetXaxis()->SetTitleOffset(titleOffsetX);
  graph->GetYaxis()->SetTitle(yAxisName.c_str());
  graph->GetYaxis()->SetTitleFont(22);
  graph->GetYaxis()->SetLabelFont(22);
  graph->GetYaxis()->SetLabelSize(labelSizeY);
  graph->GetYaxis()->SetTitleSize(titleSizeY);
  graph->GetYaxis()->SetTitleOffset(titleOffsetY);
  graph->SetMarkerColor(markerColor);
  graph->SetLineColor(markerColor);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
  graph->SetLineWidth(2);
}

void SetHistoStyle(TH1* histo, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName, std::string graphName,
		   float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

  assert(histo);

  histo->SetTitle(graphName.c_str());
  histo->GetXaxis()->SetTitle(xAxisName.c_str());
  histo->GetXaxis()->SetTitleFont(22);
  histo->GetXaxis()->SetLabelFont(22);
  histo->GetXaxis()->SetLabelSize(labelSizeX);
  histo->GetXaxis()->SetTitleSize(titleSizeX);
  histo->GetXaxis()->SetTitleOffset(titleOffsetX);
  histo->GetYaxis()->SetTitle(yAxisName.c_str());
  histo->GetYaxis()->SetTitleFont(22);
  histo->GetYaxis()->SetLabelFont(22);
  histo->GetYaxis()->SetLabelSize(labelSizeY);
  histo->GetYaxis()->SetTitleSize(titleSizeY);
  histo->GetYaxis()->SetTitleOffset(titleOffsetY);
  histo->SetMarkerColor(markerColor);
  histo->SetLineColor(markerColor);
  histo->SetMarkerStyle(markerStyle);
  histo->SetMarkerSize(markerSize);
}

void SetHistoStyle(TH2* histo, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName, std::string graphName,
		   float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

  assert(histo);
  histo->SetTitle(graphName.c_str());
  histo->GetXaxis()->SetTitle(xAxisName.c_str());
  histo->GetXaxis()->SetTitleFont(22);
  histo->GetXaxis()->SetLabelFont(22);
  histo->GetXaxis()->SetLabelSize(labelSizeX);
  histo->GetXaxis()->SetTitleSize(titleSizeX);
  histo->GetXaxis()->SetTitleOffset(titleOffsetX);
  histo->GetYaxis()->SetTitle(yAxisName.c_str());
  histo->GetYaxis()->SetTitleFont(22);
  histo->GetYaxis()->SetLabelFont(22);
  histo->GetYaxis()->SetLabelSize(labelSizeY);
  histo->GetYaxis()->SetTitleSize(titleSizeY);
  histo->GetYaxis()->SetTitleOffset(titleOffsetY);
  histo->SetMarkerColor(markerColor);
  histo->SetLineColor(markerColor);
  histo->SetMarkerStyle(markerStyle);
  histo->SetMarkerSize(markerSize);
}

void CompareHist(TCanvas* myCanvas, std::string Pad1Name, std::string PadPullName, TH1* compareHist1, TH1* compareHist2, TH1* pullHist){
  myCanvas->cd();
  TPad *pad1 = new TPad(Form("%s", Pad1Name.c_str()), Form("%s", Pad1Name.c_str()), 0.0, 0.35, 1, 1.);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad1->Draw();
  TPad *pad2Pull = new TPad(Form("%s", PadPullName.c_str()), Form("%s", PadPullName.c_str()), 0.0, 0 , 1, 0.35);
  pad2Pull->SetTopMargin(0.00001);
  pad2Pull->SetBottomMargin(0.1);
  pad2Pull->SetBorderMode(0);
  pad2Pull->Draw();
  pad1->cd();
  compareHist1->Draw("E");
  compareHist2->Draw("EHISTSAME");
  pad1->Modified();
  pad2Pull->cd();
  pullHist->Draw("EHISTSAME");
  pad2Pull->Modified();

}


void CompareHistWithLegend(TCanvas* myCanvas, std::string Pad1Name, std::string PadPullName, TH1* compareHist1, TH1* compareHist2, std::string lengenForComparisonHist1, std::string lengenForComparisonHist2, TH1* pullHist){
  gStyle->SetErrorX(0);
  myCanvas->cd();
  TLegend* legend = new TLegend(0.15, 0.75 , 0.29, 0.90);
  legend->SetBorderSize(1);

  legend->SetFillStyle(1);
  legend->AddEntry(compareHist1, Form("Compare_%s", lengenForComparisonHist1.c_str()), "lp");
  legend->AddEntry(compareHist2, Form("Compare_%s", lengenForComparisonHist2.c_str()), "lp");
  legend->Draw();
  TPad *pad1 = new TPad(Form("%s", Pad1Name.c_str()), Form("%s", Pad1Name.c_str()), 0.0, 0.35, 1, 1.);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad1->Draw();
  TPad *pad2Pull = new TPad(Form("%s", PadPullName.c_str()), Form("%s", PadPullName.c_str()), 0.0, 0 , 1, 0.35);
  pad2Pull->SetTopMargin(0.00001);
  pad2Pull->SetBottomMargin(0.1);
  pad2Pull->SetBorderMode(0);
  pad2Pull->Draw();
  pad1->cd();
  compareHist1->Draw("E");
  compareHist2->Draw("E SAME");
  pad1->Modified();
  pad2Pull->cd();
  pullHist->Draw("E");
  DrawLineOnCanva(pad2Pull, 0.0, 1.0, 116, 1.0, kRed, 7); //Could be modified.
  pad2Pull->Modified();

}

void DrawLineOnCanva(TCanvas* myCanvas, double xStart, double yStart, double xEnd, double yEnd, Color_t lineColor, Style_t lineStyle){
  myCanvas->cd();
  TLine *line  = new TLine(xStart, yStart, xEnd, yEnd);
  line->SetLineColor(lineColor);
  line->SetLineStyle(lineStyle);
  line->Draw();
}

void DrawLineOnCanva(TVirtualPad* myPad, double xStart, double yStart, double xEnd, double yEnd, Color_t lineColor, Style_t lineStyle){
  myPad->cd();
  TLine *line  = new TLine(xStart, yStart, xEnd, yEnd);
  line->SetLineColor(lineColor);
  line->SetLineStyle(lineStyle);
  line->Draw();
}

void SetCanvasStyle(TCanvas* myCanvas, bool logx, bool logy, bool logz) {

  assert(myCanvas);

  myCanvas->SetBorderSize(2);
  myCanvas->SetTickx();
  myCanvas->SetTicky();
  if (logx)
    myCanvas->SetLogx();
  if (logy)
    myCanvas->SetLogy();
  if (logz)
    myCanvas->SetLogz();
}

void SetCanvasStyle(TVirtualPad* myPad, bool logx, bool logy, bool logz) {

  assert(myPad);

  myPad->SetBorderSize(2);
  if (logx)
    myPad->SetLogx();
  if (logy)
    myPad->SetLogy();
  if (logz)
    myPad->SetLogz();
}

void SetMomentsStyle(Utilities::MomentsGraphs moments, Style_t markerStyle, Color_t markerColor, Size_t markerSize, std::string xAxisName, std::string yAxisName,
		     std::string graphName, float labelSizeX, float titleSizeX, float titleOffsetX, float labelSizeY, float titleSizeY, float titleOffsetY) {

  assert(moments.meanGraph);

  moments.meanGraph->SetTitle(graphName.c_str());
  moments.meanGraph->GetXaxis()->SetTitle(xAxisName.c_str());
  moments.meanGraph->GetXaxis()->SetTitleFont(22);
  moments.meanGraph->GetXaxis()->SetLabelFont(22);
  moments.meanGraph->GetXaxis()->SetLabelSize(labelSizeX);
  moments.meanGraph->GetXaxis()->SetTitleSize(titleSizeX);
  moments.meanGraph->GetXaxis()->SetTitleOffset(titleOffsetX);
  moments.meanGraph->GetYaxis()->SetTitle(yAxisName.c_str());
  moments.meanGraph->GetYaxis()->SetTitleFont(22);
  moments.meanGraph->GetYaxis()->SetLabelFont(22);
  moments.meanGraph->GetYaxis()->SetLabelSize(labelSizeY);
  moments.meanGraph->GetYaxis()->SetTitleSize(titleSizeY);
  moments.meanGraph->GetYaxis()->SetTitleOffset(titleOffsetY);
  moments.meanGraph->SetMarkerColor(markerColor);
  moments.meanGraph->SetLineColor(markerColor);
  moments.meanGraph->SetMarkerStyle(markerStyle);
  moments.meanGraph->SetMarkerSize(markerSize);
  moments.meanGraph->SetLineWidth(2);

  moments.meanGraphErrors->SetTitle(graphName.c_str());
  moments.meanGraphErrors->GetXaxis()->SetTitle(xAxisName.c_str());
  moments.meanGraphErrors->GetXaxis()->SetTitleFont(22);
  moments.meanGraphErrors->GetXaxis()->SetLabelFont(22);
  moments.meanGraphErrors->GetXaxis()->SetLabelSize(labelSizeX);
  moments.meanGraphErrors->GetXaxis()->SetTitleSize(titleSizeX);
  moments.meanGraphErrors->GetXaxis()->SetTitleOffset(titleOffsetX);
  moments.meanGraphErrors->GetYaxis()->SetTitle(yAxisName.c_str());
  moments.meanGraphErrors->GetYaxis()->SetTitleFont(22);
  moments.meanGraphErrors->GetYaxis()->SetLabelFont(22);
  moments.meanGraphErrors->GetYaxis()->SetLabelSize(labelSizeY);
  moments.meanGraphErrors->GetYaxis()->SetTitleSize(titleSizeY);
  moments.meanGraphErrors->GetYaxis()->SetTitleOffset(titleOffsetY);
  moments.meanGraphErrors->SetMarkerColor(markerColor);
  moments.meanGraphErrors->SetLineColor(markerColor);
  moments.meanGraphErrors->SetMarkerStyle(markerStyle);
  moments.meanGraphErrors->SetMarkerSize(markerSize);
  moments.meanGraphErrors->SetLineWidth(2);
  moments.meanMinusRmsGraph->SetMarkerStyle(markerStyle);
  moments.meanMinusRmsGraph->SetMarkerColor(markerColor);
  moments.meanMinusRmsGraph->SetMarkerSize(markerSize);
  moments.meanPlusRmsGraph->SetMarkerStyle(markerStyle);
  moments.meanPlusRmsGraph->SetMarkerColor(markerColor);
  moments.meanPlusRmsGraph->SetMarkerSize(markerSize);
}

