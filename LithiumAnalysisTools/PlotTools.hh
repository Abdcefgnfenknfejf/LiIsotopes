#ifndef PlotTools_hh
#define PlotTools_hh

#include "Utilities.hh"

#include <map>
#include <vector>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TVirtualPad.h>

void SetGraphStyle(TGraphAsymmErrors* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		   std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		   float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void SetGraphStyle(TGraphErrors* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		   std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		   float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void SetGraphStyle(TGraph* graph, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		   std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		   float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void SetHistoStyle(TH1* histo, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		   std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		   float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void SetHistoStyle(TH2* histo, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		   std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		   float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void SetCanvasStyle(TCanvas* canvas, bool logx = false, bool logy = false, bool logz = false);

void SetCanvasStyle(TVirtualPad* pad, bool logx = false, bool logy = false, bool logz = false);

void SetMomentsStyle(Utilities::MomentsGraphs moments, Style_t markerStyle, Color_t markerColor, Size_t markerSize = 1,
		     std::string xAxisName = "", std::string yAxisName = "", std::string graphName = "", float labelSizeX = 0.04, float titleSizeX = 0.04,
		     float titleOffsteX = 1, float labelSizeY = 0.04, float titleSizeY = 0.04, float titleOffsetY = 1);

void CompareHist(TCanvas* myCanvas, std::string Pad1Name, std::string PadPullName, TH1* compareHist1, TH1* compareHist2, TH1* pullHist);
void CompareHistWithLegend(TCanvas* myCanvas, std::string Pad1Name, std::string PadPullName, TH1* compareHist1, TH1* compareHist2, std::string lengenForComparisonHist1, std::string lengenForComparisonHist2, TH1* pullHist);

void DrawLineOnCanva(TCanvas* myCanvas, double xStart, double yStart, double xEnd, double yEnd, Color_t lineColor, Style_t lineStyle);
void DrawLineOnCanva(TVirtualPad* myPad, double xStart, double yStart, double xEnd, double yEnd, Color_t lineColor, Style_t lineStyle);

#endif
