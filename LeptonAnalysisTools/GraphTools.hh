#ifndef GraphTools_hh
#define GraphTools_hh

#include "BinningDefinition.hh"

class TGraph;
class TGraphAsymmErrors;

TGraphAsymmErrors* CorrectBinCentersLaffertyWyatt(TGraphAsymmErrors* spectrum, const Binning::Definition& binning, double gamma);
TGraphAsymmErrors* MultiplyByRigidityOrEnergy(TGraphAsymmErrors* spectrum, double power);

template<typename T>
void RemoveGraphPointsBelowEnergy(T* graph, double xMin) {

  int points = graph->GetN();
  double x = 0.0;
  double y = 0.0;
  graph->GetPoint(0, x, y);
  while (x < xMin && points > 0) {
    graph->RemovePoint(0);
    points = graph->GetN();
    graph->GetPoint(0, x, y);
  }
}


template<typename T>
void RemoveGraphPointsAboveEnergy(T* graph, double xMax) {

  int points = graph->GetN();
  double x = 0.0f;
  double y = 0.0f;
  graph->GetPoint(points - 1, x, y);
  while (x > xMax) {
    graph->RemovePoint(points - 1);
    points = graph->GetN();
    graph->GetPoint(points - 1, x, y);
  }  
}

template<typename T>
void TruncateGraphEnergyRange(T* graph, double xMin, double xMax) {

  if (xMax <= xMin)
    return;
  RemoveGraphPointsBelowEnergy(graph, xMin);
  RemoveGraphPointsAboveEnergy(graph, xMax);
}

template<typename T>
void RemoveGraphPointsInRange(T* graph, double xMin, double xMax) {

  if (xMax <= xMin)
    return;
  int points = graph->GetN();
  double x = 0.0;
  double y = 0.0;
  graph->GetPoint(points - 1, x, y);
  while (points > 0) {
    if (xMin < x && x < xMax) {
      graph->RemovePoint(points - 1);
      points = graph->GetN();
    } else {
      --points;
    }
    graph->GetPoint(points - 1, x, y);
  }
}

template<typename T>
void AdjustGraphYAxisRange(T* graph, double lowFactor = 0.85, double highFactor = 1.15) {

  double yMin = 1e9;
  double yMax = -1e9;
  for (int i = 0; i < graph->GetN(); ++i) {
    double x = 0.0;
    double y = 0.0;
    graph->GetPoint(i, x, y);
    if (y < yMin)
      yMin = y;
    if (y > yMax)
      yMax = y;
  }
  graph->SetMinimum(yMin * lowFactor);
  graph->SetMaximum(yMax * highFactor);
}

void RemoveHorizontalErrorsFromGraph(TGraphAsymmErrors* graph);
void DivideFluxUncertaintyByBinWidth(TGraph* fluxUncertainty);
void DivideFluxByBinWidth(TGraphAsymmErrors* flux);

#endif
