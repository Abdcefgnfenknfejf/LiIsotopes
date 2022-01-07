#ifndef SmoothingTools_hh
#define SmoothingTools_hh

#include <cmath>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TGraph.h>

#include "GraphTools.hh"

template<typename T>
T* SmoothGraphLowess(T* graph, double xFrom = -1.0, double xTo = -1.0, double span = 0.67, int iterations = 3) {

  TGraphSmooth* smoother = new TGraphSmooth;
  TGraph* rawSmoothedGraph = smoother->SmoothLowess(graph, "", span, iterations);

  // Construct final graph from smoothed part and non-smoothed rest.
  int size = graph->GetN();
  T* smoothedGraph = new T(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    double ySmoothed = 0.0;
    if (xFrom == xTo || (x > xFrom && x < xTo))
      ySmoothed = rawSmoothedGraph->Eval(x);
    else
      ySmoothed = y;

    smoothedGraph->SetPoint(point, x, ySmoothed);
  }

  return smoothedGraph;
}

template<>
TGraph* SmoothGraphLowess(TGraph* graph, double xFrom, double xTo, double span, int iterations) {

  TGraphSmooth* smoother = new TGraphSmooth;
  TGraph* rawSmoothedGraph = smoother->SmoothLowess(graph, "", span, iterations);

  // Construct final graph from smoothed part and non-smoothed rest.
  int size = graph->GetN();
  TGraph* smoothedGraph = new TGraph(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    double ySmoothed = 0.0;
    if (xFrom == xTo || (x > xFrom && x < xTo))
      ySmoothed = rawSmoothedGraph->Eval(x);
    else
      ySmoothed = y;

    smoothedGraph->SetPoint(point, x, ySmoothed);
  }

  return smoothedGraph;
}

template<typename T>
T* SmoothGraph(T* graph, double xFrom = -1.0, double xTo = -1.0, double bass = 0.0, double span = 0.0) {

  TGraphSmooth* smoother = new TGraphSmooth;
  TGraph* rawSmoothedGraph = smoother->SmoothSuper(graph, "", bass, span);

  // Construct final graph from smoothed part and non-smoothed rest.
  int size = graph->GetN();
  T* smoothedGraph = new T(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    double ySmoothed = 0.0;
    if (xFrom == xTo || (x > xFrom && x < xTo))
      ySmoothed = rawSmoothedGraph->Eval(x);
    else
      ySmoothed = y;

    smoothedGraph->SetPoint(point, x, ySmoothed);
  }

  return smoothedGraph;
}

template<>
TGraph* SmoothGraph(TGraph* graph, double xFrom, double xTo, double bass, double span) {

  TGraphSmooth* smoother = new TGraphSmooth;
  TGraph* rawSmoothedGraph = smoother->SmoothSuper(graph, "", bass, span);

  // Construct final graph from smoothed part and non-smoothed rest.
  int size = graph->GetN();
  TGraph* smoothedGraph = new TGraph(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    double ySmoothed = 0.0;
    if (xFrom == xTo || (x > xFrom && x < xTo))
      ySmoothed = rawSmoothedGraph->Eval(x);
    else
      ySmoothed = y;

    smoothedGraph->SetPoint(point, x, ySmoothed);
  }

  return smoothedGraph;
}

TGraph* SmoothLogarithmicGraph(TGraph* graph, double xFrom = -1.0, double xTo = -1.0, double bass = 0.0, double span = 0.0) {

  static double sLogOffset = 12.0; // Transform log values to assure that y-values are positive after transformation.

  int size = graph->GetN();
  TGraph* logGraph = new TGraph(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    graph->GetPoint(point, x, y);
    if (y == 0)
      y = 1;
    logGraph->SetPoint(point, x, std::log10(y) + sLogOffset);
  }

  TGraph* logGraphSmoothed = SmoothGraph(logGraph, xFrom, xTo, bass, span);
  TGraph* graphSmoothed = new TGraph(size);
  for (int point = 0; point < size; ++point) {
    double x, y;
    logGraphSmoothed->GetPoint(point, x, y);
    y = std::pow(10, y - sLogOffset);
    if (y > 1)
      y = 1;
    graphSmoothed->SetPoint(point, x, y);
  }

  delete logGraph;
  delete logGraphSmoothed;
  return graphSmoothed;
}

class FitRange {
public:
  FitRange(double _start, double _stop, double _smoothingParameter = 0.67)
    : start(_start)
    , stop(_stop)
    , smoothingParameter(_smoothingParameter) {

  }

  double start { 0.0 };
  double stop { 0.0 };
  double smoothingParameter { 0.0 };
};

TGraph* SmoothGraphInMultipleRanges(TGraph* graph, const std::vector<FitRange>& ranges) {

  TGraph* graphSmoothed = new TGraph(graph->GetN());
  std::vector<TGraph*> parts;
  for (const auto& range : ranges) {
    if (range.stop <= range.start)
      continue;

    TGraph* part = static_cast<TGraph*>(graph->Clone());
    TruncateGraphEnergyRange(part, range.start, range.stop);
    parts.emplace_back(SmoothGraphLowess(part, range.start, range.stop, range.smoothingParameter));
  }

  assert(parts.size() > 0);

  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);

    if (x <= ranges.front().start)
      y = graph->Eval(x);
    else {
      for (TGraph* part : parts) {
        if (x <= part->GetX()[part->GetN() - 1]) {
          y = part->Eval(x);
          break;
        }
      }
    }

    graphSmoothed->SetPoint(i, x, y);
  }

  return graphSmoothed;
}

#endif // SmoothingTools_hh
