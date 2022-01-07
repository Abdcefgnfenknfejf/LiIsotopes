#ifndef SplineFitTools_hh
#define SplineFitTools_hh


#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSpline.h>

template<typename T>
class TGraphAsTF1 {
public:
  TGraphAsTF1(T* graph)
    : fGraph(graph) {

  }

  double operator() (double* x, double *) {

    return fGraph->Eval(x[0]);
  }

private:
  T* fGraph;
};

TSpline5* FitTGraphAsymmErrorsUsingSplineReturningSpline(TGraphAsymmErrors* graph) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();
  const double xStart = binning.Min();
  const double xEnd = binning.Max();
  const double xDelta = xEnd - xStart;
  const double epsilon = xDelta * 1e-5;

  TGraphAsTF1<TGraphAsymmErrors>* graphAsTF1 = new TGraphAsTF1<TGraphAsymmErrors>(graph);
  TF1* fitFunction = new TF1("splineFitFunction", graphAsTF1, xStart + epsilon, xEnd - epsilon, 0);
  TSpline5* spline5 = new TSpline5("FifthOrderSpline", const_cast<double*>(binning.Bins().data()), fitFunction, int(binning.NumberOfBins() + 1), "b1e1b2e2",
                                   fitFunction->Derivative(xStart), fitFunction->Derivative(xEnd),
                                   (fitFunction->Derivative(xStart + epsilon) - fitFunction->Derivative(xStart)) / epsilon,
                                   (fitFunction->Derivative(xEnd) - fitFunction->Derivative(xEnd - epsilon)) / epsilon);
  return spline5;
}

TGraphAsymmErrors* FitTGraphAsymmErrorsUsingSpline(TGraphAsymmErrors* graph) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();
  const double xStart = binning.Min();
  const double xEnd = binning.Max();
  const double xDelta = xEnd - xStart;
  const double epsilon = xDelta * 1e-5;

  TGraphAsTF1<TGraphAsymmErrors>* graphAsTF1 = new TGraphAsTF1<TGraphAsymmErrors>(graph);
  TF1* fitFunction = new TF1("splineFitFunction", graphAsTF1, xStart + epsilon, xEnd - epsilon, 0);
  TSpline5* spline5 = new TSpline5("FifthOrderSpline", const_cast<double*>(binning.Bins().data()), fitFunction, int(binning.NumberOfBins() + 1), "b1e1b2e2",
                                   fitFunction->Derivative(xStart), fitFunction->Derivative(xEnd),
                                   (fitFunction->Derivative(xStart + epsilon) - fitFunction->Derivative(xStart)) / epsilon,
                                   (fitFunction->Derivative(xEnd) - fitFunction->Derivative(xEnd - epsilon)) / epsilon);

  TGraphAsymmErrors* graphWithErrors = new TGraphAsymmErrors(graph->GetN());
  for (int i = 0; i < graph->GetN(); ++i) {
    double x = graph->GetX()[i];
    double binContent = spline5->Eval(x);
    graphWithErrors->SetPoint(i, x, binContent);
    graphWithErrors->SetPointError(i, 0, 0, graph->GetErrorYlow(i), graph->GetErrorYhigh(i));
  }

  return graphWithErrors;
}

TGraphErrors* FitTGraphErrorsUsingSpline(TGraphErrors* graph) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();
  const double xStart = binning.Min();
  const double xEnd = binning.Max();
  const double xDelta = xEnd - xStart;
  const double epsilon = xDelta * 1e-5;

  TGraphAsTF1<TGraphErrors>* graphAsTF1 = new TGraphAsTF1<TGraphErrors>(graph);
  TF1* fitFunction = new TF1("splineFitFunction", graphAsTF1, xStart + epsilon, xEnd - epsilon, 0);
  TSpline5* spline5 = new TSpline5("FifthOrderSpline", const_cast<double*>(binning.Bins().data()), fitFunction, int(binning.NumberOfBins() + 1), "b1e1b2e2",
                                   fitFunction->Derivative(xStart), fitFunction->Derivative(xEnd),
                                   (fitFunction->Derivative(xStart + epsilon) - fitFunction->Derivative(xStart)) / epsilon,
                                   (fitFunction->Derivative(xEnd) - fitFunction->Derivative(xEnd - epsilon)) / epsilon);

  TGraphErrors* graphWithErrors = new TGraphErrors(graph->GetN());
  for (int i = 0; i < graph->GetN(); ++i) {
    double x = graph->GetX()[i];
    double binContent = spline5->Eval(x);
    graphWithErrors->SetPoint(i, x, binContent);
    graphWithErrors->SetPointError(i, 0, graph->GetErrorY(i));
  }

  return graphWithErrors;
}

#endif
