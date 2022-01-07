#ifndef AveragingTools_hh
#define AveragingTools_hh

#include <vector>

#include <TGraphAsymmErrors.h>
#include <TMath.h>

void CalculateAverageFlux(std::vector<TGraphAsymmErrors*>& graphs, TGraphAsymmErrors* meanGraph, TGraphAsymmErrors* rmsGraph, TGraphAsymmErrors* meanUncertaintyGraph) {

  std::vector<double> yData, yDataFiltered, yUncertaintyLow, yUncertaintyHigh, yUncertaintyLowFiltered, yUncertaintyHighFiltered;
  for (int point = 0; point < graphs[0]->GetN(); ++point) {
    yData.clear();
    yUncertaintyLow.clear();
    yUncertaintyHigh.clear();
    unsigned int jMin = 0;
    unsigned int jMax = 0;
    double yMin = 1e99;
    double yMax = -1e99;
    for (unsigned int graphIndex = 0; graphIndex < graphs.size(); ++graphIndex) {
      double y = graphs[graphIndex]->GetY()[point];
      yData.push_back(y);
      yUncertaintyLow.push_back(graphs[graphIndex]->GetErrorYlow(point));
      yUncertaintyHigh.push_back(graphs[graphIndex]->GetErrorYhigh(point));

      if (y < yMin) {
        yMin = y;
        jMin = graphIndex;
      }

      if (y > yMax) {
        yMax = y;
        jMax = graphIndex;
      }
    }
    double median = TMath::Median(yData.size(), yData.data());
    yDataFiltered.clear();
    yUncertaintyLowFiltered.clear();
    yUncertaintyHighFiltered.clear();
    if (yData.size() < 5) {
      for (unsigned int i = 0; i < yData.size(); ++i) {
        if (std::abs(yData[i]) < 1e-6 && std::abs(yUncertaintyLow[i] - 1.0) < 1e-6 && std::abs(yUncertaintyHigh[i] - 1.0) < 1e-6)
          continue;
        yDataFiltered.push_back(yData[i]);
        yUncertaintyLowFiltered.push_back(yUncertaintyLow[i]);
        yUncertaintyHighFiltered.push_back(yUncertaintyHigh[i]);
      }
    } else {
      for (unsigned int i = 0; i < yData.size(); ++i) {
        // Remove outliers
        if (i == jMin || i == jMax)
          continue;
        if (std::abs(yData[i]) < 1e-6 && std::abs(yUncertaintyLow[i] - 1.0) < 1e-6 && std::abs(yUncertaintyHigh[i] - 1.0) < 1e-6)
          continue;
        yDataFiltered.push_back(yData[i]);
        yUncertaintyLowFiltered.push_back(yUncertaintyLow[i]);
        yUncertaintyHighFiltered.push_back(yUncertaintyHigh[i]);
      }
    }

    double meanUncertaintyLow = TMath::Mean(yUncertaintyLowFiltered.size(), yUncertaintyLowFiltered.data());
    double meanUncertaintyHigh = TMath::Mean(yUncertaintyHighFiltered.size(), yUncertaintyHighFiltered.data());
    double rms = TMath::RMS(yDataFiltered.size(), yDataFiltered.data());

    double x = graphs[0]->GetX()[point];
    meanGraph->SetPoint(point, x, median);
    meanGraph->SetPointError(point, 0, 0, rms, rms);
    rmsGraph->SetPoint(point, x, rms);

    // FIXME: Don't flatten the error here into a symmetric one.
    meanUncertaintyGraph->SetPoint(point,x, TMath::Sqrt(0.5 * (meanUncertaintyLow * meanUncertaintyLow + meanUncertaintyHigh * meanUncertaintyHigh)));
  }
}

#endif
