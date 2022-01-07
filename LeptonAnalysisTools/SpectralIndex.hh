#ifndef SpectralIndex_hh
#define SpectralIndex_hh

#include <functional>
#include <vector>
#include "ParticleId.hh"

class TF1;
class TGraph;
class TGraphErrors;

namespace Modelling {
class HistogramDataset;
}

namespace Binning {
class Definition;
}

// All methods return -gamma.
TGraphErrors* CalculateSpectralIndexGraphUsingSlidingWindowFit(Modelling::HistogramDataset* ds, std::function<int(double)> nPointsEachSide, int graphColor);
TGraphErrors* CalculateSpectralIndexGraphUsingSlidingWindowFit(Modelling::HistogramDataset* ds, const std::vector<double>& energyIntervals, int graphColor);

TGraph* CalculateSpectralIndexGraphUsingModulatedPowerLawFit(Modelling::HistogramDataset* ds, ParticleId::Species particleType, const Binning::Definition& binning, double startFitEnergy, double stopFitEnergy);

#endif
