#include "MigrationRandomGenerator.hh"

//#define DEBUG_LEVEL 0
#define INFO_OUT_TAG "MigrationRandomGenerator"
#include <debugging.hh>

#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "MigrationParameterization.hh"

MigrationRandomGenerator::MigrationRandomGenerator(const Binning::Definition& binning, const MigrationParameterization* param) :
  fEnergyBinning(binning),
  fParam(param) {

  // cache integrals of parameterization
  for (unsigned int ibin = 1; ibin <= fEnergyBinning.NumberOfBins(); ++ibin) {

    double Etrue = fEnergyBinning.BinCenterLog(ibin);

    TF1* fc = param->MakeFunction(Etrue);

    unsigned int N = fEnergyBinning.NumberOfBins();
    fIntegralMap[ibin] = std::vector<double>(N + 3, 0.0);
    auto& v = fIntegralMap[ibin];
    v[1] = param->UnderflowProbability(Etrue); // pUnderflow
    for (unsigned int jbin = 1; jbin <= N; ++jbin) {

      double Erec = fEnergyBinning.BinCenterLog(jbin);

      // 0: 0.0
      // 1: underflow
      if (param->IsPdf())
        v[jbin+1] = v[jbin] + fc->Integral(fEnergyBinning.LowEdge(jbin), fEnergyBinning.UpEdge(jbin));
      else
        v[jbin+1] = v[jbin] + fc->Eval(Erec);
      DEBUG_OUT << "Etrue= " << Etrue << " jbin= " << jbin << ", cdf= " << v[jbin] << std::endl;
    }

    v[N+2] = v[N+1] + param->OverflowProbability(Etrue); // pOverflow

    delete fc;

    // normalize
    double norm = v.back();
    for (unsigned int i = 0; i < v.size(); ++i) {
      v[i] /= norm;
    }
  }
}


double MigrationRandomGenerator::GetRandomSmearedEnergy(double Etruth) const {

  // FIXME smear within bin

  if (!fEnergyBinning.IsInRange(Etruth))
    return 0.0;

  const auto& cdf_vector = fIntegralMap.find(fEnergyBinning.FindBin(Etruth))->second;

  if (!fRandom)
    fRandom = new TRandom3(0);

  Double_t r = fRandom->Rndm();
  auto ibin = TMath::BinarySearch(cdf_vector.size(), cdf_vector.data(), r);

  if (ibin < 1){
    DEBUG_OUT_L(2) << "E= " << Etruth << " ran " << r << " -> ibin " << ibin << std::endl;
    return fEnergyBinning.Min() - 0.0001;
  }

  if (ibin > int(fEnergyBinning.NumberOfBins())){
    DEBUG_OUT_L(2) << "E= " << Etruth << " ran " << r << " -> ibin " << ibin << std::endl;
    return fEnergyBinning.Max() + 0.0001;
  }

  double Esmeared = fEnergyBinning.BinCenterLog(ibin);
  DEBUG_OUT_L(2) << "E= " << Etruth << " ran " << r << " -> ibin " << ibin << ", Esmeared= " << Esmeared << std::endl;

  return Esmeared;
}

void MigrationRandomGenerator::SetRandomGenerator(TRandom* rand) {

  if (fRandom)
    delete fRandom;

  fRandom = rand;
}
