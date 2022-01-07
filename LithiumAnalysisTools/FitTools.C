#include "FitTools.hh"

#include <TFitResult.h>

#define INFO_OUT_TAG "FitTools"
#include "debugging.hh"

int IsGoodFitResult(TFitResultPtr fitResult, bool expectMinosErrors) {

  if (fitResult->IsEmpty())
    return -1;

  if (fitResult->IsZombie())
    return -2;

  if (!fitResult->IsValid())
    return -3;

  if (std::isnan(fitResult->Edm()))
    return -4;

  if (std::isnan(fitResult->Chi2()))
    return -5;

  if (fitResult->Chi2() <= 0.0)
     return -6;

  if (fitResult->CovMatrixStatus() != 3) {
    if (fitResult->CovMatrixStatus() == 2)
      WARN_OUT << "Minuit: Covariance Matrix forced positive." << std::endl;
    return  -7;
  }

  for (unsigned int i = 0; i < fitResult->NPar(); ++i) {
    if (fitResult->IsParameterFixed(i))
      continue;

    if (expectMinosErrors && !fitResult->HasMinosError(i)) {
      if (fitResult->CovMatrixStatus() == 2)
        WARN_OUT << "Minuit: Parameter " << i << " " << fitResult->ParName(i) << " has no MINOS Error." << std::endl;
      return -8;
    }

    if (std::isnan(fitResult->ParError(i)))
      return -9;

    if (fitResult->ParError(i) <= 0)
      return -10;
  }

  return 0;
}
