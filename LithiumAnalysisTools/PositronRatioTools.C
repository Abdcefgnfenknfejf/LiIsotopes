#include "PositronRatioTools.hh"

#include <cmath>

void CalculatePositronElectronRatio(double& positronElectronRatio, double& positronElectronRatioUncertainty,
                                    double electrons, double electronsUncertainty,
                                    double positrons, double positronsUncertainty) {

  if (electrons > 0 && positrons > 0) {
    positronElectronRatio = positrons / electrons;
    positronElectronRatioUncertainty = positronElectronRatio * std::sqrt(std::pow(positronsUncertainty / positrons, 2) +
                                                                         std::pow(electronsUncertainty / electrons, 2));
  } else {
    positronElectronRatio = 0.0;
    positronElectronRatioUncertainty = 0.0;
  }
}

Utilities::Quantity CalculatePositronElectronRatio(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons) {

  double positronElectronRatio = 0.0;
  double positronElectronRatioUncertainty = 0.0;
  CalculatePositronElectronRatio(positronElectronRatio, positronElectronRatioUncertainty,
                                 electrons.value, electrons.uncertainty, positrons.value, positrons.uncertainty);

  return { positronElectronRatio, positronElectronRatioUncertainty };
}

void CalculatePositronElectronRatio(double& positronElectronRatio, double& positronElectronRatioUncertainty,
                                    double electrons, double electronsUncertainty,
                                    double positrons, double positronsUncertainty,
                                    double electronPositronCorrelation) {

  if (electrons > 0 && positrons > 0) {
    positronElectronRatio = positrons / electrons;

    double regularTerm = std::pow(positronsUncertainty / positrons, 2) + std::pow(electronsUncertainty / electrons, 2);
    double correlationTerm = 2.0 * electronPositronCorrelation * electronsUncertainty * positronsUncertainty / electrons / positrons;

    if (regularTerm >= correlationTerm)
      positronElectronRatioUncertainty = positronElectronRatio * std::sqrt(regularTerm - correlationTerm);
    else {
      positronElectronRatio = 0.0;
      positronElectronRatioUncertainty = 0.0;
    }
  } else {
    positronElectronRatio = 0.0;
    positronElectronRatioUncertainty = 0.0;
  }
}

Utilities::Quantity CalculatePositronElectronRatio(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons, double electronPositronCorrelation) {

  double positronElectronRatio = 0.0;
  double positronElectronRatioUncertainty = 0.0;
  CalculatePositronElectronRatio(positronElectronRatio, positronElectronRatioUncertainty,
                                 electrons.value, electrons.uncertainty, positrons.value, positrons.uncertainty, electronPositronCorrelation);

  return { positronElectronRatio, positronElectronRatioUncertainty };
}

void CalculatePositronFraction(double& positronFraction, double& positronFractionUncertainty,
                               double electrons, double electronsUncertainty,
                               double positrons, double positronsUncertainty) {

  if (electrons > 0 && positrons > 0) {
    double positronElectronRatio = 0.0;
    double positronElectronRatioUncertainty = 0.0;
    CalculatePositronElectronRatio(positronElectronRatio, positronElectronRatioUncertainty, electrons, electronsUncertainty, positrons, positronsUncertainty);

    positronFraction = 1.0 - 1.0 / (positronElectronRatio + 1.0);
    positronFractionUncertainty = positronElectronRatioUncertainty / std::pow(positronElectronRatio + 1.0, 2);
  } else {
    positronFraction = 0.0;
    positronFractionUncertainty = 0.0;
  }
}

Utilities::Quantity CalculatePositronFraction(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons) {

 double positronFraction = 0.0;
 double positronFractionUncertainty = 0.0;
 CalculatePositronFraction(positronFraction, positronFractionUncertainty,
                           electrons.value, electrons.uncertainty, positrons.value, positrons.uncertainty);

 return { positronFraction, positronFractionUncertainty };
}

void CalculatePositronFraction(double& positronFraction, double& positronFractionUncertainty,
                               double electrons, double electronsUncertainty,
                               double positrons, double positronsUncertainty,
                               double electronPositronCorrelation) {

  if (electrons > 0 && positrons > 0) {
    double positronElectronRatio = 0.0;
    double positronElectronRatioUncertainty = 0.0;
    CalculatePositronElectronRatio(positronElectronRatio, positronElectronRatioUncertainty, electrons, electronsUncertainty, positrons, positronsUncertainty, electronPositronCorrelation);

    positronFraction = 1.0 - 1.0 / (positronElectronRatio + 1.0);
    positronFractionUncertainty = positronElectronRatioUncertainty / std::pow(positronElectronRatio + 1.0, 2);
  } else {
    positronFraction = 0.0;
    positronFractionUncertainty = 0.0;
  }
}

Utilities::Quantity CalculatePositronFraction(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons, double electronPositronCorrelation) {

 double positronFraction = 0.0;
 double positronFractionUncertainty = 0.0;
 CalculatePositronFraction(positronFraction, positronFractionUncertainty,
                           electrons.value, electrons.uncertainty, positrons.value, positrons.uncertainty, electronPositronCorrelation);

 return { positronFraction, positronFractionUncertainty };
}
