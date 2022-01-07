#ifndef PositronRatioTools_hh
#define PositronRatioTools_hh

#include "Quantity.hh"

void CalculatePositronElectronRatio(double& positronElectronRatio, double& positronElectronRatioUncertainty,
                                    double electrons, double electronsUncertainty,
                                    double positrons, double positronsUncertainty);

Utilities::Quantity CalculatePositronElectronRatio(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons);

void CalculatePositronElectronRatio(double& positronElectronRatio, double& positronElectronRatioUncertainty,
                                    double electrons, double electronsUncertainty,
                                    double positrons, double positronsUncertainty,
                                    double electronPositronCorrelation);

Utilities::Quantity CalculatePositronElectronRatio(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons, double electronPositronCorrelation);

void CalculatePositronFraction(double& positronFraction, double& positronFractionUncertainty,
                               double electrons, double electronsUncertainty,
                               double positrons, double positronsUncertainty);

Utilities::Quantity CalculatePositronFraction(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons);

void CalculatePositronFraction(double& positronFraction, double& positronFractionUncertainty,
                               double electrons, double electronsUncertainty,
                               double positrons, double positronsUncertainty,
                               double electronPositronCorrelation);

Utilities::Quantity CalculatePositronFraction(const Utilities::Quantity& electrons, const Utilities::Quantity& positrons, double electronPositronCorrelation);

#endif
