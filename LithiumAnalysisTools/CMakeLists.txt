set(source_files
  AnalysisSettings.C
  DynamicSelector.C
  EcalBDTEfficiencyParameterization.C
  EcalChiSquareLateralCalibration.C
  EfficiencyTools.C
  EffectiveMeasuringTime.C
  ElectronModelTwoPowerLawsWithSolarModulation.C
  ElectronAndPositronFluxModel.C
  ElectronPositronModelWithSmoothDoubleBreak.C
  ElectronPositronDoublePowerLawApproximation.C
  EnergyScaleTree.C
  FitTools.C
  FluxModel.C
  GraphTools.C
  MagneticField.C
  MCGeneratorTree.C
  TrdEstimatorTools.C
  LeptonAnalysisTree.C
  LeptonTreePreselection.C
  PositronRatioTools.C
  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff.C
  SolarFlares.C
  SimplePowerLawModel.C
  SpectralIndex.C
  SystematicUncertainties.C
  MeasuringTimeParameterization.C
  ModulatedPowerLawModel.C
  TemplateFitTree.C
  TemplateFitTools.C
  TemplateFit2DFitFunction.C
  TrdLikelihoodRatioCalibration.C
  MultiTrackChargeConfusionStudyTree.C
  RooFitTools.C
  NucleiAnalysisTree.C
  EfficiencyTemplateFits.C
)

set(header_files
  AveragingTools.hh
  AnalysisSettings.hh
  CanvasTools.hh
  DynamicSelector.hh
  EcalBDTEfficiencyParameterization.hh
  EcalChiSquareLateralCalibration.hh
  EfficiencyTools.hh
  ChargeConfusionTools.hh
  EffectiveMeasuringTime.hh
  ElectronModelTwoPowerLawsWithSolarModulation.hh
  ElectronAndPositronFluxModel.hh
  ElectronPositronModelWithSmoothDoubleBreak.hh
  ElectronPositronDoublePowerLawApproximation.hh
  EnergyScaleTree.hh
  FitTools.hh
  FluxModel.hh
  MagneticField.hh
  MCGeneratorTree.hh
  TrdEstimatorTools.hh
  LeptonAnalysisTree.hh
  LeptonTreePreselection.hh
  MultiTrackChargeConfusionStudyTree.hh
  SolarFlares.hh
  SimplePowerLawModel.hh
  SpectralIndex.hh
  SystematicUncertainties.hh
  MeasuringTimeParameterization.hh
  ModulatedPowerLawModel.hh
  PositronRatioTools.hh
  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff.hh
  RooFitTools.hh
  SplineFitTools.hh
  TemplateFitTree.hh
  TemplateFitTools.hh
  TemplateFit2DFitFunction.hh
  TimeTools.hh
  TrdLikelihoodRatioCalibration.hh
  TriggerEfficiency.hh
  NucleiAnalysisTree.hh
)

set(header_files_with_classdef
  EfficiencyTemplateFits.hh
)

# Work-around for gbatch support, see ACANA-24.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
endif()

# These additional ROOT libraries are required.
set(ADDITIONAL_ROOT_LIBS TreePlayer RooFitCore)

add_project_library(LithiumAnalysisTools SOURCES ${source_files} HEADERS ${header_files} HEADERS_WITH_CLASSDEF ${header_files_with_classdef} LINKDEF Linkdef.h)
target_link_libraries(LithiumAnalysisTools ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} ${ADDITIONAL_ROOT_LIBS} LithiumBinning)
