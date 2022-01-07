set(source_files
  BayesUnfoldingWithCutoff.C
  CompositeParameterizationFunction.C
  CrystalBallParameterization.C
  CompositeParameterization.C
  EnergyScaleDeviationModel.C
  EnergyScaleForwardFoldingDataset.C
  FoldedAcceptanceUnfoldingForLeptons.C
  MatrixTools.C
  MigrationRandomGenerator.C
  ModelFits.C
  UnfoldingInput.C
  UnfoldingSystematicErrors.C
  UnfoldingToyMc.C
  UnfoldingTools.C
)

set(header_files
  BayesUnfoldingWithCutoff.hh
  CompositeParameterizationFunction.hh
  CrystalBallParameterization.hh
  CompositeParameterization.hh
  EnergyScaleDeviationModel.hh
  EnergyScaleForwardFoldingDataset.hh
  FoldedAcceptanceUnfoldingForLeptons.hh
  MatrixTools.hh
  MigrationRandomGenerator.hh
  MigrationParameterization.hh
  ModelFits.hh
  UnfoldingInput.hh
  UnfoldingSystematicErrors.hh
  UnfoldingToyMc.hh
  UnfoldingTools.hh
)

add_project_library(LeptonUnfoldingTools SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(LeptonUnfoldingTools ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} RooUnfold LeptonAnalysisTools)
