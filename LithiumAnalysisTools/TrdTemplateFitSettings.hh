#ifndef TrdTemplateFitSettings_hh
#define TrdTemplateFitSettings_hh

// Time-averaged settings
static const unsigned int sComplexModelStopBinElectronSingleTrack = 57;
static const unsigned int sComplexModelStartBinProtonSingleTrack = 3;
static const unsigned int sComplexModelStopBinProtonSingleTrack = 67;
static const unsigned int sKeepProtWidthGausDeltaConstantAboveBinSingleTrack = 24;

static const unsigned int sComplexModelStopBinElectronMultiTracks = 42;
static const unsigned int sComplexModelStartBinProtonMultiTracks = 3;
static const unsigned int sComplexModelStopBinProtonMultiTracks = 67;

static const unsigned int sUseBinnedLikelihoodFitsAboveBin = 0;
static const unsigned int sUseBinnedLikelihoodFitsBelowBin = 64;

// Time-dependant settings
static const unsigned int sTimeDependantFitStopCCProtonTemplateFitAboveBin = 33;
static const unsigned int sBartelsRotationAfterWhichOnlySmallTrdOperationChangesAreExpected = 15;

#endif
