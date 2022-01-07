// ACsoft includes
#include "BinningDefinition.hh"
#include "ConfigHandler.hh"
#include "CutFactory.hh"
#include "Environment.hh"
#include "FileManagerController.hh"
#include "GlobalOptions.hh"
#include "MPIEnvironment.hh"
#include "MeasuringTime.hh"
#include "ObjectManager.hh"
#include "LithiumBinning.hh"

#define INFO_OUT_TAG "rti_reader_example"
#include "debugging.hh"

int main(int argc, char** argv) {

  // Command line option handling.
  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("rti_reader_example",
                            "Illustrates the calculation of total livetime for a given time period based on MeasuringTime class.");

  config.AddHelpExample("Case #1): Loop directly over RTI files.", "--startTime 1305853512 --endTime 1546404002");
  config.AddHelpExample("Case #2): Loop over RTI tree written by rti_writer_example.", "--RTI/TreePattern RTITree*.root --startTime 1305853512 --endTime 1527490046");

  int startTime = -1;
  config.GetValue("OPTIONS", "startTime", startTime,
                  "Start time in seconds as Unix time stamp.");

  int endTime = -1;
  config.GetValue("OPTIONS", "endTime", endTime,
                  "End time in seconds as Unix time stamp.");

  Utilities::ObjectManager objectManager(&config, "", "");
  if (!config.PerformChecksAfterOptionParsing())
    return EXIT_FAIL_CONFIG;

  if (startTime == -1 || endTime == -1)
    FATAL_OUT << "You have to specify both --startTime X / --endTime Y." << std::endl;

  // Preperations to calculate the measuring time:
  // 1. Select a geomagnetic cut-off cut, the field-of-view (here: 40 degress for both positive/negative particle hypothesis) and a safety factor
  Cuts::Cut* cutOffCut = Cuts::CreateCut("RigidityAboveGeomagneticCutoff", "40PN", 1.2);

  // 2. Open the cut config file containing at least "BadRuns" / "RTI" selectors, optionally also "TrdCalibration".
  std::string cutConfigfile = "/home/lk974211/lithiumanalysis/Configuration/NucleiAnalysisTreeWriter.cfg";
  Environment::ExpandEnvironmentVariables(cutConfigfile);
  config.Read(cutConfigfile);

#ifdef ENABLE_MPI
  if (IO::MPIEnvironment::IsMPIEnabled())
    MPI_Init(&argc, &argv);
#endif

  // NOTE: These two lines are only needed in a standalone program that does NOT loop over ACQt files.
  // If you've instantiated an IO::FileManager before, this is automatically deduced from the file list.
  IO::FileManagerController::Self()->SetRunType(AC::ISSRun);
  IO::FileManagerController::Self()->SetFirstAndLastEventTimes(startTime, endTime);

  const std::vector<std::string>& inputTreeFiles = GlobalOptions::Self()->RTITreeFileName();
  if (inputTreeFiles.empty())
    objectManager.SetPrefix("MeasuringTime_from_RTI_Files");
  else
    objectManager.SetPrefix("MeasuringTime_from_RTI_Tree");

  const Binning::Definition RigidityBinning = LithiumRigidityBinning();
  RTI::MeasuringTime measuringTime(config, objectManager, cutOffCut, RigidityBinning);
  measuringTime.ComputeMeasuringTime();

  objectManager.WriteToFile();

#ifdef ENABLE_MPI
  if (IO::MPIEnvironment::IsMPIEnabled())
    MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
