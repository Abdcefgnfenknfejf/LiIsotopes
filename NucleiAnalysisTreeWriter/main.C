#include "NucleiAnalysisTree.hh"
#include "LithiumBinning.hh"
#include "SelectionCuts.hh"
// ACsoft includes
#include "AcceptanceManager.hh"
#include "AMSGeometry.h"
#include "AnalysisEvent.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "CustomSpectrum.hh"
#include "DetectorManager.hh"
#include "Environment.hh"
#include "EventFactory.hh"
#include "FileManager.hh"
#include "FileManagerController.hh"
#include "McSpectrumScaler.hh"
#include "ObjectManager.hh"
#include "PowerSpectrum.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "TrdTracking.hh"

// This is needed if running directly on AMS root files to use 'gbatch-like' code in this tree                                                                                                                                                                                                                          
#ifdef HAVE_AMS_SUPPORT
#include "DisableWarnings.h"
#include "root.h"
#include "Tofrec02_ihep.h"
#include "MagField.h"
#include "RichCharge.h"
#include "EnableWarnings.h"
#endif



#define INFO_OUT_TAG "NucleiAnalysisTreeWriter"
#include "debugging.hh"

#include <TFile.h>
#include <TH1D.h>

int main(int argc, char** argv) {

  // Command line option handling.
  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("NucleiAnalysisTreeWriter",
                            "Illustrates the usage of IO::TreeWriter to write ROOT trees from ACQt files.");

  config.AddHelpExample("Loop over given filelist.", "--inputlist list.txt");

  std::string inputList;
  config.GetValue("OPTIONS", "inputlist", inputList,
                  "List of ACQt input files (full path).");

  std::string resultDirectory;
  config.GetValue("OPTIONS", "resultdir", resultDirectory,
                  "General directory where result files should be stored. Current directory is used if option is not specified.");

  //  bool fillAmsVariables = config.GetFlag("OPTIONS", "fillAmsVariables",
  //				 "Fill additional variables from AMSROOT files.");
  std::string suffix;
  config.GetValue("OPTIONS", "suffix", suffix,
                 "A string identifier to be used in parallel computing, to uniquely identify result files.");


  bool fillAmsVariables = config.GetFlag("OPTIONS", "fillAmsVariables",
                                         "Fill additional variables from AMSROOT files.");

  // Load & parse cut configuration file
  std::string cutConfigfile = "/home/lk974211/lithiumanalysis/Configuration/MCGeneratorTreeWriter.cfg";  //Write a flag for MC or data.. to be done.

  Environment::ExpandEnvironmentVariables(cutConfigfile);
  config.Read(cutConfigfile);
  Cuts::SelectionParser selectionParser(config);

  // Setup file manager to process ACQt data.
  IO::FileManager fileManager(&config);

  Analysis::EventFactory* eventFactory = Analysis::EventFactory::Create(&config);
  Analysis::Event event;

  //Get hSpectrumLi from InputFile
   std::string hist_SpectrumScaler = "/home/lk974211/lithiumanalysis/Data/FluxMcSpectrumScaler.root";
   TFile* fileHistograms = TFile::Open(hist_SpectrumScaler.c_str(), "read");
   assert(fileHistograms);
   assert(fileHistograms->IsOpen());

   TH1D* hist_SpectrumLi = (TH1D*)fileHistograms->Get("hLithiumFluxTimesDeltaTimeVsRigidity");
   TH1D* hist_SpectrumBe = (TH1D*)fileHistograms->Get("hBerylliumFluxTimesDeltaTimeVsRigidity");
   TH1D* hist_SpectrumCarbon = (TH1D*)fileHistograms->Get("hCarbonFluxTimesDeltaTimeVsRigidity");

   TH1D* hist_SpectrumBoron = (TH1D*)fileHistograms->Get("hBoronFluxTimesDeltaTimeVsRigidity");
   TH1D* hist_SpectrumNitrogen = (TH1D*)fileHistograms->Get("hNitrogenFluxTimesDeltaTimeVsRigidity");
   TH1D* hist_SpectrumOxygen = (TH1D*)fileHistograms->Get("hOxigenFluxTimesDeltaTimeVsRigidity");

  
  // McSpectrumScaler for MC event weights
  Utilities::McSpectrumScaler scaler(&config, resultDirectory, suffix);
  scaler.SetTargetSpectrum(ParticleId::Li6, Utilities::CustomSpectrum { hist_SpectrumLi, Utilities::KinematicVariable::Rigidity, ParticleId::Li6 } );
  scaler.SetTargetSpectrum(ParticleId::Li7, Utilities::CustomSpectrum { hist_SpectrumLi, Utilities::KinematicVariable::Rigidity, ParticleId::Li7 } );

  scaler.SetTargetSpectrum(ParticleId::Be7, Utilities::CustomSpectrum { hist_SpectrumBe, Utilities::KinematicVariable::Rigidity, ParticleId::Be7 } );
  scaler.SetTargetSpectrum(ParticleId::Be9, Utilities::CustomSpectrum { hist_SpectrumBe, Utilities::KinematicVariable::Rigidity, ParticleId::Be9 } );
  scaler.SetTargetSpectrum(ParticleId::Be10, Utilities::CustomSpectrum { hist_SpectrumBe, Utilities::KinematicVariable::Rigidity, ParticleId::Be10 } );

  scaler.SetTargetSpectrum(ParticleId::B10, Utilities::CustomSpectrum { hist_SpectrumBoron, Utilities::KinematicVariable::Rigidity, ParticleId::B10} );
  scaler.SetTargetSpectrum(ParticleId::B11, Utilities::CustomSpectrum { hist_SpectrumBoron, Utilities::KinematicVariable::Rigidity, ParticleId::B11} );
  scaler.SetTargetSpectrum(ParticleId::C12, Utilities::CustomSpectrum { hist_SpectrumCarbon, Utilities::KinematicVariable::Rigidity, ParticleId::C12 } );
  scaler.SetTargetSpectrum(ParticleId::N14, Utilities::CustomSpectrum { hist_SpectrumNitrogen, Utilities::KinematicVariable::Rigidity, ParticleId::N14 } );
  scaler.SetTargetSpectrum(ParticleId::O16, Utilities::CustomSpectrum { hist_SpectrumOxygen, Utilities::KinematicVariable::Rigidity, ParticleId::O16} );



  eventFactory->RegisterMcSpectrumScaler(&scaler);
  
  // 'AuxiliaryObjectManager' holds all auxiliary histograms / selectors created while processing the ACQt files.
  // NOTE: You should NOT write a TTree together with other histograms etc. in the ROOT file. You most likely
  // want to merge your histograms / selectors from batch jobs, but not the trees. That's why it's a good idea
  // in general to split up in two files: one for holding the tree, one for the rest.
  Utilities::ObjectManager auxiliaryObjectManager("AuxiliaryObjectManager", &config, resultDirectory, suffix);
  auxiliaryObjectManager.SetPrefix("NucleiAnalysis_Auxiliary");

// Read geometry config file for acceptance manager.
  static Acceptance::AcceptanceManager* acceptanceManager = nullptr;
  std::string geometryConfigFile;

  //FixMe: how to choose a right configuration file for acceptance manager?
  if(!acceptanceManager) {
   acceptanceManager = new Acceptance::AcceptanceManager();
   geometryConfigFile = "/home/lk974211/lithiumanalysis/Configuration/BeFluxesGeometry.cfg";
   Environment::ExpandEnvironmentVariables(geometryConfigFile);
  }

  eventFactory->RegisterAcceptanceManager(acceptanceManager);

  //FixMe: do i have to generate histograms in the tree writer?
  //Store MC generated histograms for acceptance calculation
  const Binning::Definition& binning = LithiumRigidityBinning();
   const int maxSpecies = 12;
   TH1D* generatedHistogram[maxSpecies];
   ParticleId::Species particleSpecies[maxSpecies] = {ParticleId::Proton, ParticleId::Alpha, ParticleId::Li6, ParticleId::Li7, ParticleId::Be7, ParticleId::Be9, ParticleId::Be10,
                                                     ParticleId::B10, ParticleId::B11, ParticleId::C12, ParticleId::N14, ParticleId::O16};
   for (int species = 0; species < maxSpecies; species++) {
     generatedHistogram[species] = Make<TH1D>(Form("GeneratedEvents_%d", ParticleId::Id(particleSpecies[species])), "Generated events", binning);
     // Fill generated events sampled uniformly from the whole semisphere, not just the selected acceptance.                                                                                               
     scaler.FillHistogram(particleSpecies[species], generatedHistogram[species],Utilities::KinematicVariable::Momentum);
     auxiliaryObjectManager.Add(generatedHistogram[species], Form("McSpectrum_%d", ParticleId::Id(particleSpecies[species])));
   }

  if (!config.PerformChecksAfterOptionParsing())
    return EXIT_FAIL_CONFIG;

  if (!fileManager.ReadFileList(inputList))
    return EXIT_FAIL_FILEMANAGER;

  // Construct tree manager which will manage the output file to hold the resulting tree.
  // Construct tree.                                                                                                                                                                                                                                                                                                        
  NucleiAnalysisTree tree;
  tree.SetFillAmsVariables(fillAmsVariables);

  IO::TreeWriter treeWriter(&tree, IO::TreeOptions::DontWriteInMemoryBranches);
  std::string treeFileName = Utilities::ObjectManager::MakeStandardRootFileName(resultDirectory, "NucleiAnalysisTree", suffix);
  treeWriter.Initialize(treeFileName);

  // Load cut selector(s).
  Cuts::Selector* preselectionBadRunsRTI = auxiliaryObjectManager.Add(selectionParser.GetSelector("PreselectionBadRunsRTI"));
  Cuts::Selector* IonSelection = auxiliaryObjectManager.Add(selectionParser.GetSelector("IonSelection"));
  Cuts::Selector* McParticleInAmsAcceptanceSelection = auxiliaryObjectManager.Add(selectionParser.GetSelector("McParticleInAmsAcceptanceSelection"));
  
  // Begin event loop.
  INFO_OUT_ON_MASTER << "Looping over " << fileManager.GetEntries() << " events..." << std::endl;

#ifdef HAVE_AMS_SUPPORT
  if (fillAmsVariables) {
    const char* amsDataDir = getenv("AMSDataDir");
    if (!amsDataDir)
      FATAL_OUT << "AMSDataDir environment variable is not defined. Please set it to (typically): /cvmfs/ams.cern.ch/Offline/AMSDataDir" << std::endl;
    
    std::stringstream fieldMapName;
    fieldMapName << amsDataDir << "/v5.00/" << "MagneticFieldMapPermanent_NEW_FULL.bin";
    MagField::GetPtr()->Read(fieldMapName.str().c_str());
  }
#endif

  static int sProductionSteps = Analysis::CreateSplineTrack;
  bool firstEvent = true;
  while (fileManager.GetNextEvent()) {
    // Initialize the AcceptanceManager, after the DetectorManager received the run type (needed for MPI).                                                                                                                                                                                                                   

    if (firstEvent) {
      acceptanceManager->InitSetup(geometryConfigFile);
      firstEvent = false;
    }

    fileManager.DumpEventLoopProgress(20000);

    eventFactory->SetupEmptyEvent(event);
    eventFactory->CreateParticles(event);

    if (!preselectionBadRunsRTI->Passes(event))
      continue;

    if (!event.IsMC()) {
      if (!IonSelection->Passes(event))
        continue;
    }

    if(!McParticleInAmsAcceptanceSelection->Passes(event))
      continue;

    Analysis::Particle* particle = (Analysis::Particle*)event.PrimaryParticle();
    bool isHighZ = false;
    if (particle && particle->HasTrackerTrack())
      isHighZ = (particle->TrackerCharge() > 2.4);

    eventFactory->TrdTrackingObject()->SetMode(isHighZ ? Analysis::HighZMode : Analysis::NormalMode);

    // eventFactory->PerformTrdTracking(event);                                                                                                                                                                                                                                                                            
    // eventFactory->PerformTrdVertexFinding(event);                                                                                                                                                                                                                                                                         
    eventFactory->FillParticles(event, sProductionSteps);

#ifdef HAVE_AMS_SUPPORT
    if (fillAmsVariables) {
      EcalShowerR::enableAutomaticRecoveryOfBrokenBuilds = false;
      EcalShowerR::enableNewDeadCellTreatment = false;
      TofRecH::RebuildBetaHInReadHeader = false;
      RichPMTCalib::loadPmtCorrections = false;
      RichRingR::loadPmtCorrections = false;
      RichRingR::loadChargeUniformityCorrection = false;

      fileManager.AssociatedAMSEvent();
    }
#endif

    treeWriter.Fill(event);
  }

  // Print Preselection statistics
  preselectionBadRunsRTI->PrintSummary();
  IonSelection->PrintSummary();
  McParticleInAmsAcceptanceSelection->PrintSummary();
 
  // Finish writing tree file.
  treeWriter.Finish();

  // Write auxiliary output file.
  auxiliaryObjectManager.WriteToFile();

  return EXIT_SUCCESS;
}
