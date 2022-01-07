#include "MagneticField.hh"

#define INFO_OUT_TAG "MagneticField"
#include "debugging.hh"

#ifdef HAVE_AMS_SUPPORT

#ifndef _PGTRACK_
#define _PGTRACK_
#endif

#include "DisableWarnings.h"
#include "MagField.h"
#include "EnableWarnings.h"

void InitializeMagneticField() {

  const char* amsDataDir = getenv("AMSDataDir");
  if (!amsDataDir)
    FATAL_OUT << "AMSDataDir environment variable is not defined. Please set it to (typically): /cvmfs/ams.cern.ch/Offline/AMSDataDir" << std::endl;

  std::stringstream magneticFieldMapName;
  magneticFieldMapName << amsDataDir << "/v5.00/" << "MagneticFieldMapPermanent_NEW_FULL.bin";
  MagField::GetPtr()->Read(magneticFieldMapName.str().c_str());
}

void MagneticField(float* location, float *B) {

  MagField::GetPtr()->GuFld(location, B);

  // GuFld returns kilogauss.
  B[0] *= 0.1;
  B[1] *= 0.1;
  B[2] *= 0.1;
}
#else
void InitializeMagneticField() {

  FATAL_OUT << "You have to recompile with AMS ROOT support to be able to use this function. Aborting!" << std::endl;
}

void MagneticField(float*, float*) {

  FATAL_OUT << "You have to recompile with AMS ROOT support to be able to use this function. Aborting!" << std::endl;
}
#endif
