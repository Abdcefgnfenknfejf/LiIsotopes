set(source_files_lib NucleiAnalysisSelection.C)
set(header_files_lib_with_classdef NucleiAnalysisSelection.hh)
add_project_library(NucleiAnalysisSelectionLibrary SOURCES ${source_files_lib} HEADERS_WITH_CLASSDEF ${header_files_lib_with_classdef} LINKDEF Linkdef.h LIBS NucleiAnalysisTree LithiumBinning SelectionCuts)

set(source_files_exe main.C)
add_project_executable(NucleiAnalysisSelection SOURCES ${source_files_exe} LIBS NucleiAnalysisSelectionLibrary NucleiAnalysisTree LithiumBinning  SelectionCuts) 
