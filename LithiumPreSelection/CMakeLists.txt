set(source_files_lib LithiumPreSelection.C)
set(header_files_lib_with_classdef LithiumPreSelection.hh)
add_project_library(LithiumPreSelectionLibrary SOURCES ${source_files_lib} HEADERS_WITH_CLASSDEF ${header_files_lib_with_classdef} LINKDEF Linkdef.h LIBS NucleiAnalysisTree SelectionCuts)

set(source_files_exe main.C)
add_project_executable(LithiumPreSelection SOURCES ${source_files_exe} LIBS NucleiAnalysisTree LithiumBinning LithiumPreSelectionLibrary)
