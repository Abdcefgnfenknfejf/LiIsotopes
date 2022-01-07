set(source_files main.C)
add_project_executable(NucleiAnalysisTreeWriter SOURCES ${source_files} HEADERS ${header_files} LIBS NucleiAnalysisTree LithiumBinning SelectionCuts)
