set(source_files LithiumTemplateFitFunction.C LithiumTemplateFitter.C)
set(header_files LithiumTemplateFitFunction.hh LithiumTemplateFitter.hh)

add_project_library(LiTemplateFitter SOURCES ${source_files} HEADERS ${header_files} HEADERS_WITH_CLASSDEF ${header_files_lib_with_classdef} LINKDEF Linkdef.h)
target_link_libraries(LiTemplateFitter ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} LithiumAnalysisTools TreePlayer)
