set(source_files RooGaussDoubleSidedExp.C
		 RooGaussExp.C
		 TemplateParametrisation.C)
set(header_files RooGaussDoubleSidedExp.hh
		 RooGaussExp.hh
		 TemplateParametrisation.hh)

set(header_files_with_classdef RooGaussExp.hh RooGaussDoubleSidedExp.hh TemplateParametrisation.hh)

add_project_library(AnalyticTemplateParametrisation SOURCES ${source_files} HEADERS ${header_files} HEADERS_WITH_CLASSDEF ${header_files_with_classdef} LINKDEF Linkdef.h)
target_link_libraries(AnalyticTemplateParametrisation ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} TreePlayer RooFitCore)
