#ifndef Linkdef_h
#define Linkdef_h

#ifdef __CINT__

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedefs;

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RooGaussExp+;
#pragma link C++ class RooGaussDoubleSidedExp+;

#endif

#endif // Linkdef_h
