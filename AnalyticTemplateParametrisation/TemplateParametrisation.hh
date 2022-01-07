#ifndef TemplateParametrisation_hh
#define TemplateParametrisation_hh

#include "Utilities.hh"
#include <string>
#include "ParticleId.hh"
#include "LithiumBinning.hh"

class TH1;
class TH1D;
class TH2;
class TGraphAsymmErrors;
class TMultiGraph;
class TAttLine;


double MyGaus(double* R, double* par);
double function_MassTemplateWidth(double* xEkinN, double* par);
double AsyGaus(double* x, double* par);
double GausPlusAsyGausFunction(double *x, double* par);
double DoubleGausFunction(double *x, double* par);
double TwoIsotopesMassFunction(double *x, double* par);
double TwoIsotopesMassFunction_FixFactor(double *x, double* par);

void SetGraphStyle(TGraph* graphs, Style_t markerStyle, Color_t markerColor, Size_t markerSize);
void SetFunctionStyle(TF1* Fuction, Style_t lineStyle, Color_t lineColor, Size_t lineWidth);
void SetHistStyle(TH1* Hist, Style_t MarkerStyle, Color_t MarkerColor, Color_t lineColor, Size_t MarkerSize);
void SetGraphAxis(TMultiGraph* multigraph_SigmaPrimary);
void SetResidualPlot(TH1* hist, TH1* pullHistogram, TF1* Function);
void format_line(TAttLine* line,int col,int sty);

#endif
