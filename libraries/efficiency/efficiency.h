#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include <sstream>
#include <iomanip>

#include "clas12reader.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

using namespace clas12;

Double_t signal(Double_t *x, Double_t *par);

Double_t mmiss_signal_gauss(Double_t *x, Double_t *par);

double * hist_projections_backsub(TCanvas * can, TH2D * hist2d, int num_hist, bool subtract_bk, char v);

#endif
