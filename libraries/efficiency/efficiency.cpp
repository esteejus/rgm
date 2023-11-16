#include "efficiency.h"



// efficiency constants
double Mlow = 0.85;
double Mhigh = 1.05;
double Mdisp_lo = 0.5;
double Mdisp_hi = 1.5;






Double_t signal(Double_t *x, Double_t *par) { // height, mean, width
  return par[0]*exp(-pow((x[0]-par[1]),2.)/(2*pow(par[2],2.))); 
}
Double_t mmiss_signal_gauss(Double_t *x, Double_t *par) {
  return signal(x,par) + signal(x,&par[3]);
}




// input: 2d histogram of missing mass vs either momentum or theta
// output: create 1d histograms in intervals of momentum or theta, fit to 2 Gaussians (optional subtract background)
// output: interval of missing mass peak from Mlow to Mhigh
// if subtract_bk = False, the missing mass histograms will be fit to the sum of 2 Gaussians (signal + background)
// if subtract_bk = True, the background fit will be subtracted from the missing mass histogram, and the signal fit displayed
double * hist_projections_backsub(TCanvas * can, TH2D * hist2d, int num_hist, bool subtract_bk, char v)
{
  double p_start_val[num_hist];
  double x_min = hist2d->GetXaxis()->GetXmin();
  double x_max = hist2d->GetXaxis()->GetXmax();
  double dp = (x_max-x_min)/num_hist;
  double * S = new double[num_hist];
  // plot and fit each graph
  for (int i=0; i<num_hist; i++)
  {
    // set momentum/theta interval for this missing mass projection
    p_start_val[i] = x_min + i*dp;
    int bin1 = hist2d->GetXaxis()->FindBin(p_start_val[i]);
    int bin2 = hist2d->GetXaxis()->FindBin(p_start_val[i]+dp);
    // make projection for x interval
    can->cd(i+1);
    TH1D * proj = hist2d->ProjectionY("",bin1,bin2,"d");

    // create name of missing mass histogram for current momentum/theta interval
    std::ostringstream sObj1, sObj2;
    std::string leftTitle = "Missing Mass in ("; std::string midTitle = ",";
    std::string rightTitle;
    if (v=='p')
    {
      rightTitle = ") GeV/c";
      sObj1 << std::fixed << std::setprecision(3) << p_start_val[i];
      sObj2 << std::fixed << std::setprecision(3) << p_start_val[i] + dp;
    }
    else if (v=='t')
    {
      rightTitle = ") deg";
      sObj1 << std::fixed << std::setprecision(0) << p_start_val[i];
      sObj2 << std::fixed << std::setprecision(0) << p_start_val[i] + dp;
    }
    else
    {
      std::cout << "Invalid projection variable for missing mass\n";
    }
    std::string result = leftTitle + sObj1.str() + midTitle + sObj2.str() + rightTitle;
    proj->SetTitle(result.c_str());

    // fit histogram to Gaussian (signal) + Gaussian (background)
    TF1 * cfit = new TF1("cfit",mmiss_signal_gauss,Mdisp_lo,Mdisp_hi,6);
    cfit->SetLineColor(kMagenta);
    cfit->SetParameters(10000,0.94,0.05,1200,1.3,0.1);
    cfit->SetParLimits(0,0,1000000);
    cfit->SetParLimits(1,0.92,1.01);
    cfit->SetParLimits(2,0.02,0.1);
    cfit->SetParLimits(3,0,1000000);
    cfit->SetParLimits(4,1.1,1.8);
    cfit->SetParLimits(5,0.05,0.3);
    proj->Fit("cfit","QN");
    if (subtract_bk)  {cfit->Draw("same");}
    // separate signal/background
    // fit background
    Double_t par[6];
    cfit->GetParameters(par);
    TF1 * bkfit = new TF1("backfit",signal,Mdisp_lo,Mdisp_hi,3);
    //bkfit->SetLineColor(kBlue);
    bkfit->SetParameters(&par[3]);
    if (subtract_bk)  {bkfit->Draw("same");}
    // fit signal
    TF1 * sgfit = new TF1("sgfit",signal,Mdisp_lo,Mdisp_hi,3);
    sgfit->SetLineColor(kGreen);
    sgfit->SetParameters(par);
    // background subtraction!
    if (subtract_bk)
    {
      proj->Draw("same");
      cfit->SetLineColor(kGreen); cfit->Draw("same");
      sgfit->SetLineColor(kRed);  sgfit->Draw("same");
      bkfit->SetLineColor(kBlue); bkfit->Draw("same");
    }
    if (!subtract_bk)
    {
      cfit->SetLineColor(kGreen);  cfit->Draw("same");
      bkfit->SetLineColor(kBlue);  bkfit->Draw("same");
    }
    // find background-subtracted signal

    TH1D * proj_sub = (TH1D*)proj->Clone();
    if (subtract_bk) {proj_sub->Add(bkfit,-1);}
    S[i] = proj_sub->Integral(proj_sub->GetXaxis()->FindBin(Mlow),proj_sub->GetXaxis()->FindBin(Mhigh));
    //S[i] = sgfit->Integral(Mlow,Mhigh);
        /*std::ostringstream sObj3;
        sObj3 << std::fixed << std::setprecision(3) << S[i];
        std::string result = sObj3.str();// leftTitle + sObj1.str() + midTitle + sObj2.str() + rightTitle + ;
        proj->SetTitle(result.c_str());*/
    //std::cout << "integral = " << proj->Integral(Mlow,Mhigh) << '\n';
    //std::cout << "integral and error = " << proj->IntegralAndError(Mlow,Mhigh) << '\n';
  }
  return S;
}

