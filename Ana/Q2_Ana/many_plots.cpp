#include "many_plots.h"

many_plots::many_plots()
{
}

many_plots::many_plots(std::string temp_name, std::string temp_title, double xmin, double xmax)
{
  h_ep = new TH1D((temp_name+"_ep").c_str(),(temp_title+";"+temp_title+";Counts").c_str(),100,xmin,xmax);
  h_epp = new TH1D((temp_name+"_epp").c_str(),(temp_title+";"+temp_title+";Counts").c_str(),100,xmin,xmax);

  for(int i=0; i<10; i++){
    h_ep_Q2bin[i] = new TH1D((temp_name+"_ep_Q2bin"+std::to_string(i)).c_str(),(temp_title+";"+temp_title+";Counts").c_str(),50,xmin,xmax);
    h_epp_Q2bin[i] = new TH1D((temp_name+"_epp_Q2bin"+std::to_string(i)).c_str(),(temp_title+";"+temp_title+";Counts").c_str(),50,xmin,xmax);
  }

}

many_plots::~many_plots()
{
}

void many_plots::Fill_hist_set(bool is_epp, double Q2, double x)
{
  h_ep->Fill(x);
  h_ep_Q2bin[binQ2(Q2)]->Fill(x);
  if(is_epp){
    h_epp->Fill(x);
    h_epp_Q2bin[binQ2(Q2)]->Fill(x);
  }  
}

void many_plots::Write_hist_set(TFile *f, char fileName[100], TCanvas * myCanvas)
{
  f->cd();
  h_ep->Write();
  h_epp->Write();
  for(int i=0; i<10; i++){
    h_ep_Q2bin[i]->Write();
    h_epp_Q2bin[i]->Write();
  }

  myCanvas->Divide(4,3);
  myCanvas->cd(1);
  h_ep->Scale(1/h_ep->Integral());
  h_ep->SetLineColor(1);
  h_ep->Draw();
  h_epp->Scale(1/h_epp->Integral());
  h_epp->SetLineColor(2);
  h_epp->Draw("SAME");
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_ep_Q2bin[i]->Scale(1/h_ep_Q2bin[i]->Integral());
    h_ep_Q2bin[i]->SetLineColor(1);
    h_ep_Q2bin[i]->Draw();
    h_epp_Q2bin[i]->Scale(1/h_epp_Q2bin[i]->Integral());
    h_epp_Q2bin[i]->SetLineColor(2);
    h_epp_Q2bin[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


}

int many_plots::binQ2(double q2){
  if(q2<1.65){return 0;}
  else if(q2<1.80){return 1;}
  else if(q2<1.95){return 2;}
  else if(q2<2.10){return 3;}
  else if(q2<2.25){return 4;}
  else if(q2<2.40){return 5;}
  else if(q2<2.70){return 6;}
  else if(q2<3.00){return 7;}
  else if(q2<3.50){return 8;}
  else{return 9;}
}
