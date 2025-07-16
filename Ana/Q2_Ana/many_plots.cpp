#include "many_plots.h"

many_plots::many_plots()
{
}

many_plots::many_plots(std::string temp_name, std::string temp_title, double xmin, double xmax)
{
  h_ep = new TH1D((temp_name+"_ep").c_str(),(temp_title+" (e,e'p);"+temp_title+";Counts").c_str(),100,xmin,xmax);
  h_epp = new TH1D((temp_name+"_epp").c_str(),(temp_title+" (e,e'pp);"+temp_title+";Counts").c_str(),100,xmin,xmax);

  char temp_Q2[100];
  double b[11]= {1.50,1.65,1.80,1.95,2.10,2.25,2.40,2.70,3.00,3.50,5.00};
  
  for(int i=0; i<10; i++){
    sprintf(temp_Q2,"%f - %f",b[i],b[i+1]);
    h_ep_Q2bin[i] = new TH1D((temp_name+"_ep_Q2bin"+std::to_string(i)).c_str(),temp_Q2,50,xmin,xmax);
    h_epp_Q2bin[i] = new TH1D((temp_name+"_epp_Q2bin"+std::to_string(i)).c_str(),temp_Q2,25,xmin,xmax);
  }

}

many_plots::many_plots(std::string temp_name, std::string temp_title, TFile * inFile)
{
  h_ep = (TH1D*)inFile->Get((temp_name+"_ep").c_str());
  h_epp = (TH1D*)inFile->Get((temp_name+"_epp").c_str());

  char temp_Q2[100];
  double b[11]= {1.50,1.65,1.80,1.95,2.10,2.25,2.40,2.70,3.00,3.50,5.00};
  
  for(int i=0; i<10; i++){
    sprintf(temp_Q2,"%f - %f",b[i],b[i+1]);
    h_ep_Q2bin[i] = (TH1D*)inFile->Get((temp_name+"_ep_Q2bin"+std::to_string(i)).c_str());
    h_epp_Q2bin[i] = (TH1D*)inFile->Get((temp_name+"_epp_Q2bin"+std::to_string(i)).c_str());
  }

}

many_plots::~many_plots()
{
}

void many_plots::Fill_hist_set(bool is_epp, double Q2, double x, double wep, double wepp)
{
  h_ep->Fill(x,wep);
  h_ep_Q2bin[binQ2(Q2)]->Fill(x,wep);
  if(is_epp){
    h_epp->Fill(x,wepp);
    h_epp_Q2bin[binQ2(Q2)]->Fill(x,wepp);
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
  /*
  myCanvas->Divide(4,3);
  myCanvas->cd(1);
  h_ep->Scale(h_epp->Integral()/h_ep->Integral());
  h_ep->SetLineColor(2);
  h_ep->Draw();
  h_epp->SetLineColor(1);
  h_epp->Draw("SAME");
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_ep_Q2bin[i]->Scale(1/h_ep_Q2bin[i]->Integral());
    h_ep_Q2bin[i]->SetLineColor(2);
    h_ep_Q2bin[i]->Draw();
    h_epp_Q2bin[i]->Scale(1/(2*h_epp_Q2bin[i]->Integral()));
    h_epp_Q2bin[i]->SetLineColor(1);
    h_epp_Q2bin[i]->Draw("SAME");
  }
  for(int i=0; i<10; i++){
    myCanvas->cd(2);
    h_ep_Q2bin[i]->SetLineColor(i+1);
    h_ep_Q2bin[i]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */

}

void many_plots::Write_hist_set_epp(TFile *f, char fileName[100], TCanvas * myCanvas)
{
  f->cd();
  h_epp->Write();
  for(int i=0; i<10; i++){
    h_epp_Q2bin[i]->Write();
  }
  /*
  myCanvas->Divide(4,3);
  myCanvas->cd(1);
  h_epp->SetLineColor(1);
  h_epp->Draw("SAME");
    for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_epp_Q2bin[i]->SetLineColor(1);
    h_epp_Q2bin[i]->Draw("SAME");
  }
  myCanvas->cd(2);
  for(int i=0; i<10; i++){
    TH1D * h_epp_clone = (TH1D*) h_epp_Q2bin[i]->Clone();
    h_epp_clone->SetLineColor(i);
    h_epp_clone->Scale(1/h_epp_Q2bin[i]->Integral());
    h_epp_clone->Draw("SAME");
  }

  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */

}

void many_plots::Write_ratio_set(TFile *f, char fileName[100], TCanvas * myCanvas)
{
  f->cd();
  TH1D * h_epp_clone = (TH1D*) h_epp->Clone();
  h_epp_clone->Divide(h_ep);
  h_epp_clone->Write();

  TH1D * h_epp_clone_Q2bin[10];
  for(int i=0; i<10; i++){
    h_epp_clone_Q2bin[i] = (TH1D*) h_epp_Q2bin[i]->Clone();
    h_epp_clone_Q2bin[i]->Divide(h_ep_Q2bin[i]);
    h_epp_clone_Q2bin[i]->Write();
  }

  myCanvas->Divide(4,3);
  myCanvas->cd(1);
  h_epp_clone->Draw();  
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_epp_clone_Q2bin[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
}

void many_plots::Draw_hist_set(char fileName[100], TCanvas * myCanvas)
{
  myCanvas->cd(1);
  h_ep->SetLineColor(1);
  h_ep->Draw();
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_ep_Q2bin[i]->SetLineColor(1);
    h_ep_Q2bin[i]->Draw();
  }

}

void many_plots::Draw_hist_set_same(char fileName[100], TCanvas * myCanvas, many_plots scaler)
{
  double s = scaler.getScale_ep()/h_ep->Integral();
  myCanvas->cd(1);
  h_ep->Scale(s);
  h_ep->SetLineColor(2);
  h_ep->Draw("SAME");
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    double si = scaler.getScale_epQ2(i)/h_ep_Q2bin[i]->Integral();
    h_ep_Q2bin[i]->Scale(si);
    h_ep_Q2bin[i]->SetLineColor(2);
    h_ep_Q2bin[i]->Draw("SAME");
  }

}

void many_plots::Draw_hist_set_epp(char fileName[100], TCanvas * myCanvas)
{
  myCanvas->cd(1);
  h_epp->SetLineColor(1);
  h_epp->Draw();
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    h_epp_Q2bin[i]->SetLineColor(1);
    h_epp_Q2bin[i]->Draw();
  }

}

void many_plots::Draw_hist_set_epp_same(char fileName[100], TCanvas * myCanvas, many_plots scaler)
{
  double s = scaler.getScale_epp()/h_epp->Integral();
  myCanvas->cd(1);
  h_epp->Scale(s);
  h_epp->SetLineColor(2);
  h_epp->Draw("SAME");
  for(int i=0; i<10; i++){
    myCanvas->cd(i+3);
    double si = scaler.getScale_eppQ2(i)/h_epp_Q2bin[i]->Integral();
    h_epp_Q2bin[i]->Scale(si);
    h_epp_Q2bin[i]->SetLineColor(2);
    h_epp_Q2bin[i]->Draw("SAME");
  }

}

double many_plots::getScale_ep(){
  return h_ep->Integral();  
}

double many_plots::getScale_epp(){
  return h_epp->Integral();  
}

double many_plots::getScale_epQ2(int i){
  return h_ep_Q2bin[i]->Integral();  
}

double many_plots::getScale_eppQ2(int i){
  return h_epp_Q2bin[i]->Integral();  
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
