#include <cstdlib>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;

int main(int argc, char ** argv)
{
  if (argc!=5)
    {
      cout << "Wrong number of arguments. Instead use \n compare </path/to/output/root/file> </path/to/output/pdf/file> </path/to/input/simulation/root/file> </path/to/input/data/root/file> \n";
      exit(-2);
    }


  ///////////////////////
  ////input root files///                                                                            ///////////////////////       
                                  
  TFile * sim_file  = new TFile (argv[3]); //simulation  ("golden run") root file                                                      
  TFile * data_file = new TFile  (argv[4]); //data root file                                                                           
  //////////////////////////////////////////
  ///grab histograms from simulation file///                                                         //////////////////////////////////////////

  TH1D * sim_xB_SRC = (TH1D*)sim_file->Get("xB_SRC");
  TH1D * sim_pmiss_SRC = (TH1D*)sim_file->Get("pmiss_SRC");
  TH1D * sim_mmiss_SRC = (TH1D*)sim_file->Get("mmiss_SRC");

  TH1D * sim_p_2_AllRec = (TH1D*)sim_file->Get("p_2_AllRec");
  TH1D * sim_count_AllRec = (TH1D*)sim_file->Get("count_AllRec");

  TH1D * sim_p_2_Rec = (TH1D*)sim_file->Get("p_2_Rec");
  TH1D * sim_p_rel_Rec = (TH1D*)sim_file->Get("p_rel_Rec");
  TH1D * sim_p_cm_Rec = (TH1D*)sim_file->Get("p_cm_Rec");
  TH1D * sim_p_t_cm_Rec = (TH1D*)sim_file->Get("p_t_cm_Rec");
  TH1D * sim_p_y_cm_Rec =(TH1D*)sim_file->Get("p_y_cm_Rec");
  TH1D * sim_p_x_cm_Rec =(TH1D*)sim_file->Get("p_x_cm_Rec");
  TH1D * sim_theta_rel_Rec =(TH1D*)sim_file->Get("theta_rel_Rec");

  ////////////////////////////////////
  ///grab histograms from data file///                             
  ////////////////////////////////////
  TH1D * data_xB_SRC = (TH1D*)data_file->Get("xB_SRC");
  TH1D * data_pmiss_SRC = (TH1D*)data_file->Get("pmiss_SRC");
  TH1D * data_mmiss_SRC = (TH1D*)data_file->Get("mmiss_SRC");

  TH1D * data_p_2_AllRec = (TH1D*)data_file->Get("p_2_AllRec");
  TH1D * data_count_AllRec = (TH1D*)data_file->Get("count_AllRec");

  TH1D * data_p_2_Rec = (TH1D*)data_file->Get("p_2_Rec");
  TH1D * data_p_rel_Rec = (TH1D*)data_file->Get("p_rel_Rec");
  TH1D * data_p_cm_Rec =(TH1D*)data_file->Get("p_cm_Rec");
  TH1D * data_p_t_cm_Rec =(TH1D*)data_file->Get("p_t_cm_Rec");
  TH1D * data_p_y_cm_Rec =(TH1D*)data_file->Get("p_y_cm_Rec");
  TH1D * data_p_x_cm_Rec =(TH1D*)data_file->Get("p_x_cm_Rec");
  TH1D * data_theta_rel_Rec =(TH1D*)data_file->Get("theta_rel_Rec");
                                                                        

  ///////////////////
  ////output files///                                                                                ///////////////////      
                               
  TFile * outfile = new TFile (argv[1],"RECREATE");
  char * pdffile = argv[2];

  
  ////////////////////////////////////////
  ////calculating normalization factors///
  ////////////////////////////////////////

  //Lead SRC Checks
  double n_bins_xB_SRC = sim_xB_SRC->GetXaxis()->GetNbins();
  double sim_xB_SRC_norm = sim_xB_SRC->Integral(1,n_bins_xB_SRC);
  double data_xB_SRC_norm = data_xB_SRC->Integral(1,n_bins_xB_SRC);

  double n_bins_pmiss_SRC = sim_pmiss_SRC->GetXaxis()->GetNbins();
  double sim_pmiss_SRC_norm = sim_pmiss_SRC->Integral(1,n_bins_pmiss_SRC);
  double data_pmiss_SRC_norm = data_pmiss_SRC->Integral(1,n_bins_pmiss_SRC);

  double n_bins_mmiss_SRC = sim_mmiss_SRC->GetXaxis()->GetNbins();
  double sim_mmiss_SRC_norm = sim_mmiss_SRC->Integral(1,n_bins_mmiss_SRC);
  double data_mmiss_SRC_norm = data_mmiss_SRC->Integral(1,n_bins_mmiss_SRC);

  //Recoil Nucleons
  double n_bins_p_2_AllRec = sim_p_2_AllRec->GetXaxis()->GetNbins();
  double sim_p_2_AllRec_norm = sim_p_2_AllRec->Integral(1,n_bins_p_2_AllRec);
  double data_p_2_AllRec_norm = data_p_2_AllRec->Integral(1,n_bins_p_2_AllRec);

  double n_bins_count_AllRec = sim_count_AllRec->GetXaxis()->GetNbins();
  double sim_count_AllRec_norm = sim_count_AllRec->Integral(1,n_bins_count_AllRec);
  double data_count_AllRec_norm = data_count_AllRec->Integral(1,n_bins_count_AllRec);


  //Recoil SRC Nucleons
  double n_bins_p_2_Rec = sim_p_2_Rec->GetXaxis()->GetNbins();
  double sim_p_2_Rec_norm = sim_p_2_Rec->Integral(1,n_bins_p_2_Rec);
  double data_p_2_Rec_norm = data_p_2_Rec->Integral(1,n_bins_p_2_Rec);

  double n_bins_p_rel_Rec = sim_p_rel_Rec->GetXaxis()->GetNbins();
  double sim_p_rel_Rec_norm = sim_p_rel_Rec->Integral(1,n_bins_p_rel_Rec);
  double data_p_rel_Rec_norm = data_p_rel_Rec->Integral(1,n_bins_p_rel_Rec);

  double n_bins_p_cm_Rec = sim_p_cm_Rec->GetXaxis()->GetNbins();
  double sim_p_cm_Rec_norm=sim_p_cm_Rec->Integral(1,n_bins_p_cm_Rec);
  double data_p_cm_Rec_norm=data_p_cm_Rec->Integral(1,n_bins_p_cm_Rec);

  double n_bins_p_t_cm_Rec = sim_p_t_cm_Rec->GetXaxis()->GetNbins();
  double sim_p_t_cm_Rec_norm=sim_p_t_cm_Rec->Integral(1,n_bins_p_t_cm_Rec);
  double data_p_t_cm_Rec_norm=data_p_t_cm_Rec->Integral(1,n_bins_p_t_cm_Rec);

  double n_bins_p_y_cm_Rec = sim_p_y_cm_Rec->GetXaxis()->GetNbins();
  double sim_p_y_cm_Rec_norm=sim_p_y_cm_Rec->Integral(1,n_bins_p_y_cm_Rec);
  double data_p_y_cm_Rec_norm=data_p_y_cm_Rec->Integral(1,n_bins_p_y_cm_Rec);

  double n_bins_p_x_cm_Rec = sim_p_x_cm_Rec->GetXaxis()->GetNbins();
  double sim_p_x_cm_Rec_norm=sim_p_x_cm_Rec->Integral(1,n_bins_p_x_cm_Rec);
  double data_p_x_cm_Rec_norm=data_p_x_cm_Rec->Integral(1,n_bins_p_x_cm_Rec);

  double n_bins_theta_rel_Rec = sim_theta_rel_Rec->GetXaxis()->GetNbins();
  double sim_theta_rel_Rec_norm=sim_theta_rel_Rec->Integral(1,n_bins_theta_rel_Rec);
  double data_theta_rel_Rec_norm=data_theta_rel_Rec->Integral(1,n_bins_theta_rel_Rec);

 
  /////////////////////////////////////////////////////////////
  /////re-scaling simulation distributions to incoming data////
  /////////////////////////////////////////////////////////////
  
  //Lead SRC Proton Checks
  sim_xB_SRC->Scale(data_xB_SRC_norm/sim_xB_SRC_norm);
  sim_pmiss_SRC->Scale(data_pmiss_SRC_norm/sim_pmiss_SRC_norm);
  sim_mmiss_SRC->Scale(data_mmiss_SRC_norm/sim_mmiss_SRC_norm);

  //Recoil Nucleons
  sim_p_2_AllRec->Scale(data_p_2_AllRec_norm/sim_p_2_AllRec_norm);
  sim_count_AllRec->Scale(data_count_AllRec_norm/sim_count_AllRec_norm);

  //Recoil SRC Nucleons
  sim_p_2_Rec->Scale(data_p_2_Rec_norm/sim_p_2_Rec_norm);
  sim_p_rel_Rec->Scale(data_p_rel_Rec_norm/sim_p_rel_Rec_norm);
  sim_p_cm_Rec->Scale(data_p_cm_Rec_norm/sim_p_cm_Rec_norm);
  sim_p_t_cm_Rec->Scale(data_p_t_cm_Rec_norm/sim_p_t_cm_Rec_norm);
  sim_p_y_cm_Rec->Scale(data_p_y_cm_Rec_norm/sim_p_y_cm_Rec_norm);
  sim_p_x_cm_Rec->Scale(data_p_x_cm_Rec_norm/sim_p_x_cm_Rec_norm);
  sim_theta_rel_Rec->Scale(data_theta_rel_Rec_norm/sim_theta_rel_Rec_norm);

  /////////////////////////////
  ////create the output PDF////
  /////////////////////////////


  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);

  char fileName[100];
  sprintf(fileName,"%s[",pdffile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdffile);

  
  /////////////////////////////
  /////Draw Histograms/////////
  /////////////////////////////



  //Lead SRC Proton Checks
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD}) Cuts");
  text.DrawLatex(0.2,0.7,"1.5 < Q^{2} [GeV]");
  text.DrawLatex(0.2,0.6,"0.3 [GeV] < p_{miss}");
  text.DrawLatex(0.2,0.5,"0.84 [GeV] < m_{mmiss} < 1.04 [GeV]");
  text.DrawLatex(0.2,0.4,"0.62 < |p|/|q| < 0.96");
  text.DrawLatex(0.2,0.3,"1.2 < x_{B}");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_xB_SRC->SetLineColor(1);
  data_xB_SRC->SetMarkerColor(1);
  data_xB_SRC->Draw();
  sim_xB_SRC->SetLineColor(2);
  sim_xB_SRC->SetMarkerColor(2);
  sim_xB_SRC->Draw("SAME");
  TLegend *legend1 = new TLegend(0.11,0.8,0.3,0.9);
  legend1->SetTextSize(.03);
  legend1->SetHeader("Legend","C");
  legend1->AddEntry(data_xB_SRC,"Data","lep");
  legend1->AddEntry(sim_xB_SRC, "Simulation","lep");
  legend1->Draw();
  myCanvas->cd(2);
  data_pmiss_SRC->SetLineColor(1);
  data_pmiss_SRC->SetMarkerColor(1); 
  data_pmiss_SRC->Draw();
  sim_pmiss_SRC->SetLineColor(2);
  sim_pmiss_SRC->SetMarkerColor(2);
  sim_pmiss_SRC->Draw("SAME");
  myCanvas->cd(3);
  data_mmiss_SRC->SetLineColor(1);
  data_mmiss_SRC->SetMarkerColor(1);
  data_mmiss_SRC->Draw();
  sim_mmiss_SRC->SetLineColor(2);
  sim_mmiss_SRC->SetMarkerColor(2);
  sim_mmiss_SRC->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  //Recoil Nucleons
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC}) Cuts");
  text.DrawLatex(0.2,0.7,"Second Proton Detected");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_p_2_AllRec->SetLineColor(1);
  data_p_2_AllRec->SetMarkerColor(1);
  data_p_2_AllRec->Draw();
  sim_p_2_AllRec->SetLineColor(2);
  sim_p_2_AllRec->SetMarkerColor(2);
  sim_p_2_AllRec->Draw("SAME");
  TLegend *legend2 = new TLegend(0.11,0.8,0.3,0.9);
  legend2->SetTextSize(.03);
  legend2->SetHeader("Legend","C");
  legend2->AddEntry(data_p_2_AllRec,"Data","lep");
  legend2->AddEntry(sim_p_2_AllRec, "Simulation","lep");
  legend2->Draw();
  myCanvas->cd(2);
  data_count_AllRec->SetLineColor(1);
  data_count_AllRec->SetMarkerColor(1);
  data_count_AllRec->Draw();
  sim_count_AllRec->SetLineColor(2);
  sim_count_AllRec->SetMarkerColor(2);
  sim_count_AllRec->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  //Recoil SRC Nucleons

  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC},p_{Rec}) Cuts");
  text.DrawLatex(0.2,0.7,"0.35 [GeV] < p_{Rec}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_p_2_Rec->SetLineColor(1);
  data_p_2_Rec->SetMarkerColor(1);
  data_p_2_Rec->Draw();
  sim_p_2_Rec->SetLineColor(2);
  sim_p_2_Rec->SetMarkerColor(2);
  sim_p_2_Rec->Draw("SAME");
  TLegend *legend3 = new TLegend(0.11,0.8,0.3,0.9);
  legend3->SetTextSize(.03);
  legend3->SetHeader("Legend","C");
  legend3->AddEntry(data_p_2_Rec,"Data","lep");
  legend3->AddEntry(sim_p_2_Rec, "Simulation","lep");
  legend3->Draw();
  myCanvas->cd(2);
  data_p_rel_Rec->SetLineColor(1);
  data_p_rel_Rec->SetMarkerColor(1);
  data_p_rel_Rec->Draw();
  sim_p_rel_Rec->SetLineColor(2);
  sim_p_rel_Rec->SetMarkerColor(2);
  sim_p_rel_Rec->Draw("SAME");
  myCanvas->cd(3);
  data_p_cm_Rec->SetLineColor(1);
  data_p_cm_Rec->SetMarkerColor(1);
  data_p_cm_Rec->Draw();
  sim_p_cm_Rec->SetLineColor(2);
  sim_p_cm_Rec->SetMarkerColor(2);
  sim_p_cm_Rec->Draw("SAME");
  myCanvas->cd(4);
  data_p_t_cm_Rec->SetLineColor(1);
  data_p_t_cm_Rec->SetMarkerColor(1);
  data_p_t_cm_Rec->Draw();
  sim_p_t_cm_Rec->SetLineColor(2);
  sim_p_t_cm_Rec->SetMarkerColor(2);
  sim_p_t_cm_Rec->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_p_y_cm_Rec->SetLineColor(1);
  data_p_y_cm_Rec->SetMarkerColor(1);
  data_p_y_cm_Rec->Draw();
  sim_p_y_cm_Rec->SetLineColor(2);
  sim_p_y_cm_Rec->SetMarkerColor(2);
  sim_p_y_cm_Rec->Draw("SAME");
  TLegend *legend4 = new TLegend(0.11,0.8,0.3,0.9);
  legend4->SetTextSize(.03);
  legend4->SetHeader("Legend","C");
  legend4->AddEntry(data_p_y_cm_Rec,"Data","lep");
  legend4->AddEntry(sim_p_y_cm_Rec, "Simulation","lep");
  legend4->Draw();
  myCanvas->cd(2);
  data_p_x_cm_Rec->SetLineColor(1);
  data_p_x_cm_Rec->SetMarkerColor(1);
  data_p_x_cm_Rec->Draw();
  sim_p_x_cm_Rec->SetLineColor(2);
  sim_p_x_cm_Rec->SetMarkerColor(2);
  sim_p_x_cm_Rec->Draw("SAME");
  myCanvas->cd(3);
  data_theta_rel_Rec->SetLineColor(1);
  data_theta_rel_Rec->SetMarkerColor(1);
  data_theta_rel_Rec->Draw();
  sim_theta_rel_Rec->SetLineColor(2);
  sim_theta_rel_Rec->SetMarkerColor(2);
  sim_theta_rel_Rec->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdffile);
  myCanvas->Print(fileName,"pdf");
  
  outfile->Close();

    

}
