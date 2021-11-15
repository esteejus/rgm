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
  if (argc != 8 and argc != 5)
    {
      cout << "Wrong number of arguments. To use default simulations, enter \n compare  <Ebeam (GeV)> <Z> <N> <CTOF = 0, FTOF = 1> </path/to/output/root/file> </path/to/output/pdf/file> </path/to/input/data/root/file> \n\n To directly point to an input simulation file, enter \n compare </path/to/output/root/file> </path/to/output/pdf/file> </path/to/input/data/root/file> </path/to/input/simulation/root/file>\n\n";
      cout << "Valid Target Options:\n\n";
      cout << "Target  Z   N \n";
      cout << "--------------\n";
      cout << "d       1   1 \n";
      cout << "He      2   2 \n";
      cout << "C       6   6 \n";
      cout << "Ar     18  22 \n";
      cout << "Ca40   20  20 \n";
      cout << "Ca48   20  28 \n";
      cout << "Sn     50  70 \n";
      exit(-2);
    }

  //declaring variables to contain inputs
  TFile *sim_file, *data_file,  *outfile;
  char * pdffile; 

  if (argc==8)
    {
      double E_beam = atof(argv[1]);
      int Z = atoi(argv[2]);
      int N = atoi(argv[3]);
      int det = atoi(argv[4]);


      ///////////////////////
      ////input root files///
      ///////////////////////       
      
      
      data_file = new TFile  (argv[7]); //data root file                                     
      cout << "My data file is " << data_file << "\n";                                    

      // grabbing simulation root file
      if ((E_beam == 4.) or (E_beam == 4) or (E_beam == 4.0)) //4GeV
	{
	  if ((Z == 1) and (N == 1)) //deuterium
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/d4_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/d4_epn_ftof.root");
		}
	    }
	  else if ((Z == 6) and (N == 6)) //carbon
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/c4_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/c4_epn_ftof.root");
		}
	      
	    }
	  else if ((Z == 18) and (Z == 22)) //argon
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ar4_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ar4_epn_ftof.root");
		}
	    }
	  else
	    {
	      cout << "Please choose a valid target/beam combination \n";
	    }
	}
      else if ((E_beam == 6.) or (E_beam == 6) or (E_beam == 6.0)) //6GeV
	{
	  if ((Z == 1) and (N == 1))//deuterium
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/d6_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/d6_epn_ftof.root");
		}	    
	    }
	  else if ((Z == 2) and (N == 2)) //helium
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/he6_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/he6_epn_ftof.root");
		}	    
	    }
	  else if ( ((Z == 6) and (N == 6)) or ((Z == 20) and (N == 28)) or ((Z == 50) and (N == 70)) ) //carbon, calcium 48, tin
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/c6_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/c6_epn_ftof.root");
		}	     
	    }
	  else if ((Z == 20) and (N == 20)) //calcium 40
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ca406_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  == new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ca406_epn_ftof.root");
		}	    
	    }
	  else if ((Z == 18) and (N == 22)) //argon
	    {
	      if (det == 0) //ctof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ar6_epn_ctof.root"); 
		}
	      if (det == 1) //ftof
		{
		  sim_file  = new TFile ("/u/group/clas12/users/rgm/rgm_software/Monitoring/sim_files/epn/ar6_epn_ftof.root");
		}	    
	    }
	  else
	    {
	      cout << "Please choose a valid target \n";
	    }
	}
      else 
	{
	  cout << "Please choose 4 or 6 GeV as your beam energy \n";
	}
      cout << "my sim file is " << sim_file << "\n";

      ///////////////////
      ////output files///                                                                         
      ///////////////////      
      
      outfile = new TFile (argv[5],"RECREATE");
      pdffile = argv[6];
    

    }
  else if (argc == 5)
    {
      ///////////////////////
      ////input root files///
      ///////////////////////       
      
      data_file = new TFile  (argv[3]); //data root file                                     
      cout << "My data file is " << data_file << "\n";                                    

      sim_file  = new TFile (argv[4]); //simulation  ("golden run") root file        
      cout << "My sim file is " << sim_file <<"\n";


      ///////////////////
      ////output files///                                                                       
      ///////////////////      
      
      outfile = new TFile (argv[1],"RECREATE");
      pdffile = argv[2];


    }


 
  //////////////////////////////////////////
  ///grab histograms from simulation file///                                              
  //////////////////////////////////////////

  TH1D * sim_xB_SRC = (TH1D*)sim_file->Get("xB_SRC");
  TH1D * sim_pmiss_SRC = (TH1D*)sim_file->Get("pmiss_SRC");
  TH1D * sim_mmiss_SRC = (TH1D*)sim_file->Get("mmiss_SRC");

  TH1D * sim_p_2_AllRec = (TH1D*)sim_file->Get("p_2_AllRec");
  TH1D * sim_chiSq_rec_AllRec = (TH1D*)sim_file->Get("chiSq_rec_AllRec");
  TH1D * sim_count_AllRec = (TH1D*)sim_file->Get("count_AllRec");

  TH1D * sim_p_2_Rec = (TH1D*)sim_file->Get("p_2_Rec");
  TH1D * sim_p_rel_Rec = (TH1D*)sim_file->Get("p_rel_Rec");
  TH1D * sim_p_cm_Rec = (TH1D*)sim_file->Get("p_cm_Rec");
  TH1D * sim_p_t_cm_Rec = (TH1D*)sim_file->Get("p_t_cm_Rec");
  TH1D * sim_p_y_cm_Rec =(TH1D*)sim_file->Get("p_y_cm_Rec");
  TH1D * sim_p_x_cm_Rec =(TH1D*)sim_file->Get("p_x_cm_Rec");
  TH1D * sim_theta_rel_Rec =(TH1D*)sim_file->Get("theta_rel_Rec");
  TH1D * sim_nbeta = (TH1D*)sim_file->Get("n beta");
  TH1D * sim_tofm = (TH1D*)sim_file->Get("tof_m");
  
  if (!sim_nbeta)
    {
      cout << "sim_nbeta couldn't be grabbed\n";
    }
  TH1D * sim_cos0 = (TH1D*)sim_file->Get("cos0");

  TH1D * sim_cos0_cut = (TH1D*)sim_file->Get("cos0_cut");
  if (!sim_cos0_cut)
    {
      cout << "sim_cos0_cut couldn't be grabbed\n";
      exit(-2);
    }
  ////////////////////////////////////
  ///grab histograms from data file///                             
  ////////////////////////////////////
  TH1D * data_xB_SRC = (TH1D*)data_file->Get("xB_SRC");
  TH1D * data_pmiss_SRC = (TH1D*)data_file->Get("pmiss_SRC");
  TH1D * data_mmiss_SRC = (TH1D*)data_file->Get("mmiss_SRC");

  TH1D * data_p_2_AllRec = (TH1D*)data_file->Get("p_2_AllRec");
  TH1D * data_chiSq_rec_AllRec = (TH1D*)data_file->Get("chiSq_rec_AllRec");
  TH1D * data_count_AllRec = (TH1D*)data_file->Get("count_AllRec");

  TH1D * data_p_2_Rec = (TH1D*)data_file->Get("p_2_Rec");
  TH1D * data_p_rel_Rec = (TH1D*)data_file->Get("p_rel_Rec");
  TH1D * data_p_cm_Rec =(TH1D*)data_file->Get("p_cm_Rec");
  TH1D * data_p_t_cm_Rec =(TH1D*)data_file->Get("p_t_cm_Rec");
  TH1D * data_p_y_cm_Rec =(TH1D*)data_file->Get("p_y_cm_Rec");
  TH1D * data_p_x_cm_Rec =(TH1D*)data_file->Get("p_x_cm_Rec");
  TH1D * data_theta_rel_Rec =(TH1D*)data_file->Get("theta_rel_Rec");
  TH1D * data_nbeta = (TH1D*)data_file->Get("n beta");
  TH1D * data_tofm = (TH1D*)data_file->Get("tof_m");
                                                
  if (!data_nbeta)
    {
      cout << "data_nbeta couldn't be grabbed\n";
      exit(-2);
    }                        
  TH1D * data_cos0 = (TH1D*)data_file->Get("cos0");

  TH1D * data_cos0_cut = (TH1D*)data_file->Get("cos0_cut");

  if (!data_cos0_cut)
    {
      cout << "data_cos0_cut couldn't be grabbed\n";
      exit(-2);
    }
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

  double n_bins_chiSq_rec_AllRec = sim_chiSq_rec_AllRec->GetXaxis()->GetNbins();
  double sim_chiSq_rec_AllRec_norm = sim_chiSq_rec_AllRec->Integral(1,n_bins_chiSq_rec_AllRec);
  double data_chiSq_rec_AllRec_norm = data_chiSq_rec_AllRec->Integral(1,n_bins_chiSq_rec_AllRec);

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

  double n_bins_nbeta = sim_nbeta->GetXaxis()->GetNbins();
  double sim_nbeta_norm=sim_nbeta->Integral(1,n_bins_nbeta);
  double data_nbeta_norm=data_nbeta->Integral(1,n_bins_nbeta);

  double n_bins_tofm = sim_tofm->GetXaxis()->GetNbins();
  double sim_tofm_norm=sim_tofm->Integral(1,n_bins_tofm);
  double data_tofm_norm=data_tofm->Integral(1,n_bins_tofm);

  //Deuterium Momentum Analysis
  double n_bins_cos0 = sim_cos0->GetXaxis()->GetNbins();
  double sim_cos0_norm = sim_cos0->Integral(1,n_bins_cos0);
  double data_cos0_norm = sim_cos0->Integral(1,n_bins_cos0);

  //Deuterium with cuts
  double n_bins_cos0_cut = sim_cos0_cut->GetXaxis()->GetNbins();
  double sim_cos0_cut_norm = sim_cos0_cut->Integral(1,n_bins_cos0_cut);
  double data_cos0_cut_norm = sim_cos0_cut->Integral(1,n_bins_cos0_cut);

 
  /////////////////////////////////////////////////////////////
  /////re-scaling simulation distributions to incoming data////
  /////////////////////////////////////////////////////////////
  
  //Lead SRC Proton Checks
  sim_xB_SRC->Scale(data_xB_SRC_norm/sim_xB_SRC_norm);
  sim_pmiss_SRC->Scale(data_pmiss_SRC_norm/sim_pmiss_SRC_norm);
  sim_mmiss_SRC->Scale(data_mmiss_SRC_norm/sim_mmiss_SRC_norm);

  //Recoil Nucleons
  sim_p_2_AllRec->Scale(data_p_2_AllRec_norm/sim_p_2_AllRec_norm);
  sim_chiSq_rec_AllRec->Scale(data_chiSq_rec_AllRec_norm/sim_chiSq_rec_AllRec_norm);
  sim_count_AllRec->Scale(data_count_AllRec_norm/sim_count_AllRec_norm);

  //Recoil SRC Nucleons
  sim_p_2_Rec->Scale(data_p_2_Rec_norm/sim_p_2_Rec_norm);
  sim_p_rel_Rec->Scale(data_p_rel_Rec_norm/sim_p_rel_Rec_norm);
  sim_p_cm_Rec->Scale(data_p_cm_Rec_norm/sim_p_cm_Rec_norm);
  sim_p_t_cm_Rec->Scale(data_p_t_cm_Rec_norm/sim_p_t_cm_Rec_norm);
  sim_p_y_cm_Rec->Scale(data_p_y_cm_Rec_norm/sim_p_y_cm_Rec_norm);
  sim_p_x_cm_Rec->Scale(data_p_x_cm_Rec_norm/sim_p_x_cm_Rec_norm);
  sim_theta_rel_Rec->Scale(data_theta_rel_Rec_norm/sim_theta_rel_Rec_norm);
  sim_nbeta->Scale(data_nbeta_norm/sim_nbeta_norm);
  sim_tofm->Scale(data_tofm_norm/sim_tofm_norm);

  //Deuterium Momentum Analysis
  sim_cos0->Scale(data_cos0_norm/sim_cos0_norm);
  
  //Deuterium with Cuts
  sim_cos0_cut->Scale(data_cos0_cut_norm/sim_cos0_cut_norm);


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
  TLegend *legend1 = new TLegend(0.11,0.7,0.3,0.9);
  legend1->SetTextSize(.04);
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
  TLegend *legend2 = new TLegend(0.11,0.7,0.3,0.9);
  legend2->SetTextSize(.04);
  legend2->SetHeader("Legend","C");
  legend2->AddEntry(data_p_2_AllRec,"Data","lep");
  legend2->AddEntry(sim_p_2_AllRec, "Simulation","lep");
  legend2->Draw();
  myCanvas->cd(2);
  data_chiSq_rec_AllRec->SetLineColor(1);
  data_chiSq_rec_AllRec->SetMarkerColor(1);
  data_chiSq_rec_AllRec->Draw();
  sim_chiSq_rec_AllRec->SetLineColor(2);
  sim_chiSq_rec_AllRec->SetMarkerColor(2);
  sim_chiSq_rec_AllRec->Draw("SAME");
  myCanvas->cd(3);
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
  TLegend *legend3 = new TLegend(0.11,0.7,0.3,0.9);
  legend3->SetTextSize(.04);
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

  myCanvas->Divide(3,2);
  myCanvas->cd(1);
  data_p_y_cm_Rec->SetLineColor(1);
  data_p_y_cm_Rec->SetMarkerColor(1);
  data_p_y_cm_Rec->Draw();
  sim_p_y_cm_Rec->SetLineColor(2);
  sim_p_y_cm_Rec->SetMarkerColor(2);
  sim_p_y_cm_Rec->Draw("SAME");
  TLegend *legend4 = new TLegend(0.11,0.7,0.3,0.9);
  legend4->SetTextSize(.04);
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
  myCanvas->cd(4);
  data_nbeta->SetLineColor(1);
  data_nbeta->SetMarkerColor(1);
  data_nbeta->Draw();
  sim_nbeta->SetLineColor(1);
  sim_nbeta->SetMarkerColor(1);
  sim_nbeta->Draw("SAME");
  myCanvas->cd(5);
  data_tofm->SetLineColor(1);
  data_tofm->SetMarkerColor(1);
  data_tofm->Draw();
  sim_tofm->SetLineColor(1);
  sim_tofm->SetMarkerColor(1);
  sim_tofm->Draw("SAME");  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  //Deuterium only                                                                                 
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts");
  text.DrawLatex(0.2,0.8,"Deuterium only");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_cos0->SetLineColor(1);
  data_cos0->SetMarkerColor(1);
  data_cos0->Draw();
  sim_cos0->SetLineColor(2);
  sim_cos0->SetMarkerColor(2);
  sim_cos0->Draw("SAME");
  TLegend *legend5 = new TLegend(0.11,0.7,0.3,0.9);
  legend5->SetTextSize(.04);
  legend5->SetHeader("Legend","C");
  legend5->AddEntry(data_p_2_Rec,"Data","lep");
  legend5->AddEntry(sim_p_2_Rec, "Simulation","lep");
  legend5->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //Deuterium with cuts                                                                            
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts");
  text.DrawLatex(0.2,0.8,"Deuterium only");
  text.DrawLatex(0.2,0.7,"cos(#theta_{pmiss,pneutron}>0.95");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  data_cos0_cut->SetLineColor(1);
  data_cos0_cut->SetMarkerColor(1);
  data_cos0_cut->Draw();
  sim_cos0_cut->SetLineColor(2);
  sim_cos0_cut->SetMarkerColor(2);
  sim_cos0_cut->Draw("SAME");
  TLegend *legend6 = new TLegend(0.11,0.7,0.3,0.9);
  legend6->SetTextSize(.04);
  legend6->SetHeader("Legend","C");
  legend6->AddEntry(data_p_2_Rec,"Data","lep");
  legend6->AddEntry(sim_p_2_Rec, "Simulation","lep");
  legend6->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  sprintf(fileName,"%s]",pdffile);
  myCanvas->Print(fileName,"pdf");

  outfile->Close();

    

}
