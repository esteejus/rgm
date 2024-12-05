#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "reweighter.h"
#include "TGraphErrors.h"


using namespace std;
using namespace clas12;

const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mass_pi = db->GetParticle(-211)->Mass();
double mD = 1.8756;

double beam_E = 5.98;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

//Parameters
double A[2]= {-0.000695124,-0.000355869};
double B[2]= {0.00182181,7.77933e-05};
double C[3]= {0.000266541,0.424055,49.068};



void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double func(double x, double a, double b, double c){
  return a + b*x + (c/x); 
}

double funcAB(double x, double a, double b){
  return a + b*x; 
}

double funcC(double x, double a, double b, double c){
  return a - b/(x-c); 
}

double Corr(double mom, double theta){
  double a = funcAB(theta,A[0],A[1]);
  double b = funcAB(theta,B[0],B[1]);
  double c = funcC(theta,C[0],C[1],C[2]);
  double correction = func(mom,a,b,c);
  return correction;
}

double getExp(TLorentzVector balance_ptr, TLorentzVector par){
  double theta_bpar = balance_ptr.Vect().Angle(par.Vect());
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_pi*mass_pi) - balance_ptr.M2() - par.M2()) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*par.M2();
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char fileName[100]){
  int ctr = 0;
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    //proj->Rebin(2);
    //Now preform a guassian fit
    if(proj->GetEntries()<15){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.3,0.3,3);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,mode);
    gFit->SetParLimits(1,mode-0.025,mode+0.025);
    gFit->SetParameter(2,0.01);
    gFit->SetParLimits(2,0.001,0.1);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",mode-0.05,0.05);

    /*
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    proj->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
    */

    if((gPoint == 0)){
      if(gPoint->Parameter(1)<0.07){
	g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
	g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
      }
    }
    proj->Write();
  }
}

void getFunctionMomTFRP(TH2D * h_myhist, TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double min, double max,TCanvas * myCanvas, char fileName[100]){
  getGraph(h_myhist,g_mygraph,myCanvas,fileName);

  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-0.2,+0.2);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,0,+0.01);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-0.5,+0.05);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",min,max);
      
}

void getABC(TH2D * h_myhist[16], TGraphErrors * g_mygraph[16], TF1 * f_myfunc[16], TFitResultPtr p_mypoint[16], double min, double max, TGraph * g_Pargraph[3], TF1 * f_Parfunc[3], TFitResultPtr p_Parpoint[3], TF1 * f_Combfunc[16], TCanvas * myCanvas, char fileName[100]){
  
  for(int i=0; i<16; i++){
    double theta = 4.25 + i*2.5;
    getFunctionMomTFRP(h_myhist[i],g_mygraph[i],f_myfunc[i],p_mypoint[i],min,max,myCanvas,fileName);
    for(int j=0; j<3; j++){
      if(p_mypoint[i]!=-1){
	g_Pargraph[j]->SetPoint(g_Pargraph[j]->GetN(),theta,p_mypoint[i]->Parameter(j));
      }
    }
  }
  
  for(int j=0; j<2; j++){
    f_Parfunc[j]->SetLineColor(4);
    f_Parfunc[j]->SetParameter(0,0);
    f_Parfunc[j]->SetParLimits(0,-0.01,0.01);
    f_Parfunc[j]->SetParameter(1,0);
    f_Parfunc[j]->SetParLimits(1,-0.01,0.01);
    p_Parpoint[j] = g_Pargraph[j]->Fit(f_Parfunc[j],"SrBeqn","",5,41);
  }  
  f_Parfunc[2]->SetLineColor(4);
  f_Parfunc[2]->SetParameter(0,0);
  f_Parfunc[2]->SetParLimits(0,-0.1,0.1);
  f_Parfunc[2]->SetParameter(1,0.01);
  f_Parfunc[2]->SetParLimits(1,0,10);
  f_Parfunc[2]->SetParameter(2,60);
  f_Parfunc[2]->SetParLimits(2,45,100);
  p_Parpoint[2] = g_Pargraph[2]->Fit(f_Parfunc[2],"SrBeqn","",5,41);

  for(int i=0; i<16; i++){
    double theta = 4.25 + i*2.5;
    f_Combfunc[i]->SetLineColor(4);
    f_Combfunc[i]->SetLineWidth(1);
    f_Combfunc[i]->SetParameter(0,f_Parfunc[0]->Eval(theta));
    f_Combfunc[i]->SetParameter(1,f_Parfunc[1]->Eval(theta));
    f_Combfunc[i]->SetParameter(2,f_Parfunc[2]->Eval(theta));
  }
  

  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
    myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    f_myfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(4,4);
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);    
    h_myhist[i]->Draw("colz");
    g_mygraph[i]->Draw("SAME");
    f_myfunc[i]->Draw("SAME");
    f_Combfunc[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  myCanvas->Divide(2,2);
  for(int j=0; j<3; j++){
    myCanvas->cd(j+1);    
    g_Pargraph[j]->Draw();
    f_Parfunc[j]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
    
}


void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  int isMC = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target_ptr(0,0,0,mD);
  //TLorentzVector target_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr_p1(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr_p2(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pim_ptr(0,0,0,db->GetParticle(-211)->Mass());
  reweighter newWeight(beam_E,6,6,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;


  //////////////////////////////////
  //Proton FD
  //////////////////////////////////
  TH2D * h_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TGraphErrors * g_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TF1 * f_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TFitResultPtr p_pFDmom_pFDDeltaP_thetaGroup[6][16];
  TGraph * g_pFDmom_pFDDeltaP_Pars[6][3];
  TF1 * f_pFDmom_pFDDeltaP_Pars[6][3];
  TFitResultPtr p_pFDmom_pFDDeltaP_Pars[6][3];
  TF1 * f_pFDmom_pFDDeltaP_combined_thetaGroup[6][16];
  TH2D * h_pFDmom_pFDDeltaP_corrected_thetaGroup[6][16];
  for(int k=0; k<6; k++){
    for(int i=0; i<16; i++){
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TH2D(temp_name,temp_name,25,0,4.5,200,-0.2,0.2);
      hist_list.push_back(h_pFDmom_pFDDeltaP_thetaGroup[k][i]);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetLineColor(2);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_combined_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_corrected_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_corrected_thetaGroup[k][i] = new TH2D(temp_name,temp_name,25,0,4.5,100,-0.2,0.2);
      hist_list.push_back(h_pFDmom_pFDDeltaP_corrected_thetaGroup[k][i]);
    }
    for(int j=0; j<3; j++){
      g_pFDmom_pFDDeltaP_Pars[k][j] = new TGraph();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_Pars_%d",k+1,j+1);
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetLineColor(3);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_Pars_%d",k+1,j+1);
      if(j<2){
	f_pFDmom_pFDDeltaP_Pars[k][j] = new TF1(temp_name,[&](double *x, double *p){ return funcAB(x[0],p[0],p[1]); },5,41,2);
      }
      else{
	f_pFDmom_pFDDeltaP_Pars[k][j] = new TF1(temp_name,[&](double *x, double *p){ return funcC(x[0],p[0],p[1],p[2]); },5,41,3);
      }
    }

    for(int i=0; i<16; i++){
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      h_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TH2D(temp_name,temp_name,50,0,4.5,50,-0.2,0.2);
      hist_list.push_back(h_pFDmom_pFDDeltaP_thetaGroup[k][i]);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_thetaGroup[k][i]->SetLineColor(2);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",k+1,i+1);
      f_pFDmom_pFDDeltaP_thetaGroup[k][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
      
    }
  }
  
  TH2D * h_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TGraphErrors * g_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TF1 * f_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_thetaGroup[16];
  TGraph * g_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_Pars[3];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_combined_thetaGroup[16];
  TH2D * h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[16];
  for(int i=0; i<16; i++){
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TH2D(temp_name,temp_name,25,0,4.5,200,-0.2,0.2);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_thetaGroup[i]);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TGraphErrors();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetLineColor(2);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_combined_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_corrected_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i] = new TH2D(temp_name,temp_name,25,0,4.5,200,-0.2,0.2);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]);
  }
  for(int j=0; j<3; j++){
    g_pFDmom_pFDDeltaP_int_Pars[j] = new TGraph();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetLineColor(3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    if(j<2){
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcAB(x[0],p[0],p[1]); },5,41,2);
    }
    else{
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcC(x[0],p[0],p[1],p[2]); },5,41,3);
    }
  }
  

  TH2D * h_pFDtheta_pFDDeltaP[6];
  TGraphErrors * g_pFDtheta_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDtheta_pFDDeltaP_sector_%d",j+1);
    h_pFDtheta_pFDDeltaP[j] = new TH2D(temp_name,temp_name,50,0,50,50,-0.2,0.2);
    hist_list.push_back(h_pFDtheta_pFDDeltaP[j]);
    g_pFDtheta_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDtheta_pFDDeltaP_sector_%d",j+1);
    g_pFDtheta_pFDDeltaP[j]->SetName(temp_name);
    g_pFDtheta_pFDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pFDphi_pFDDeltaP[6];
  TGraphErrors * g_pFDphi_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDphi_pFDDeltaP_sector_%d",j+1);
    h_pFDphi_pFDDeltaP[j] = new TH2D(temp_name,temp_name,50,-35,35,50,-0.2,0.2);
    hist_list.push_back(h_pFDphi_pFDDeltaP[j]);
    g_pFDphi_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDphi_pFDDeltaP_sector_%d",j+1);
    g_pFDphi_pFDDeltaP[j]->SetName(temp_name);
    g_pFDphi_pFDDeltaP[j]->SetLineColor(2);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next() && counter <100000000000)
    {

      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight;// * newWeight.get_weight_ep(c12->mcparts());
	wepp = original_weight;// * newWeight.get_weight_epp(c12->mcparts());
      }

      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pims = clasAna.getByPid(-211);
      if(electrons.size() == 1)
	{

	  SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double vtz_e = electrons[0]->par()->getVz();
	  
	  int sector_e = electrons[0]->getSector()-1;
	  double mom_e = el.P();
	  double theta_e = el.Theta() * 180 / M_PI;
	  double phi_e = el.Phi() * 180 / M_PI;
	  double shift_e = 7.5;
	  shift_e += (sector_e==0)?0:(sector_e==1)?60:(sector_e==2)?120:(sector_e==3 && phi_e>0)?180:(sector_e==3 && phi_e<0)?-180:(sector_e==4)?-120:(sector_e==5)?-60:0;
	  phi_e -= shift_e;


	  if(vtz_e<-5.5){continue;}
	  if(vtz_e>0){continue;}

	  TVector3 pmc_e(c12->mcparts()->getPx(0),c12->mcparts()->getPy(0),c12->mcparts()->getPz(0));
	  TVector3 pmc_p(c12->mcparts()->getPx(1),c12->mcparts()->getPy(1),c12->mcparts()->getPz(1));

	  double eDeltaP = pmc_e.Mag() - el.P();
	  int thetaeGroup = (theta_e - 10.0)/3.0;
	  //Fill histograms
	  if(protons.size() != 1){continue;}
	  int index_p1 = 0;
	  SetLorentzVector(proton_ptr_p1,protons[index_p1]);
	  double beta_p1 = protons[index_p1]->par()->getBeta();
	  double mom_p1 = proton_ptr_p1.P();
	  double momT_p1 = proton_ptr_p1.Vect().Perp();
	  double theta_p1 = proton_ptr_p1.Theta() * 180 / M_PI;
	  double phi_p1 = proton_ptr_p1.Phi() * 180 / M_PI;
	  double vtz_p1 = protons[index_p1]->par()->getVz();
	  if(fabs(vtz_e-vtz_p1)>1.5){continue;}
	  double p1DeltaP = pmc_p.Mag() - proton_ptr_p1.P();
	  
	  if(protons[index_p1]->getRegion()!=FD){continue;}
	  int sector_p1 = protons[index_p1]->getSector()-1;
	  double shift_p1 = -22;
	  int thetapFDGroup = (theta_p1 - 3.0)/2.5;
	  double correction = Corr(mom_p1,theta_p1);
	  
	  shift_p1 += (sector_p1==0)?0:(sector_p1==1)?60:(sector_p1==2)?120:(sector_p1==3 && phi_p1>0)?180:(sector_p1==3 && phi_p1<0)?-180:(sector_p1==4)?-120:(sector_p1==5)?-60:0;
	  phi_p1 -= shift_p1;
	  //Proton FD

	  if((thetapFDGroup>=0) && (thetapFDGroup<16)){
	    h_pFDmom_pFDDeltaP_thetaGroup[sector_p1][thetapFDGroup]->Fill(mom_p1,p1DeltaP,wep);	      
	    h_pFDmom_pFDDeltaP_corrected_thetaGroup[sector_p1][thetapFDGroup]->Fill(mom_p1+correction,p1DeltaP-correction,wep);	      
	    h_pFDmom_pFDDeltaP_int_thetaGroup[thetapFDGroup]->Fill(mom_p1,p1DeltaP,wep);
	  }
	  
	  h_pFDtheta_pFDDeltaP[sector_p1]->Fill(theta_p1,p1DeltaP,wep);
	  h_pFDphi_pFDDeltaP[sector_p1]->Fill(phi_p1,p1DeltaP,wep);
	  
	}
    }
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);
  
  getABC(h_pFDmom_pFDDeltaP_int_thetaGroup,g_pFDmom_pFDDeltaP_int_thetaGroup,f_pFDmom_pFDDeltaP_int_thetaGroup,p_pFDmom_pFDDeltaP_int_thetaGroup,0.3,5.5,g_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_Pars,p_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_combined_thetaGroup,myCanvas,fileName);

  for(int k = 0; k < 6; k++){
    getABC(h_pFDmom_pFDDeltaP_thetaGroup[k],g_pFDmom_pFDDeltaP_thetaGroup[k],f_pFDmom_pFDDeltaP_thetaGroup[k],p_pFDmom_pFDDeltaP_thetaGroup[k],0.3,5.5,g_pFDmom_pFDDeltaP_Pars[k],f_pFDmom_pFDDeltaP_Pars[k],p_pFDmom_pFDDeltaP_Pars[k],f_pFDmom_pFDDeltaP_combined_thetaGroup[k],myCanvas,fileName);
  }
  
  myCanvas->Divide(2,2);
  for(int j=0; j<3; j++){
    myCanvas->cd(j+1);    
    g_pFDmom_pFDDeltaP_int_Pars[j]->Draw();
    for(int k=0; k<6; k++){
      g_pFDmom_pFDDeltaP_Pars[k][j]->SetLineColor(k+2);
      g_pFDmom_pFDDeltaP_Pars[k][j]->Draw("SAME");
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(3,3);  
  for(int i=0; i<16; i++){
    myCanvas->cd(i+1);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i]->SetLineWidth(3);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i]->Draw();
    for(int k=0; k<6; k++){
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->SetLineWidth(1);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->SetLineColor(k+2);
      f_pFDmom_pFDDeltaP_combined_thetaGroup[k][i]->Draw("SAME");
    }
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  for(int sector = 0; sector<6; sector++){
    myCanvas->Divide(4,4);
    for(int i=0; i<16; i++){
      myCanvas->cd(i+1);    
      h_pFDmom_pFDDeltaP_corrected_thetaGroup[sector][i]->Draw("colz");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear(); 
  }
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  cout<<"A=["<< p_pFDmom_pFDDeltaP_int_Pars[0]->Parameter(0)
      <<","<< p_pFDmom_pFDDeltaP_int_Pars[0]->Parameter(1)
      <<"]"<<endl
      <<"B=["<< p_pFDmom_pFDDeltaP_int_Pars[1]->Parameter(0)
      <<","<< p_pFDmom_pFDDeltaP_int_Pars[1]->Parameter(1)
      <<"]"<<endl
      <<"C=["<< p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(0)
      <<","<< p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(1)
      <<","<< p_pFDmom_pFDDeltaP_int_Pars[2]->Parameter(2)
      <<"]"<<endl;
  
  
  f->Close();

  return 0;
}

/*
	  TLorentzVector balance_ptr = beam + target_ptr - proton_ptr_p1 - proton_ptr_p2;
	  double theta_be = balance_ptr.Vect().Angle(el.Vect());
	  double Eb = balance_ptr.E();
	  double Pb = balance_ptr.P();
	  double K = (pim_ptr.M2() - balance_ptr.M2() - el.M2()) / 2;
	  double a = Pb*Pb*cos(theta_be)*cos(theta_be) - Eb*Eb;
	  double b = -2 * K * Pb * cos(theta_be);
	  double c = K*K - Eb*Eb*el.M2();
	  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
o	  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
*/
