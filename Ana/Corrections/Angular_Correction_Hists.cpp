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
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "reweighter.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,45};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD = 1.8756;

//double beam_E = 5.98636;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double FuncPhiDependenceFD(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

double FuncPhiDependenceCD(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C); 
}

double FuncThetaDependenceCD(double x, double A, double B){
  return A + B/(x); 
}

double FuncThetaDependenceCD2(double x, double A){
  return A ; 
}

double ErfFunc(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}


double cotan(double x){ return cos(x)/sin(x);}

double ThetaE(double ThetaP){
  double ThetaP_rad = ThetaP*M_PI/180;
  return (180/M_PI) * 2*atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaP_rad) );
}


double ThetaP(double ThetaE){
  double ThetaE_rad = ThetaE*M_PI/180;
  return (180/M_PI) * atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaE_rad/2) );
}

double SQ(double x){ return x*x;}

double ClosestPoint(double x, double ThetaEp, double ThetaPp){
  return sqrt(SQ(x-ThetaEp) + SQ(ThetaP(x)-ThetaPp));
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph){
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
    if(proj->GetEntries()<20){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-1,1,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,0.0);
    gFit->SetParLimits(1,-0.8,0.8);
    gFit->SetParameter(2,0.2);
    gFit->SetParLimits(2,0.001,1.0);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-1,1);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    proj->Write();
  }
}

void getFunction(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double guess[4],int i, int sector){
  double min = g_mygraph->GetMinimum();
  double max = g_mygraph->GetMaximum();
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,guess[0]);
  f_myfunc->SetParLimits(0,guess[0]-0.25,guess[0]+0.25);
  f_myfunc->SetParameter(1,guess[1]);
  f_myfunc->SetParLimits(1,guess[1]*0.3,guess[1]*1.0);
  f_myfunc->SetParameter(2,guess[2]);
  double min2 = (guess[2]>20)?guess[2]-10:10;
  //double max2 = ((i==0)&&(sector==4))?40:((i==0)&&(sector==5))?40:(guess[2]<540)?guess[2]+350:720;
  double max2 = ((i==0)&&(sector==4))?40:((i==0)&&(sector==5))?40:(guess[2]<180)?guess[2]+180:360;
  f_myfunc->SetParLimits(2,min2,max2);
  f_myfunc->SetParameter(3,guess[3]);
  if(i==0){
    f_myfunc->SetParLimits(3,-M_PI,M_PI);
  }
  else{
    f_myfunc->SetParLimits(3,guess[3]-(M_PI/2),guess[3]+(M_PI/2));
  }
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
      
}

void getFunctionCD(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint){
  double min = g_mygraph->GetMinimum();
  double max = g_mygraph->GetMaximum();
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-0.5,0.5);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,0.02,0.5);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-M_PI,M_PI);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-180,180);
      
}


void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfile.root \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 2)
    {
      Usage();
      return -1;
    }



  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;

  TFile * inFile = new TFile(argv[2]);

  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_pion_ptr(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector stationary_proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;


  TH2D * h_phi_corr_binSector_binTheta[6][12];
  TGraphErrors * g_phi_corr_binSector_binTheta[6][12];
  TF1 * f_phi_corr_binSector_binTheta[6][12];
  TFitResultPtr p_phi_corr_binSector_binTheta[6][12];
  TF1 * f_phi_corr_Combined_binSector_binTheta[6][12];
  for(int j=1; j<=6; j++){
    for(int i=0; i<12; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction;Counts",j,min,max);
      h_phi_corr_binSector_binTheta[j-1][i] = (TH2D*)inFile->Get(temp_name);
      hist_list.push_back(h_phi_corr_binSector_binTheta[j-1][i]);

      g_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_phi_corr_sector_%d_theta_%d",j,i);
      g_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

      sprintf(temp_name,"f_phi_corr_sector_%d_theta_%d",j,i);
      f_phi_corr_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);
      
      f_phi_corr_Combined_binSector_binTheta[j-1][i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceFD(x[0],p[0],p[1],p[2],p[3]); },-40,40,4);

    }
  }

  TH2D * h_phi_corr_binThetaCD[6];
  TGraphErrors * g_phi_corr_binThetaCD[6];
  TF1 * f_phi_corr_binThetaCD[6];
  TFitResultPtr p_phi_corr_binThetaCD[6];
  TF1 * f_phi_corr_Combined_binThetaCD[6];
  for(int i=0; i<6; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"phi_corr_theta_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< #theta < %d);#phi;Correction;Counts",min,max);
    h_phi_corr_binThetaCD[i] = (TH2D*)inFile->Get(temp_name);
    hist_list.push_back(h_phi_corr_binThetaCD[i]);
    
    g_phi_corr_binThetaCD[i] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_theta_%d",i);
    g_phi_corr_binThetaCD[i]->SetName(temp_name);
    
    sprintf(temp_name,"f_phi_corr_theta_%d",i);
    f_phi_corr_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,3);    

    sprintf(temp_name,"f_phi_corr_Combined_theta_%d",i);
    f_phi_corr_Combined_binThetaCD[i] = new TF1(temp_name,[&](double *x, double *p){ return FuncPhiDependenceCD(x[0],p[0],p[1],p[2]); },-180,180,3);    
  }

  
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
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

  ///////////////////////////////////
  TFitResultPtr p_params[6][4];
  char* names[4] = {"C_{a}","C_{b}","C_{c}","C_{d}"};
  for(int j = 0; j < 6; j++){
    TGraph * params[4];
    TF1 * f_params[4];
    double guess[4];
    for(int k = 0; k < 4; k++){
      params[k] = new TGraph();
      sprintf(temp_name,"params_%d",k);
      f_params[k] = new TF1(temp_name,[&](double *x, double *p){ return ErfFunc(x[0],p[0],p[1],p[2],p[3]);},5,50,4);
      f_params[k]->SetParameter(2,20);
      f_params[k]->SetParLimits(2,6,100);
      f_params[k]->SetParameter(3,25);
      f_params[k]->SetParLimits(3,1,50);
      if((k==0)||(k==1)){
	f_params[k]->SetParameter(0,0);
	f_params[k]->SetParLimits(0,-1,1);
	f_params[k]->SetParameter(1,0);
	f_params[k]->SetParLimits(1,-1,1);
      }
      else if(k==2){
	f_params[k]->SetParameter(0,0);
	f_params[k]->SetParLimits(0,-1000,1000);
	f_params[k]->SetParameter(1,100);
	f_params[k]->SetParLimits(1,10,1000);
      }
      else if(k==3){
	f_params[k]->SetParameter(0,0);
	f_params[k]->SetParLimits(0,-30,30);
	f_params[k]->SetParameter(1,0);
	f_params[k]->SetParLimits(1,-10,10);
      }

      if((k==0)&&(j==0)){
	f_params[k]->SetParameter(2,20);
	f_params[k]->SetParLimits(2,10,100);
      }
      if((k==2)&&(j==4)){
	f_params[k]->SetParameter(2,5);
	f_params[k]->SetParLimits(2,3,6);
      }
      if((k==2)&&(j==3)){
	f_params[k]->SetParameter(2,5);
	f_params[k]->SetParLimits(2,3,8);
      }
      
    }

    for(int i = 0; i < 12; i++){
      getGraph(h_phi_corr_binSector_binTheta[j][i],g_phi_corr_binSector_binTheta[j][i]);

      ///////////////
      if(i==0){
	guess[0] = 0;
	guess[1] = 0.25;
	guess[2] = 15;
	guess[3] = 0;
      }
      getFunction(g_phi_corr_binSector_binTheta[j][i],f_phi_corr_binSector_binTheta[j][i],p_phi_corr_binSector_binTheta[j][i],guess,i,j+1);
      ///////////////
      //f_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      for(int k = 0; k < 4; k++){
	params[k]->SetPoint(params[k]->GetN(),x,p_phi_corr_binSector_binTheta[j][i]->Parameter(k));
	guess[k]=p_phi_corr_binSector_binTheta[j][i]->Parameter(k);
      }
    }
    for(int k = 0; k < 4; k++){
      p_params[j][k] = params[k]->Fit(f_params[k],"SrBeqn","",5,45);
    }
    for(int i = 0; i < 12; i++){
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      f_phi_corr_Combined_binSector_binTheta[j][i]->SetLineColor(4);
      for(int k = 0; k < 4; k++){      
	f_phi_corr_Combined_binSector_binTheta[j][i]->SetParameter(k,f_params[k]->Eval(x));
      }
    }
    //b/b/b/b
    myCanvas->Divide(4,3);
    for(int i = 0; i < 12; i++){
      myCanvas->cd(i+1);
      h_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_phi_corr_binSector_binTheta[j][i]->SetLineColor(2);
      g_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      //f_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    myCanvas->Divide(2,2);
    for(int k = 0; k < 4; k++){
      myCanvas->cd(k+1);      
      params[k]->SetLineColor(3); 
      sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names[k],names[k]);
      params[k]->SetTitle(temp_title);
      params[k]->Draw();
      f_params[k]->SetLineColor(4);      
      f_params[k]->Draw("SAME");      
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    myCanvas->Divide(4,3);
    for(int i = 0; i < 12; i++){
      myCanvas->cd(i+1);
      h_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      g_phi_corr_binSector_binTheta[j][i]->SetLineColor(2);
      g_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
      f_phi_corr_Combined_binSector_binTheta[j][i]->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  ///////////////////////////////////
  TGraph * params_CD[3];
  TF1 * f_params_CD[3];
  TFitResultPtr p_params_CD[3];
  for(int k = 0; k < 3; k++){
    params_CD[k] = new TGraph();
    sprintf(temp_name,"params_CD_%d",k);
    if(k<2){
      f_params_CD[k] = new TF1(temp_name,[&](double *x, double *p){ return FuncThetaDependenceCD(x[0],p[0],p[1]);},30,70,2);
    }
    else{
      f_params_CD[k] = new TF1(temp_name,[&](double *x, double *p){ return FuncThetaDependenceCD2(x[0],p[0]);},30,70,1);
    }
    f_params_CD[k]->SetParameter(0,0);
    f_params_CD[k]->SetParLimits(0,-10,10);
    if(k<2){
      f_params_CD[k]->SetParameter(1,0);
      f_params_CD[k]->SetParLimits(1,-10,10);
    }
  }
  for(int i = 0; i < 6; i++){
    //myCanvas->cd(i+1);
    //h_phi_corr_binThetaCD[i]->Draw("colz");
    getGraph(h_phi_corr_binThetaCD[i],g_phi_corr_binThetaCD[i]);
    //g_phi_corr_binThetaCD[i]->Draw("SAME");
    getFunctionCD(g_phi_corr_binThetaCD[i],f_phi_corr_binThetaCD[i],p_phi_corr_binThetaCD[i]);
    //f_phi_corr_binThetaCD[i]->Draw("SAME");
    
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    for(int k = 0; k < 3; k++){
      params_CD[k]->SetPoint(params_CD[k]->GetN(),x,p_phi_corr_binThetaCD[i]->Parameter(k));
    }
  }
  
  for(int k = 0; k < 3; k++){
    //myCanvas->cd(k+7);      
    p_params_CD[k] = params_CD[k]->Fit(f_params_CD[k],"SrBeqn","",35,60);
    //f_params_CD[k]->Draw();      
    //params_CD[k]->Draw("SAME");
    for(int l = 0; l < 5; l++){
      //cout<<p_params_CD[k]->Parameter(l)<<endl;
    }
  }
  
  
  for(int i = 0; i < 6; i++){
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    f_phi_corr_Combined_binThetaCD[i]->SetLineColor(4);
    for(int k = 0; k < 3; k++){      
      f_phi_corr_Combined_binThetaCD[i]->SetParameter(k,f_params_CD[k]->Eval(x));
    }
  }
  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->SetLineColor(2);
    g_phi_corr_binThetaCD[i]->SetLineWidth(1);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  char* names_CD[3] = {"C_{a}","C_{b}","C_{c}"};
  myCanvas->Divide(2,2);
  for(int k = 0; k < 3; k++){
    myCanvas->cd(k+1);      
    sprintf(temp_title,"%s vs. #theta #circ;#theta #circ;%s",names_CD[k],names_CD[k]);
    f_params_CD[k]->SetLineColor(4);
    f_params_CD[k]->SetTitle(temp_title);      
    f_params_CD[k]->Draw();      
    params_CD[k]->SetLineColor(3);
    params_CD[k]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    g_phi_corr_binThetaCD[i]->Draw("SAME");
    f_phi_corr_Combined_binThetaCD[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  cout<<"params_FD="<<endl<<"{";
  for(int j = 0; j < 6; j++){
    cout<<"{";
    for(int k = 0; k < 4; k++){
      cout<<"{";
      for(int l = 0; l < 4; l++){
	if(l!=0){cout<<",";}
	cout<<p_params[j][k]->Parameter(l);	
      }
      cout<<"}";
      if(k!=3){cout<<","<<endl;}
    }
    cout<<"}";
    if(j!=5){cout<<","<<endl;}
  }
  cout<<"}"<<endl;

  cout<<"params_CD="<<endl<<"{";
  for(int k = 0; k < 3; k++){
    cout<<"{";
    for(int l = 0; l < 2; l++){
      if(l!=0){cout<<",";}
      if((k==2)&&(l==1)){cout<<"0";}
      else{cout<<p_params_CD[k]->Parameter(l);}      	
    }
    cout<<"}";
    if(k!=2){cout<<","<<endl;}
  }
  cout<<"}"<<endl;  
   
  return 0;
}


