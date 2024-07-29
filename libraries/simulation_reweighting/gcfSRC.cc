#include <iostream>
#include "gcfSRC.hh"
#include "universal_functions/AV18.hh"
#include "universal_functions/N2LO.hh"
#include "universal_functions/N3LO.hh"
#include "universal_functions/N2LO_12.hh"
#include "universal_functions/AV4Pc.hh"
#include "universal_functions/NV2_1a.hh"
#include "universal_functions/AV18_deut.hh"

gcfSRC::gcfSRC(int thisZ, int thisN, char* uType)
{
  fill_arrays();
  set_Interaction(uType);
  Z = thisZ;
  N = thisN;
  A=Z+N;
  set_Contacts();
}  

gcfSRC::gcfSRC(int thisZ, int thisN, NNModel uType)
{
  fill_arrays();
  set_Interaction(uType);
  Z = thisZ;
  N = thisN;
  A=Z+N;
  set_Contacts();
}

gcfSRC::~gcfSRC()
{
}


void gcfSRC::randomize_Contacts(TRandom3* myRand)
{
  Cpp0 += myRand->Gaus(0.,d_Cpp0);
  Cpn0 += myRand->Gaus(0.,d_Cpn0);
  Cnn0 += myRand->Gaus(0.,d_Cnn0);
  Cpn1 += myRand->Gaus(0.,d_Cpn1);
}

void gcfSRC::set_Interaction(char* thisPType){

  if(strcmp(thisPType,"AV18")==0 or strcmp(thisPType,"1")==0)
    set_Interaction(AV18);
  else if (strcmp(thisPType,"N2LO")==0 or strcmp(thisPType,"N2LO10")==0 or strcmp(thisPType,"N2LO_10")==0 or strcmp(thisPType,"2")==0)
    set_Interaction(N2LO_10);
  else if (strcmp(thisPType,"N2LO12")==0 or strcmp(thisPType,"N2LO_12")==0 or strcmp(thisPType,"4")==0)
    set_Interaction(N2LO_12);
  else if (strcmp(thisPType,"N3LO")==0 or strcmp(thisPType,"N3LO600")==0 or strcmp(thisPType,"N3LO_600")==0 or strcmp(thisPType,"3")==0)
    set_Interaction(N3LO_600);
  else if (strcmp(thisPType,"AV4Pc")==0 or strcmp(thisPType,"AV4")==0 or strcmp(thisPType,"5")==0)
    set_Interaction(AV4Pc);
  else if (strcmp(thisPType,"NV")==0 or strcmp(thisPType,"NV2_1a")==0 or strcmp(thisPType,"6")==0)
    set_Interaction(NV2_1a);
  else if (strcmp(thisPType,"AV18_deut")==0 or strcmp(thisPType,"7")==0)
    set_Interaction(AV18_deut);
  else{
    std::cerr <<"You are using an interaction not in the library. \n Aborting...\n";
  exit(-2);
  }
  
}

void gcfSRC::set_Interaction(NNModel thisPType){
  u = thisPType;
  if (u == AV18){
    std::cout <<"You are using the AV18 interaction.\n";
  }
  else if (u == N2LO_10){
    std::cout <<"You are using the N2L0 interaction calculated with 1.0 fm cutoff.\n";
  }
  else if (u == N3LO_600){
    std::cout <<"You are using the N3L0 interaction\n";
  }
  else if (u == N2LO_12){
    std::cout <<"You are using the N2L0 interaction calculated with 1.2 fm cutoff.\n";
  }
  else if(u == AV4Pc){
    std::cout <<"You are using the AV4' interaction.\n";
  }
  else if(u == NV2_1a){
    std::cout <<"You are using the NV2+Ia interaction.\n";
  }
  else if (u == AV18_deut){
    std::cout <<"You are using the AV18 interaction valid at all momentum ranges.\n";
  }
  else{
    std::cerr <<"You are using an interaction not in the library. \n Aborting...\n";
  exit(-2);
  }
  
}

void gcfSRC::set_Cpp0(double newCpp0){

  Cpp0 = newCpp0;
  
}

void gcfSRC::set_Cnn0(double newCnn0){

  Cnn0 = newCnn0;
  
}

void gcfSRC::set_Cpn0(double newCpn0){

  Cpn0 = newCpn0;
  
}

void gcfSRC::set_Cpn1(double newCpn1){

  Cpn1 = newCpn1;
  
}

NNModel gcfSRC::get_InteractionType(){

  return u;
  
}

double gcfSRC::get_Cpp0(){

  return Cpp0;
  
}

double gcfSRC::get_Cnn0(){

  return Cnn0;
  
}

double gcfSRC::get_Cpn0(){

  return Cpn0;
  
}

double gcfSRC::get_Cpn1(){

  return Cpn1;
  
}

double gcfSRC::get_d_Cpp0(){

  return d_Cpp0;
  
}

double gcfSRC::get_d_Cnn0(){

  return d_Cnn0;
  
}

double gcfSRC::get_d_Cpn0(){

  return d_Cpn0;
  
}

double gcfSRC::get_d_Cpn1(){

  return d_Cpn1;
  
}

int gcfSRC::get_Z()
{
  
  return Z;
  
}

int gcfSRC::get_N()
{
  
  return N;
  
}


double gcfSRC::get_S(double k_rel, int l_type, int r_type){
  if(l_type==r_type){
    return ((l_type==pCode) ? (get_pp(k_rel)) : (get_nn(k_rel)));
   }
  else{
    return get_pn(k_rel);
  }
}

double gcfSRC::get_pp(double k_rel)
{
  return 2. * Cpp0 * get_phiSq(phiSq_pp0[u],k_rel); // The 2 comes from contact definition
}

double gcfSRC::get_nn(double k_rel)
{
  return 2. * Cnn0 * get_phiSq(phiSq_nn0[u],k_rel); 
}

double gcfSRC::get_pn(double k_rel)
{
  return get_pn0(k_rel) + get_pn1(k_rel);
}

double gcfSRC::get_pn0(double k_rel)
{
  return Cpn0 * get_phiSq(phiSq_pn0[u],k_rel);
}

double gcfSRC::get_pn1(double k_rel)
{
  return Cpn1 * get_phiSq(phiSq_pn1[u],k_rel);
}

double gcfSRC::get_phiSq(double *phiPtr, double k_rel)
{

  double bin = k_rel / GeVfm / 0.1;

  if (bin < 0.)
    return 0.;
  if (bin < 1.)
    return bin * phiPtr[0];
  if (bin > 100.)
    return 0.;
  
  int b = bin;
  double x = bin - b;
  return (x*phiPtr[b] + (1.-x)*phiPtr[b-1]) / pow(GeVfm,3) / 100.;
  
}

void gcfSRC::set_Contacts()
{
    if (set_Contacts_SS_k())
    {
      std::cout << "You are using k-space contact values from the Scale and Scheme paper.\n";
    }
  else if (set_Contacts_SS_r())
    {
      std::cout << "You are using r-space contact values from the Scale and Scheme paper.\n";
    }
  else if (set_Contacts_deut())
    {
      std::cout << "You are using a deuteron momentum distribution, and therefore have no contact dependence.\n";
    }
  else if (set_Contacts_EG2())
    {
      std::cout << "You are using contact ratios from fits to EG2 data. You must be truly desperate...\n";
    }
  else
    {
      std::cerr << "You selected a nucleus with Z=" << Z << " and with N=" << N << ".\n"
		<< "This combination of interaction and nucleus does not have contacts in the library. Aborting...\n";
      exit(-2);
    }
    
}

bool gcfSRC::set_Contacts_SS_r()
{
  if ((Z==1) && (N==1))
    {
      if (u==AV18)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.898 * ( A * 0.5);
	  d_Cpn1=0.080 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=1.186 * ( A * 0.5);
	  d_Cpn1=0.034 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.664 * ( A * 0.5);
	  d_Cpn1=0.009 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.141 * ( A * 0.5);
	  d_Cpn1=0.010 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=3.878 * ( A * 0.5);
	  d_Cpn1=0.390 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==1) && (N==2))
    {
      if (u==AV18)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.549 * ( A * 0.5);
	  d_Cnn0=0.055 * ( A * 0.5);
	  Cpn0=0.295 * ( A * 0.5);
	  d_Cpn0=0.119 * ( A * 0.5);
	  Cpn1=6.246 * ( A * 0.5);
	  d_Cpn1=0.856 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.544 * ( A * 0.5);
	  d_Cnn0=0.055 * ( A * 0.5);
	  Cpn0=0.270 * ( A * 0.5);
	  d_Cpn0=0.027 * ( A * 0.5);
	  Cpn1=1.533 * ( A * 0.5);
	  d_Cpn1=0.154 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.282 * ( A * 0.5);
	  d_Cnn0=0.032 * ( A * 0.5);
	  Cpn0=0.140 * ( A * 0.5);
	  d_Cpn0=0.015 * ( A * 0.5);
	  Cpn1=6.237 * ( A * 0.5);
	  d_Cpn1=0.718 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.277 * ( A * 0.5);
	  d_Cnn0=0.033 * ( A * 0.5);
	  Cpn0=0.139 * ( A * 0.5);
	  d_Cpn0=0.017 * ( A * 0.5);
	  Cpn1=5.980 * ( A * 0.5);
	  d_Cpn1=0.695 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==2) && (N==1))
    {
      if (u==AV18)
	{
	  Cpp0=0.536 * ( A * 0.5);
	  d_Cpp0=0.054 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.161 * ( A * 0.5);
	  d_Cpn0=0.092 * ( A * 0.5);
	  Cpn1=6.851 * ( A * 0.5);
	  d_Cpn1=0.822 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.515 * ( A * 0.5);
	  d_Cpp0=0.052 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.265 * ( A * 0.5);
	  d_Cpn0=0.027 * ( A * 0.5);
	  Cpn1=1.527 * ( A * 0.5);
	  d_Cpn1=0.153 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.259 * ( A * 0.5);
	  d_Cpp0=0.032 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.135 * ( A * 0.5);
	  d_Cpn0=0.016 * ( A * 0.5);
	  Cpn1=6.113 * ( A * 0.5);
	  d_Cpn1=0.719 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.239 * ( A * 0.5);
	  d_Cpp0=0.041 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.127 * ( A * 0.5);
	  d_Cpn0=0.022 * ( A * 0.5);
	  Cpn1=5.705 * ( A * 0.5);
	  d_Cpn1=0.870 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.310 * ( A * 0.5);
	  d_Cpp0=0.031 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.234 * ( A * 0.5);
	  d_Cpn0=0.053 * ( A * 0.5);
	  Cpn1=1.637 * ( A * 0.5);
	  d_Cpn1=0.573 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==2) && (N==2))
    {
      if (u==AV18)
	{
	  Cpp0=0.567 * ( A * 0.5);
	  d_Cpp0=0.057 * ( A * 0.5);
	  Cnn0=0.567 * ( A * 0.5);
	  d_Cnn0=0.057 * ( A * 0.5);
	  Cpn0=0.567 * ( A * 0.5);
	  d_Cpn0=0.057 * ( A * 0.5);
	  Cpn1=11.605 * ( A * 0.5);
	  d_Cpn1=1.161 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.564 * ( A * 0.5);
	  d_Cpp0=0.057 * ( A * 0.5);
	  Cnn0=0.564 * ( A * 0.5);
	  d_Cnn0=0.057 * ( A * 0.5);
	  Cpn0=0.578 * ( A * 0.5);
	  d_Cpn0=0.059 * ( A * 0.5);
	  Cpn1=2.685 * ( A * 0.5);
	  d_Cpn1=0.272 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.243 * ( A * 0.5);
	  d_Cpp0=0.040 * ( A * 0.5);
	  Cnn0=0.243 * ( A * 0.5);
	  d_Cnn0=0.040 * ( A * 0.5);
	  Cpn0=0.253 * ( A * 0.5);
	  d_Cpn0=0.043 * ( A * 0.5);
	  Cpn1=10.508 * ( A * 0.5);
	  d_Cpn1=1.308 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.263 * ( A * 0.5);
	  d_Cpp0=0.059 * ( A * 0.5);
	  Cnn0=0.263 * ( A * 0.5);
	  d_Cnn0=0.059 * ( A * 0.5);
	  Cpn0=0.281 * ( A * 0.5);
	  d_Cpn0=0.062 * ( A * 0.5);
	  Cpn1=11.111 * ( A * 0.5);
	  d_Cpn1=2.595 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.333 * ( A * 0.5);
	  d_Cpp0=0.034 * ( A * 0.5);
	  Cnn0=0.333 * ( A * 0.5);
	  d_Cnn0=0.034 * ( A * 0.5);
	  Cpn0=0.333 * ( A * 0.5);
	  d_Cpn0=0.034 * ( A * 0.5);
	  Cpn1=9.200 * ( A * 0.5);
	  d_Cpn1=0.928 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==3) && (N==3))
    {
      if (u==AV18)
	{
	  Cpp0=0.415 * ( A * 0.5);
	  d_Cpp0=0.042 * ( A * 0.5);
	  Cnn0=0.415 * ( A * 0.5);
	  d_Cnn0=0.042 * ( A * 0.5);
	  Cpn0=0.415 * ( A * 0.5);
	  d_Cpn0=0.042 * ( A * 0.5);
	  Cpn1=10.140 * ( A * 0.5);
	  d_Cpn1=1.015 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.380 * ( A * 0.5);
	  d_Cpp0=0.038 * ( A * 0.5);
	  Cnn0=0.380 * ( A * 0.5);
	  d_Cnn0=0.038 * ( A * 0.5);
	  Cpn0=0.387 * ( A * 0.5);
	  d_Cpn0=0.039 * ( A * 0.5);
	  Cpn1=2.248 * ( A * 0.5);
	  d_Cpn1=0.225 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.173 * ( A * 0.5);
	  d_Cpp0=0.020 * ( A * 0.5);
	  Cnn0=0.173 * ( A * 0.5);
	  d_Cnn0=0.020 * ( A * 0.5);
	  Cpn0=0.180 * ( A * 0.5);
	  d_Cpn0=0.021 * ( A * 0.5);
	  Cpn1=8.434 * ( A * 0.5);
	  d_Cpn1=1.026 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.185 * ( A * 0.5);
	  d_Cpp0=0.031 * ( A * 0.5);
	  Cnn0=0.185 * ( A * 0.5);
	  d_Cnn0=0.031 * ( A * 0.5);
	  Cpn0=0.197 * ( A * 0.5);
	  d_Cpn0=0.034 * ( A * 0.5);
	  Cpn1=9.011 * ( A * 0.5);
	  d_Cpn1=1.478 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==6) && (N==6))
    {
      if (u==AV18)
	{
	  Cpp0=0.716 * ( A * 0.5);
	  d_Cpp0=0.075 * ( A * 0.5);
	  Cnn0=0.716 * ( A * 0.5);
	  d_Cnn0=0.075 * ( A * 0.5);
	  Cpn0=0.716 * ( A * 0.5);
	  d_Cpn0=0.075 * ( A * 0.5);
	  Cpn1=13.135 * ( A * 0.5);
	  d_Cpn1=1.324 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.547 * ( A * 0.5);
	  d_Cpp0=0.055 * ( A * 0.5);
	  Cnn0=0.547 * ( A * 0.5);
	  d_Cnn0=0.055 * ( A * 0.5);
	  Cpn0=0.559 * ( A * 0.5);
	  d_Cpn0=0.056 * ( A * 0.5);
	  Cpn1=2.458 * ( A * 0.5);
	  d_Cpn1=0.249 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.308 * ( A * 0.5);
	  d_Cpp0=0.033 * ( A * 0.5);
	  Cnn0=0.308 * ( A * 0.5);
	  d_Cnn0=0.033 * ( A * 0.5);
	  Cpn0=0.318 * ( A * 0.5);
	  d_Cpn0=0.034 * ( A * 0.5);
	  Cpn1=10.434 * ( A * 0.5);
	  d_Cpn1=1.044 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==8) && (N==8))
    {
      if (u==AV18)
	{
	  Cpp0=0.676 * ( A * 0.5);
	  d_Cpp0=0.072 * ( A * 0.5);
	  Cnn0=0.676 * ( A * 0.5);
	  d_Cnn0=0.072 * ( A * 0.5);
	  Cpn0=0.676 * ( A * 0.5);
	  d_Cpn0=0.072 * ( A * 0.5);
	  Cpn1=11.372 * ( A * 0.5);
	  d_Cpn1=1.158 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.658 * ( A * 0.5);
	  d_Cpp0=0.066 * ( A * 0.5);
	  Cnn0=0.658 * ( A * 0.5);
	  d_Cnn0=0.066 * ( A * 0.5);
	  Cpn0=0.675 * ( A * 0.5);
	  d_Cpn0=0.068 * ( A * 0.5);
	  Cpn1=2.910 * ( A * 0.5);
	  d_Cpn1=0.293 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.270 * ( A * 0.5);
	  d_Cpp0=0.034 * ( A * 0.5);
	  Cnn0=0.270 * ( A * 0.5);
	  d_Cnn0=0.034 * ( A * 0.5);
	  Cpn0=0.275 * ( A * 0.5);
	  d_Cpn0=0.033 * ( A * 0.5);
	  Cpn1=9.103 * ( A * 0.5);
	  d_Cpn1=1.020 * ( A * 0.5);
	  return true;
	}
    }
  else if (((Z==20) && (N==20)) || ((Z==18) && (N==22)))
    {
      if (u==AV18)
	{
	  Cpp0=0.723 * ( A * 0.5);
	  d_Cpp0=0.081 * ( A * 0.5);
	  Cnn0=0.723 * ( A * 0.5);
	  d_Cnn0=0.081 * ( A * 0.5);
	  Cpn0=0.723 * ( A * 0.5);
	  d_Cpn0=0.081 * ( A * 0.5);
	  Cpn1=11.570 * ( A * 0.5);
	  d_Cpn1=1.196 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.834 * ( A * 0.5);
	  d_Cpp0=0.084 * ( A * 0.5);
	  Cnn0=0.834 * ( A * 0.5);
	  d_Cnn0=0.084 * ( A * 0.5);
	  Cpn0=0.854 * ( A * 0.5);
	  d_Cpn0=0.086 * ( A * 0.5);
	  Cpn1=3.284 * ( A * 0.5);
	  d_Cpn1=0.339 * ( A * 0.5);
	  return true;
	}
    }
  return false;
}

bool gcfSRC::set_Contacts_SS_k()
{
  if ((Z==1) && (N==1))
    {
      if (u==AV18)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.764 * ( A * 0.5);
	  d_Cpn1=0.007 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=1.165 * ( A * 0.5);
	  d_Cpn1=0.037 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.691 * ( A * 0.5);
	  d_Cpn1=0.030 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=4.244 * ( A * 0.5);
	  d_Cpn1=0.032 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.0 * ( A * 0.5);
	  d_Cpn0=0.0 * ( A * 0.5);
	  Cpn1=3.840 * ( A * 0.5);
	  d_Cpn1=0.398 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==1) && (N==2))
    {
      if (u==AV18)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.590 * ( A * 0.5);
	  d_Cnn0=0.060 * ( A * 0.5);
	  Cpn0=0.311 * ( A * 0.5);
	  d_Cpn0=0.033 * ( A * 0.5);
	  Cpn1=6.441 * ( A * 0.5);
	  d_Cpn1=0.645 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.532 * ( A * 0.5);
	  d_Cnn0=0.054 * ( A * 0.5);
	  Cpn0=0.281 * ( A * 0.5);
	  d_Cpn0=0.029 * ( A * 0.5);
	  Cpn1=1.515 * ( A * 0.5);
	  d_Cpn1=0.152 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.521 * ( A * 0.5);
	  d_Cnn0=0.060 * ( A * 0.5);
	  Cpn0=0.255 * ( A * 0.5);
	  d_Cpn0=0.063 * ( A * 0.5);
	  Cpn1=6.885 * ( A * 0.5);
	  d_Cpn1=0.789 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.0 * ( A * 0.5);
	  d_Cpp0=0.0 * ( A * 0.5);
	  Cnn0=0.687 * ( A * 0.5);
	  d_Cnn0=0.119 * ( A * 0.5);
	  Cpn0=0.337 * ( A * 0.5);
	  d_Cpn0=0.061 * ( A * 0.5);
	  Cpn1=6.111 * ( A * 0.5);
	  d_Cpn1=1.011 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==2) && (N==1))
    {
      if (u==AV18)
	{
	  Cpp0=0.570 * ( A * 0.5);
	  d_Cpp0=0.058 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.325 * ( A * 0.5);
	  d_Cpn0=0.034 * ( A * 0.5);
	  Cpn1=6.249 * ( A * 0.5);
	  d_Cpn1=0.625 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.515 * ( A * 0.5);
	  d_Cpp0=0.052 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.264 * ( A * 0.5);
	  d_Cpn0=0.028 * ( A * 0.5);
	  Cpn1=1.466 * ( A * 0.5);
	  d_Cpn1=0.147 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.454 * ( A * 0.5);
	  d_Cpp0=0.069 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.264 * ( A * 0.5);
	  d_Cpn0=0.047 * ( A * 0.5);
	  Cpn1=6.703 * ( A * 0.5);
	  d_Cpn1=0.702 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.547 * ( A * 0.5);
	  d_Cpp0=0.079 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.307 * ( A * 0.5);
	  d_Cpn0=0.057 * ( A * 0.5);
	  Cpn1=5.717 * ( A * 0.5);
	  d_Cpn1=0.730 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.336 * ( A * 0.5);
	  d_Cpp0=0.034 * ( A * 0.5);
	  Cnn0=0.0 * ( A * 0.5);
	  d_Cnn0=0.0 * ( A * 0.5);
	  Cpn0=0.233 * ( A * 0.5);
	  d_Cpn0=0.024 * ( A * 0.5);
	  Cpn1=4.885 * ( A * 0.5);
	  d_Cpn1=0.491 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==2) && (N==2))
    {
      if (u==AV18)
	{
	  Cpp0=0.655 * ( A * 0.5);
	  d_Cpp0=0.071 * ( A * 0.5);
	  Cnn0=0.655 * ( A * 0.5);
	  d_Cnn0=0.071 * ( A * 0.5);
	  Cpn0=0.687 * ( A * 0.5);
	  d_Cpn0=0.075 * ( A * 0.5);
	  Cpn1=12.274 * ( A * 0.5);
	  d_Cpn1=1.232 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.542 * ( A * 0.5);
	  d_Cpp0=0.055 * ( A * 0.5);
	  Cnn0=0.542 * ( A * 0.5);
	  d_Cnn0=0.055 * ( A * 0.5);
	  Cpn0=0.564 * ( A * 0.5);
	  d_Cpn0=0.057 * ( A * 0.5);
	  Cpn1=2.995 * ( A * 0.5);
	  d_Cpn1=0.300 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.655 * ( A * 0.5);
	  d_Cpp0=0.093 * ( A * 0.5);
	  Cnn0=0.655 * ( A * 0.5);
	  d_Cnn0=0.093 * ( A * 0.5);
	  Cpn0=0.703 * ( A * 0.5);
	  d_Cpn0=0.096 * ( A * 0.5);
	  Cpn1=12.372 * ( A * 0.5);
	  d_Cpn1=1.372 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.851 * ( A * 0.5);
	  d_Cpp0=0.102 * ( A * 0.5);
	  Cnn0=0.851 * ( A * 0.5);
	  d_Cnn0=0.102 * ( A * 0.5);
	  Cpn0=0.934 * ( A * 0.5);
	  d_Cpn0=0.133 * ( A * 0.5);
	  Cpn1=10.446 * ( A * 0.5);
	  d_Cpn1=2.223 * ( A * 0.5);
	  return true;
	}
      else if (u==NV2_1a)
	{
	  Cpp0=0.355 * ( A * 0.5);
	  d_Cpp0=0.039 * ( A * 0.5);
	  Cnn0=0.355 * ( A * 0.5);
	  d_Cnn0=0.039 * ( A * 0.5);
	  Cpn0=0.504 * ( A * 0.5);
	  d_Cpn0=0.059 * ( A * 0.5);
	  Cpn1=10.143 * ( A * 0.5);
	  d_Cpn1=1.022 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==3) && (N==3))
    {
      if (u==AV18)
	{
	  Cpp0=0.485 * ( A * 0.5);
	  d_Cpp0=0.058 * ( A * 0.5);
	  Cnn0=0.485 * ( A * 0.5);
	  d_Cnn0=0.058 * ( A * 0.5);
	  Cpn0=0.529 * ( A * 0.5);
	  d_Cpn0=0.071 * ( A * 0.5);
	  Cpn1=10.492 * ( A * 0.5);
	  d_Cpn1=1.056 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.380 * ( A * 0.5);
	  d_Cpp0=0.038 * ( A * 0.5);
	  Cnn0=0.380 * ( A * 0.5);
	  d_Cnn0=0.038 * ( A * 0.5);
	  Cpn0=0.369 * ( A * 0.5);
	  d_Cpn0=0.039 * ( A * 0.5);
	  Cpn1=2.205 * ( A * 0.5);
	  d_Cpn1=0.222 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.501 * ( A * 0.5);
	  d_Cpp0=0.074 * ( A * 0.5);
	  Cnn0=0.501 * ( A * 0.5);
	  d_Cnn0=0.074 * ( A * 0.5);
	  Cpn0=0.540 * ( A * 0.5);
	  d_Cpn0=0.086 * ( A * 0.5);
	  Cpn1=9.444 * ( A * 0.5);
	  d_Cpn1=1.141 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_12)
	{
	  Cpp0=0.668 * ( A * 0.5);
	  d_Cpp0=0.104 * ( A * 0.5);
	  Cnn0=0.668 * ( A * 0.5);
	  d_Cnn0=0.105 * ( A * 0.5);
	  Cpn0=0.749 * ( A * 0.5);
	  d_Cpn0=0.168 * ( A * 0.5);
	  Cpn1=8.650 * ( A * 0.5);
	  d_Cpn1=1.545 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==6) && (N==6))
    {
      if (u==AV18)
	{
	  Cpp0=1.140 * ( A * 0.5);
	  d_Cpp0=0.210 * ( A * 0.5);
	  Cnn0=1.140 * ( A * 0.5);
	  d_Cnn0=0.210 * ( A * 0.5);
	  Cpn0=1.244 * ( A * 0.5);
	  d_Cpn0=0.319 * ( A * 0.5);
	  Cpn1=15.876 * ( A * 0.5);
	  d_Cpn1=1.770 * ( A * 0.5);
	  return true;
	}
      else if (u==AV4Pc)
	{
	  Cpp0=0.653 * ( A * 0.5);
	  d_Cpp0=0.067 * ( A * 0.5);
	  Cnn0=0.653 * ( A * 0.5);
	  d_Cnn0=0.067 * ( A * 0.5);
	  Cpn0=0.558 * ( A * 0.5);
	  d_Cpn0=0.069 * ( A * 0.5);
	  Cpn1=2.676 * ( A * 0.5);
	  d_Cpn1=0.272 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.870 * ( A * 0.5);
	  d_Cpp0=0.095 * ( A * 0.5);
	  Cnn0=0.870 * ( A * 0.5);
	  d_Cnn0=0.095 * ( A * 0.5);
	  Cpn0=0.988 * ( A * 0.5);
	  d_Cpn0=0.161 * ( A * 0.5);
	  Cpn1=10.643 * ( A * 0.5);
	  d_Cpn1=1.094 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==8) && (N==8))
    {
      if (u==AV4Pc)
	{
	  Cpp0=0.784 * ( A * 0.5);
	  d_Cpp0=0.084 * ( A * 0.5);
	  Cnn0=0.784 * ( A * 0.5);
	  d_Cnn0=0.084 * ( A * 0.5);
	  Cpn0=0.702 * ( A * 0.5);
	  d_Cpn0=0.086 * ( A * 0.5);
	  Cpn1=2.911 * ( A * 0.5);
	  d_Cpn1=0.307 * ( A * 0.5);
	  return true;
	}
      else if (u==N2LO_10)
	{
	  Cpp0=0.781 * ( A * 0.5);
	  d_Cpp0=0.173 * ( A * 0.5);
	  Cnn0=0.781 * ( A * 0.5);
	  d_Cnn0=0.173 * ( A * 0.5);
	  Cpn0=0.928 * ( A * 0.5);
	  d_Cpn0=0.365 * ( A * 0.5);
	  Cpn1=10.338 * ( A * 0.5);
	  d_Cpn1=1.310 * ( A * 0.5);
	  return true;
	}
    }
  else if ((Z==20) && (N==20))
    {
      if (u==AV4Pc)
	{
	  Cpp0=1.329 * ( A * 0.5);
	  d_Cpp0=0.144 * ( A * 0.5);
	  Cnn0=1.329 * ( A * 0.5);
	  d_Cnn0=0.144 * ( A * 0.5);
	  Cpn0=1.357 * ( A * 0.5);
	  d_Cpn0=0.164 * ( A * 0.5);
	  Cpn1=4.476 * ( A * 0.5);
	  d_Cpn1=0.460 * ( A * 0.5);
	  return true;
	}
    }
  return false;
}

bool gcfSRC::set_Contacts_deut()
{
  if ((Z==1) && (N==1))
    {
      if (u==AV18_deut)
	{
	  Cpp0=0.0;
	  d_Cpp0=0.0;
	  Cnn0=0.0;
	  d_Cnn0=0.0;
	  Cpn0=0.0;
	  d_Cpn0=0.0;
	  Cpn1=100.0;
	  d_Cpn1=0.0;
	  return true;
	}
    }
  return false;
}

bool gcfSRC::set_Contacts_EG2()
{
  if ((Z==6) && (N==6))
    {
      if (u==AV18)
	{
	  double R01 = 0.067;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_10)
	{
	  double R01 = 0.081;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_12)
	{
	  double R01 = 0.151;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
    }
  else if ((Z==13) && (N==14))
    {
      if (u==AV18)
	{
	  double R01 = 0.045;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_10)
	{
	  double R01 = 0.053;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_12)
	{
	  double R01 = 0.097;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
    }
  else if ((Z==26) && (N==30))
    {
      if (u==AV18)
	{
	  double R01 = 0.064;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_10)
	{
	  double R01 = 0.080;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_12)
	{
	  double R01 = 0.132;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0;
	  d_Cpp0=0.;
	  Cnn0=C0;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
    }
    else if ((Z==82) && (N==126))
    {
      if (u==AV18)
	{
	  double R01 = 0.014;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0*Z/N;
	  d_Cpp0=0.;
	  Cnn0=C0*N/Z;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_10)
	{
	  double R01 = 0.018;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0*Z/N;
	  d_Cpp0=0.;
	  Cnn0=C0*N/Z;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
      else if (u==N2LO_12)
	{
	  double R01 = 0.036;
	  double C1 = 16.;
	  double C0 = R01*C1;
	    
	  Cpp0=C0*Z/N;
	  d_Cpp0=0.;
	  Cnn0=C0*N/Z;
	  d_Cnn0=0.;
	  Cpn0=C0;
	  d_Cpn0=0.;
	  Cpn1=C1;
	  d_Cpn1=0.;
	  return true;
	}
    }
  return false;
}

void gcfSRC::fill_arrays()
{
  fill_arrays_AV18();
  fill_arrays_n2lo_local();
  fill_arrays_n3lo_nonlocal();
  fill_arrays_n2lo_12_local();
  fill_arrays_AV4Pc();
  fill_arrays_NV2_1a();
  fill_arrays_AV18_deut();
}

void gcfSRC::fill_arrays_AV18()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[AV18][i] = AV18_pp0[i];
      phiSq_nn0[AV18][i] = AV18_pp0[i];
      phiSq_pn0[AV18][i] = AV18_pn0[i];
      phiSq_pn1[AV18][i] = AV18_pn1[i];
    }
}

void gcfSRC::fill_arrays_n2lo_local()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[N2LO_10][i] = N2LO_pp0[i];
      phiSq_nn0[N2LO_10][i] = N2LO_pp0[i];
      phiSq_pn0[N2LO_10][i] = N2LO_pn0[i];
      phiSq_pn1[N2LO_10][i] = N2LO_pn1[i];
    }
}

void gcfSRC::fill_arrays_n3lo_nonlocal()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[N3LO_600][i] = N3LO_pp0[i];
      phiSq_nn0[N3LO_600][i] = N3LO_pp0[i];
      phiSq_pn0[N3LO_600][i] = N3LO_pn0[i];
      phiSq_pn1[N3LO_600][i] = N3LO_pn1[i];
    }
}

void gcfSRC::fill_arrays_n2lo_12_local()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[N2LO_12][i] = N2LO_12_pp0[i];
      phiSq_nn0[N2LO_12][i] = N2LO_12_pp0[i];
      phiSq_pn0[N2LO_12][i] = N2LO_12_pn0[i];
      phiSq_pn1[N2LO_12][i] = N2LO_12_pn1[i];
    }
}

void gcfSRC::fill_arrays_AV4Pc()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[AV4Pc][i] = AV4Pc_pp0[i];
      phiSq_nn0[AV4Pc][i] = AV4Pc_pp0[i];
      phiSq_pn0[AV4Pc][i] = AV4Pc_pn0[i];
      phiSq_pn1[AV4Pc][i] = AV4Pc_pn1[i];
    }
}

void gcfSRC::fill_arrays_NV2_1a()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[NV2_1a][i] = NV2_1a_pp0[i];
      phiSq_nn0[NV2_1a][i] = NV2_1a_pp0[i];
      phiSq_pn0[NV2_1a][i] = NV2_1a_pn0[i];
      phiSq_pn1[NV2_1a][i] = NV2_1a_pn1[i];
    }
}

void gcfSRC::fill_arrays_AV18_deut()
{
  for (int i=0 ; i<100; i++)
    {
      phiSq_pp0[AV18_deut][i] = AV18_deut_pp0[i];
      phiSq_nn0[AV18_deut][i] = AV18_deut_pp0[i];
      phiSq_pn0[AV18_deut][i] = AV18_deut_pn0[i];
      phiSq_pn1[AV18_deut][i] = AV18_deut_pn1[i];
    }
}
