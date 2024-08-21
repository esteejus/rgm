#include "TRandom3.h"

#ifndef __GCF_SRC_H__
#define __GCF_SRC_H__

const double GeVfm  =0.1973;
const int pCode = 2212;
const int nCode = 2112;
const double alpha = 0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;
const double nbGeVSq = cmSqGeVSq*1.E33;
const double Vud = 0.97417;
const double mu_p=2.79;
const double mu_n=-1.91;
const double mW = 80.379;
const double mN = 0.93892;

enum NNModel {AV18, AV4Pc, N2LO_10, N2LO_12, N3LO_600, NV2_1a, AV18_deut};

class gcfSRC
{
 public:
  gcfSRC(int thisZ, int thisN, char* uType);
  gcfSRC(int thisZ, int thisN, NNModel uType);
  ~gcfSRC();
  double get_S(double krel, int l_type, int r_type);
  double get_pp(double k_rel);
  double get_nn(double k_rel);
  double get_pn(double k_rel);
  double get_pn0(double k_rel);
  double get_pn1(double k_rel);
  NNModel get_InteractionType();
  int get_Z();
  int get_N();
  double get_Cnn0();
  double get_Cpp0();
  double get_Cpn0();
  double get_Cpn1();
  double get_d_Cnn0();
  double get_d_Cpp0();
  double get_d_Cpn0();
  double get_d_Cpn1();
  void randomize_Contacts(TRandom3* myRand);
  
  void set_Interaction(NNModel thisNNType);
  void set_Interaction(char* thisNNType);
  void set_Cpp0(double newCpp0);
  void set_Cnn0(double newCnn0);
  void set_Cpn0(double newCpn0);
  void set_Cpn1(double newCpn1);

 private:
  int Z;
  int N;
  int A;
  NNModel u;
  double phiSq_pp0[7][100];
  double phiSq_nn0[7][100];
  double phiSq_pn0[7][100];
  double phiSq_pn1[7][100];
  double Cpp0;
  double d_Cpp0;
  double Cnn0;
  double d_Cnn0;
  double Cpn0;
  double d_Cpn0;
  double Cpn1;
  double d_Cpn1;
  
  double get_phiSq(double *phiPtr, double k_rel);

  void set_Contacts();
  bool set_Contacts_SS_r();
  bool set_Contacts_SS_k();
  bool set_Contacts_deut();
  bool set_Contacts_EG2();

  void fill_arrays();
  void fill_arrays_AV18();
  void fill_arrays_n2lo_local();
  void fill_arrays_n3lo_nonlocal();
  void fill_arrays_n2lo_12_local();
  void fill_arrays_AV4Pc();
  void fill_arrays_NV2_1a();
  void fill_arrays_AV18_deut();
    
};

#endif
