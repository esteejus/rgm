//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 17 16:25:15 2021 by ROOT version 6.18/04
// from TTree gst/GENIE Summary Event Tree
// found on file: /w/hallb-scifs17exp/clas/claseg2/apapadop/apapadop_2_261GeV_SuSav2_hN.root
//////////////////////////////////////////////////////////

#ifndef gst_h
#define gst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           iev;
   Int_t           neu;
   Int_t           fspl;
   Int_t           tgt;
   Int_t           Z;
   Int_t           A;
   Int_t           hitnuc;
   Int_t           hitqrk;
   Int_t           resid;
   Bool_t          sea;
   Bool_t          qel;
   Bool_t          mec;
   Bool_t          res;
   Bool_t          dis;
   Bool_t          coh;
   Bool_t          dfr;
   Bool_t          imd;
   Bool_t          imdanh;
   Bool_t          singlek;
   Bool_t          nuel;
   Bool_t          em;
   Bool_t          cc;
   Bool_t          nc;
   Bool_t          charm;
   Bool_t          amnugamma;
   Int_t           neut_code;
   Int_t           nuance_code;
   Double_t        wght;
   Double_t        xs;
   Double_t        ys;
   Double_t        ts;
   Double_t        Q2s;
   Double_t        Ws;
   Double_t        x;
   Double_t        y;
   Double_t        t;
   Double_t        Q2;
   Double_t        W;
   Double_t        EvRF;
   Double_t        Ev;
   Double_t        pxv;
   Double_t        pyv;
   Double_t        pzv;
   Double_t        En;
   Double_t        pxn;
   Double_t        pyn;
   Double_t        pzn;
   Double_t        El;
   Double_t        pxl;
   Double_t        pyl;
   Double_t        pzl;
   Double_t        pl;
   Double_t        cthl;
   Int_t           nfp;
   Int_t           nfn;
   Int_t           nfpip;
   Int_t           nfpim;
   Int_t           nfpi0;
   Int_t           nfkp;
   Int_t           nfkm;
   Int_t           nfk0;
   Int_t           nfem;
   Int_t           nfother;
   Int_t           nip;
   Int_t           nin;
   Int_t           nipip;
   Int_t           nipim;
   Int_t           nipi0;
   Int_t           nikp;
   Int_t           nikm;
   Int_t           nik0;
   Int_t           niem;
   Int_t           niother;
   Int_t           ni;
   Int_t           pdgi[11];   //[ni]
   Int_t           resc[11];   //[ni]
   Double_t        Ei[11];   //[ni]
   Double_t        pxi[11];   //[ni]
   Double_t        pyi[11];   //[ni]
   Double_t        pzi[11];   //[ni]
   Int_t           nf;
   Int_t           pdgf[11];   //[nf]
   Double_t        Ef[11];   //[nf]
   Double_t        pxf[11];   //[nf]
   Double_t        pyf[11];   //[nf]
   Double_t        pzf[11];   //[nf]
   Double_t        pf[11];   //[nf]
   Double_t        cthf[11];   //[nf]
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Double_t        vtxt;
   Double_t        sumKEf;
   Double_t        calresp0;
   Double_t        XSec;
   Double_t        DXSec;
   UInt_t          KPS;

   // List of branches
   TBranch        *b_iev;   //!
   TBranch        *b_neu;   //!
   TBranch        *b_fspl;   //!
   TBranch        *b_tgt;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_A;   //!
   TBranch        *b_hitnuc;   //!
   TBranch        *b_hitqrk;   //!
   TBranch        *b_resid;   //!
   TBranch        *b_sea;   //!
   TBranch        *b_qel;   //!
   TBranch        *b_mec;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_dfr;   //!
   TBranch        *b_imd;   //!
   TBranch        *b_imdanh;   //!
   TBranch        *b_singlek;   //!
   TBranch        *b_nuel;   //!
   TBranch        *b_em;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_charm;   //!
   TBranch        *b_amnugamma;   //!
   TBranch        *b_neut_code;   //!
   TBranch        *b_nuance_code;   //!
   TBranch        *b_wght;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_Q2s;   //!
   TBranch        *b_Ws;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_EvRF;   //!
   TBranch        *b_Ev;   //!
   TBranch        *b_pxv;   //!
   TBranch        *b_pyv;   //!
   TBranch        *b_pzv;   //!
   TBranch        *b_En;   //!
   TBranch        *b_pxn;   //!
   TBranch        *b_pyn;   //!
   TBranch        *b_pzn;   //!
   TBranch        *b_El;   //!
   TBranch        *b_pxl;   //!
   TBranch        *b_pyl;   //!
   TBranch        *b_pzl;   //!
   TBranch        *b_pl;   //!
   TBranch        *b_cthl;   //!
   TBranch        *b_nfp;   //!
   TBranch        *b_nfn;   //!
   TBranch        *b_nfpip;   //!
   TBranch        *b_nfpim;   //!
   TBranch        *b_nfpi0;   //!
   TBranch        *b_nfkp;   //!
   TBranch        *b_nfkm;   //!
   TBranch        *b_nfk0;   //!
   TBranch        *b_nfem;   //!
   TBranch        *b_nfother;   //!
   TBranch        *b_nip;   //!
   TBranch        *b_nin;   //!
   TBranch        *b_nipip;   //!
   TBranch        *b_nipim;   //!
   TBranch        *b_nipi0;   //!
   TBranch        *b_nikp;   //!
   TBranch        *b_nikm;   //!
   TBranch        *b_nik0;   //!
   TBranch        *b_niem;   //!
   TBranch        *b_niother;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_pdgi;   //!
   TBranch        *b_resc;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_pxi;   //!
   TBranch        *b_pyi;   //!
   TBranch        *b_pzi;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_pdgf;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_pxf;   //!
   TBranch        *b_pyf;   //!
   TBranch        *b_pzf;   //!
   TBranch        *b_pf;   //!
   TBranch        *b_cthf;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_vtxt;   //!
   TBranch        *b_sumKEf;   //!
   TBranch        *b_calresp0;   //!
   TBranch        *b_XSec;   //!
   TBranch        *b_DXSec;   //!
   TBranch        *b_KPS;   //!

   gst(TTree *tree=0);
   virtual ~gst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gst_cxx
gst::gst(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/w/hallb-scifs17exp/clas/claseg2/apapadop/apapadop_2_261GeV_SuSav2_hN.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/w/hallb-scifs17exp/clas/claseg2/apapadop/apapadop_2_261GeV_SuSav2_hN.root");
      }
      f->GetObject("gst",tree);

   }
   Init(tree);
}

gst::~gst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gst::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("neu", &neu, &b_neu);
   fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
   fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
   fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
   fChain->SetBranchAddress("resid", &resid, &b_resid);
   fChain->SetBranchAddress("sea", &sea, &b_sea);
   fChain->SetBranchAddress("qel", &qel, &b_qel);
   fChain->SetBranchAddress("mec", &mec, &b_mec);
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
   fChain->SetBranchAddress("imd", &imd, &b_imd);
   fChain->SetBranchAddress("imdanh", &imdanh, &b_imdanh);
   fChain->SetBranchAddress("singlek", &singlek, &b_singlek);
   fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
   fChain->SetBranchAddress("em", &em, &b_em);
   fChain->SetBranchAddress("cc", &cc, &b_cc);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("charm", &charm, &b_charm);
   fChain->SetBranchAddress("amnugamma", &amnugamma, &b_amnugamma);
   fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
   fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
   fChain->SetBranchAddress("wght", &wght, &b_wght);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
   fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("EvRF", &EvRF, &b_EvRF);
   fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
   fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
   fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
   fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
   fChain->SetBranchAddress("En", &En, &b_En);
   fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
   fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
   fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
   fChain->SetBranchAddress("El", &El, &b_El);
   fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
   fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
   fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
   fChain->SetBranchAddress("pl", &pl, &b_pl);
   fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
   fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
   fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
   fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
   fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
   fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
   fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
   fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
   fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
   fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
   fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
   fChain->SetBranchAddress("nip", &nip, &b_nip);
   fChain->SetBranchAddress("nin", &nin, &b_nin);
   fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
   fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
   fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
   fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
   fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
   fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
   fChain->SetBranchAddress("niem", &niem, &b_niem);
   fChain->SetBranchAddress("niother", &niother, &b_niother);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
   fChain->SetBranchAddress("resc", resc, &b_resc);
   fChain->SetBranchAddress("Ei", Ei, &b_Ei);
   fChain->SetBranchAddress("pxi", pxi, &b_pxi);
   fChain->SetBranchAddress("pyi", pyi, &b_pyi);
   fChain->SetBranchAddress("pzi", pzi, &b_pzi);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
   fChain->SetBranchAddress("Ef", Ef, &b_Ef);
   fChain->SetBranchAddress("pxf", pxf, &b_pxf);
   fChain->SetBranchAddress("pyf", pyf, &b_pyf);
   fChain->SetBranchAddress("pzf", pzf, &b_pzf);
   fChain->SetBranchAddress("pf", pf, &b_pf);
   fChain->SetBranchAddress("cthf", cthf, &b_cthf);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
   fChain->SetBranchAddress("sumKEf", &sumKEf, &b_sumKEf);
   fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
   fChain->SetBranchAddress("XSec", &XSec, &b_XSec);
   fChain->SetBranchAddress("DXSec", &DXSec, &b_DXSec);
   fChain->SetBranchAddress("KPS", &KPS, &b_KPS);
   Notify();
}

Bool_t gst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gst_cxx
