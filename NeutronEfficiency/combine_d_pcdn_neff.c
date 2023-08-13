void combine_d_pcdn_neff() {

  TFile * f = new TFile("allevents_2gev_pCDn.root","READ");
  TH2D * h_det2d = (TH2D*)f->Get("det2d");
  TH2D * h_cand2d = (TH2D*)f->Get("cand2d");
  TH2D * h_neff = (TH2D*)h_det2d->Clone();
  h_neff->Divide(h_cand2d);
  h_neff->SetMaximum(0.16);
  h_neff->Draw("colz");

  TFile * fout = new TFile("neff_2gev_pCDn.root","RECREATE");
  fout->cd();
  h_neff->Write();

}
