void print_goodn_plots() {

  TFile * f = new TFile("all_goodn.root","READ");

  //char * pdfFile = new char("plots_all_goodn_pCD.pdf");
  //TString pdfFile("plots_all_goodn_pCD.pdf");
  char const * pdfFile = "plots_all_goodn_pCD.pdf";


  // proton selection
  TH2D * h_pangles = (TH2D*)f->Get("pangles");
  TH2D * h_dbeta_p_fd = (TH2D*)f->Get("dbeta_p_fd");
  TH1D * h_vzp_fd = (TH1D*)f->Get("vzp_fd");
  TH1D * h_chipid_fd = (TH1D*)f->Get("chipid_fd");
  TH2D * h_dbeta_p_cd = (TH2D*)f->Get("dbeta_p_cd");
  TH1D * h_vzp_cd = (TH1D*)f->Get("vzp_cd");
  TH1D * h_chipid_cd = (TH1D*)f->Get("chipid_cd");


  // neutron candidate selection
  TH1D * h_tof = (TH1D*)f->Get("tof");
  TH2D * h_mmiss_pn = (TH2D*)f->Get("mmiss_pn");
  TH1D * h_nangles = (TH1D*)f->Get("nangles");


  // neutron signal selection
  TH1D * h_pxminuspx2 = (TH1D*)f->Get("pxminuspx2");
  TH1D * h_pyminuspy2 = (TH1D*)f->Get("pyminuspy2");
  TH1D * h_pzminuspz2 = (TH1D*)f->Get("pzminuspz2");
  TH1D * h_pminusp2 = (TH1D*)f->Get("pminusp2");
  TH2D * h_pvsp2 = (TH2D*)f->Get("pvsp2");
  TH2D * h_pmiss_thetamiss = (TH2D*)f->Get("pmiss_thetamiss");
  TH2D * h_compare = (TH2D*)f->Get("compare");
  TH2D * h_thetapn_dpp = (TH2D*)f->Get("thetapn_dpp");
  TH1D * h_tof2 = (TH1D*)f->Get("tof2");

  // neutron background selection



  // create output PDFs
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



  // electron selection
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts:");
  double line = 0.8;
  myText->Print(fileName,"pdf");  
  myText->Clear();


  // proton selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dbeta_p_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_vzp_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_chipid_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dbeta_p_cd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_vzp_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_chipid_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





  // neutron candidate selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_tof->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // neutron candidate selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_mmiss_pn->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // neutron signal selection
  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pxminuspx2->Draw("colz");
  myCanvas->cd(2);  h_pyminuspy2->Draw("colz");
  myCanvas->cd(3);  h_pzminuspz2->Draw();
  myCanvas->cd(4);  h_pminusp2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pvsp2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1); h_pmiss_thetamiss->Draw("colz");
  TLine * l_ptmiss1 = new TLine(40,0,40,1.2);
  l_ptmiss1->SetLineColor(kRed);
  l_ptmiss1->SetLineWidth(3);
  l_ptmiss1->Draw("same");
  TLine * l_ptmiss2 = new TLine(140,0,140,1.2);
  l_ptmiss2->SetLineColor(kRed);
  l_ptmiss2->SetLineWidth(3);
  l_ptmiss2->Draw("same");
  TLine * l_ptmiss3 = new TLine(0,0.25,180,0.25);
  l_ptmiss3->SetLineColor(kRed);
  l_ptmiss3->SetLineWidth(3);
  l_ptmiss3->Draw("same");
  TLine * l_ptmiss4 = new TLine(0,0.8,180,0.8);
  l_ptmiss4->SetLineColor(kRed);
  l_ptmiss4->SetLineWidth(3);
  //l_ptmiss4->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_compare->Draw("colz");
  TLine * l_comp1 = new TLine(-0.2,0,-0.2,180);
  l_comp1->SetLineColor(kRed);
  l_comp1->SetLineWidth(3);
  l_comp1->Draw("same");
  TLine * l_comp2 = new TLine(0.2,0,0.2,180);
  l_comp2->SetLineColor(kRed);
  l_comp2->SetLineWidth(3);
  l_comp2->Draw("same");
  TLine * l_comp3 = new TLine(-2,20,2,20);
  l_comp3->SetLineColor(kRed);
  l_comp3->SetLineWidth(3);
  l_comp3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_thetapn_dpp->Draw("colz");
  TLine * l_tdpp = new TLine(-2,60,2,60);
  l_tdpp->SetLineColor(kRed);
  l_tdpp->SetLineWidth(3);
  l_tdpp->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_tof2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


}
