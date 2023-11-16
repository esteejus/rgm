void print_ppipn_plots() {

  TFile * f = new TFile("all_ppipn.root","READ");

  char const * pdfFile = "plots_all_ppipn_pCD.pdf";


  // proton selection
  TH2D * h_pangles = (TH2D*)f->Get("pangles");
  TH2D * h_dbeta_p = (TH2D*)f->Get("dbeta_p");
  TH1D * h_vzp = (TH1D*)f->Get("vzp");
  //TH1D * h_chipid = (TH1D*)f->Get("chipid");

  // pion selection
  TH2D * h_dbeta_pi = (TH2D*)f->Get("dbeta_pi");
  TH1D * h_dvz_pi = (TH1D*)f->Get("dvz_pi");
  TH2D * h_piangles = (TH2D*)f->Get("piangles");

  // neutron candidate selection
  TH1D * h_tof = (TH1D*)f->Get("tof");
  TH2D * h_mmiss_pn = (TH2D*)f->Get("mmiss_pn");
  TH1D * h_nangles = (TH1D*)f->Get("nangles");
  TH2D * h_compare = (TH2D*)f->Get("andrew");

  // neutron signal selection
  TH1D * h_pxminuspx2 = (TH1D*)f->Get("pxminuspx2");
  TH1D * h_pyminuspy2 = (TH1D*)f->Get("pyminuspy2");
  TH1D * h_pzminuspz2 = (TH1D*)f->Get("pzminuspz2");
  TH1D * h_pminusp2 = (TH1D*)f->Get("pminusp2");
  TH2D * h_pvsp2 = (TH2D*)f->Get("pvsp2");
  TH2D * h_pmiss_thetamiss = (TH2D*)f->Get("pmiss_thetamiss");
  TH1D * h_nangles2 = (TH1D*)f->Get("nangles2");
  TH2D * h_compare2 = (TH2D*)f->Get("andrew2");
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
  myCanvas->cd(1);  h_dbeta_p->Draw("colz");
  TLine * l_betap1 = new TLine(0.5,-0.2,0.5,0.2);
  l_betap1->SetLineColor(kRed);
  l_betap1->SetLineWidth(3);
  l_betap1->Draw("same");
  TLine * l_betap2 = new TLine(0,-0.05,3,-0.05);
  l_betap2->SetLineColor(kRed);
  l_betap2->SetLineWidth(3);
  l_betap2->Draw("same");
  TLine * l_betap3 = new TLine(0,0.05,3,0.05);
  l_betap3->SetLineColor(kRed);
  l_betap3->SetLineWidth(3);
  l_betap3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_vzp->Draw();
  TLine * l_vzp1 = new TLine(-5,0,-5,h_vzp->GetMaximum());
  l_vzp1->SetLineColor(kRed);
  l_vzp1->SetLineWidth(3);
  l_vzp1->Draw("same");
  TLine * l_vzp2 = new TLine(5,0,5,h_vzp->GetMaximum());
  l_vzp2->SetLineColor(kRed);
  l_vzp2->SetLineWidth(3);
  l_vzp2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pangles->Draw("colz");
  TLine * l_pangles = new TLine(-180,40,180,40);
  l_pangles->SetLineColor(kRed);
  l_pangles->SetLineWidth(3);
  l_pangles->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // pion selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dbeta_pi->Draw("colz");
  TLine * l_betapi1 = new TLine(0.5,-0.2,0.5,0.2);
  l_betapi1->SetLineWidth(3);
  l_betapi1->SetLineColor(kRed);
  l_betapi1->Draw("same");
  TLine * l_betapi2 = new TLine(0,-0.05,3,-0.05);
  l_betapi2->SetLineWidth(3);
  l_betapi2->SetLineColor(kRed);
  l_betapi2->Draw("same");
  TLine * l_betapi3 = new TLine(0,0.05,3,0.05);
  l_betapi3->SetLineWidth(3);
  l_betapi3->SetLineColor(kRed);
  l_betapi3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dvz_pi->Draw("colz");
  TLine * l_vtz1 = new TLine(-4,0,-4,h_dvz_pi->GetMaximum());
  l_vtz1->SetLineWidth(3);
  l_vtz1->SetLineColor(kRed);
  l_vtz1->Draw("same");
  TLine * l_vtz2 = new TLine(4,0,4,h_dvz_pi->GetMaximum());
  l_vtz2->SetLineWidth(3);
  l_vtz2->SetLineColor(kRed);
  l_vtz2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_piangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // neutron candidate selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_tof->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_mmiss_pn->Draw("colz");
  TLine * l_mmiss = new TLine(0,1.05,1.5,1.05);
  l_mmiss->SetLineColor(kRed);
  l_mmiss->SetLineWidth(3);
  l_mmiss->Draw("same");
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
  TLine * l_ptmiss1 = new TLine(0,0.25,180,0.25);
  l_ptmiss1->SetLineColor(kRed);
  l_ptmiss1->SetLineWidth(3);
  l_ptmiss1->Draw("same");
  TLine * l_ptmiss2 = new TLine(0,1.25,180,1.25);
  l_ptmiss2->SetLineColor(kRed);
  l_ptmiss2->SetLineWidth(3);
  l_ptmiss2->Draw("same");
  TLine * l_ptmiss3 = new TLine(40,0,40,1.3);
  l_ptmiss3->SetLineColor(kRed);
  l_ptmiss3->SetLineWidth(3);
  l_ptmiss3->Draw("same");
  TLine * l_ptmiss4 = new TLine(140,0,140,1.3);
  l_ptmiss4->SetLineColor(kRed);
  l_ptmiss4->SetLineWidth(3);
  l_ptmiss4->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_compare->Draw("colz");
  TLine * l_compare1 = new TLine(-3,50,1,50);
  l_compare1->SetLineColor(kRed);
  l_compare1->SetLineWidth(3);
  l_compare1->Draw("same");
  TLine * l_compare2 = new TLine(-0.5,0,-0.5,180);
  l_compare2->SetLineColor(kRed);
  l_compare2->SetLineWidth(3);
  l_compare2->Draw("same");
  TLine * l_compare3 = new TLine(0.5,0,0.5,180);
  l_compare3->SetLineColor(kRed);
  l_compare3->SetLineWidth(3);
  l_compare3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nangles2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_compare2->Draw("colz");
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
