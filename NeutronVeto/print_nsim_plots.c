void print_nsim_plots() {

  //TFile *f = new TFile("test_p.root","READ");
  TFile * f = new TFile("all_nsim.root","READ");


  //char const * pdfFile = "plots_all_goodp_sim.pdf";
  char const * pdfFile = "plots_all_goodn_sim.pdf";


  // neutron candidate selection
  TH2D * h_nangles = (TH2D*)f->Get("nangles");
  TH2D * h_dpp_edep = (TH2D*)f->Get("dpp_edep");
  TH2D * h_edep_tof = (TH2D*)f->Get("edep_tof");
  TH2D * h_cos0 = (TH2D*)f->Get("cos0");
  TH1D * h_pxminuspx = (TH1D*)f->Get("pxminuspx");
  TH1D * h_pyminuspy = (TH1D*)f->Get("pyminuspy");
  TH1D * h_pzminuspz = (TH1D*)f->Get("pzminuspz");
  TH1D * h_pminusp = (TH1D*)f->Get("pminusp");
  TH2D * h_pvsp = (TH2D*)f->Get("pvsp");
  TH2D * h_dpp = (TH2D*)f->Get("dpp");
  TH1D * h_energy = (TH1D*)f->Get("energy");

  TH2D * h_nangles2 = (TH2D*)f->Get("nangles2");
  TH2D * h_dpp_edep2 = (TH2D*)f->Get("dpp_edep");
  TH2D * h_edep_tof2 = (TH2D*)f->Get("edep_tof");
  TH2D * h_cos02 = (TH2D*)f->Get("cos02");
  TH1D * h_pxminuspx2 = (TH1D*)f->Get("pxminuspx2");
  TH1D * h_pyminuspy2 = (TH1D*)f->Get("pyminuspy2");
  TH1D * h_pzminuspz2 = (TH1D*)f->Get("pzminuspz2");
  TH1D * h_pminusp2 = (TH1D*)f->Get("pminusp2");
  TH2D * h_pvsp2 = (TH2D*)f->Get("pvsp2");
  TH2D * h_dpp2 = (TH2D*)f->Get("dpp2");
  TH1D * h_energy2 = (TH1D*)f->Get("energy2");



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


  // neutron candidate selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nangles->Draw("colz");
  TLine * line_nangles1 = new TLine(-180,40,180,40);
  line_nangles1->SetLineColor(kRed);
  line_nangles1->SetLineWidth(3);
  line_nangles1->Draw("same");
  TLine * line_nangles2 = new TLine(-180,140,180,140);
  line_nangles2->SetLineColor(kRed);
  line_nangles2->SetLineWidth(3);
  line_nangles2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dpp_edep->Draw("colz");
  TLine * line_edep = new TLine(1.5,-0.4,1.5,0.4);
  line_edep->SetLineColor(kRed);
  line_edep->SetLineWidth(3);
  line_edep->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_edep_tof->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_cos0->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pxminuspx->Draw();
  myCanvas->cd(2);  h_pyminuspy->Draw();
  myCanvas->cd(3);  h_pzminuspz->Draw();
  myCanvas->cd(4);  h_pminusp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pminusp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


/*  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pminusp->Draw();
  myCanvas->cd(2);  h_pminusp_p->Draw();
  myCanvas->cd(3);  h_pminusp_theta->Draw();
  myCanvas->cd(4);  h_pminusp_phi->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
*/


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dpp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_energy->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // neutron/proton selection
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nangles2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dpp_edep2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_edep_tof2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_cos02->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pxminuspx2->Draw();
  myCanvas->cd(2);  h_pyminuspy2->Draw();
  myCanvas->cd(3);  h_pzminuspz2->Draw();
  myCanvas->cd(4);  h_pminusp2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_dpp2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_energy2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


}
