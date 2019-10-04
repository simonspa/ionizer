
// rlq mfp.root
// .x drawf.C
{
  xmfpvse->Draw("p");
  gPad->SetLogy();
  TF1 * f1 = new TF1( "f1", "0.26*pow(10,x)/(0.15+pow(10,x))", -3, 3 );
  f1->Draw("same");
}
