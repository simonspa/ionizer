
// Daniel Pitzl, Sep 2019
// delta ray emission angle vs energy
// root -l emission.C

//------------------------------------------------------------------------------
//void emission()
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleOffset( 1.4, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(111111);
  gStyle->SetStatFormat("8.6g"); // more digits, default is 6.4g
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize(1); // no 'shadow'

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.90);

  gStyle->SetPalette(1); // rainbow colors

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );

  c1.Print( "emission.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // plot:

  TH1F * he = new
    TH1F( "he", "delta ray emission angle;delta ray energy [keV];emission angle [deg]",
	  999, 0, 999 );
  //he->SetNdivisions( -605, "x" );
  he->SetStats(kFALSE); // no statistics
  he->SetMinimum(0);
  he->SetMaximum(91);
  //he->SetNdivisions( -505, "y" );
  he->Draw();

  TLegend * lgnd = new TLegend( 0.50, 0.20, 0.93, 0.40 );

  TF1 * fMorris = new
    TF1( "fMorris", "57.296*acos(sqrt(x/(x+1022)))", 0, 999 );
  fMorris->SetNpx(1000);
  fMorris->SetLineWidth(2);
  fMorris->SetLineColor(kBlue+0);
  fMorris->Draw("samec"); // without axis option: overlay
  lgnd->AddEntry( fMorris, "Morris", "l" );

  TF1 * fGeant4h = new
    TF1( "fGeant4h", "57.296*acos(sqrt(x/(x+1022)*([0]+1022)/[0]))", 0, 999 );
  fGeant4h->SetParameter( 0, 30*1000 ); // incident kinetic energy [keV]
  fGeant4h->SetNpx(1000);
  fGeant4h->SetLineWidth(2);
  fGeant4h->SetLineColor(kBlack);
  fGeant4h->Draw("samec"); // without axis option: overlay
  lgnd->AddEntry( fGeant4h, "Geant4 (30 MeV)", "l" );

  TF1 * fGeant4 = new
    TF1( "fGeant4", "57.296*acos(sqrt(x/(x+1022)*([0]+1022)/[0]))", 0, 999 );
  fGeant4->SetParameter( 0, 3*1000 ); // incident kinetic energy [keV]
  fGeant4->SetNpx(1000);
  fGeant4->SetLineWidth(2);
  fGeant4->SetLineColor(kRed);
  fGeant4->Draw("samec"); // without axis option: overlay
  lgnd->AddEntry( fGeant4, "Geant4 (3 MeV)", "l" );

  TF1 * fBari = new // mccovpritms-orig.f
    TF1( "fBari", "57.296*asin(sqrt(1-x/[0]))", 0, 999 );
  fBari->SetParameter( 0, 3*1000 ); // incident kinetic energy [keV]
  fBari->SetNpx(1000);
  fBari->SetLineWidth(2);
  fBari->SetLineColor(kGreen+2);
  fBari->Draw("samec"); // without axis option: overlay
  lgnd->AddEntry( fBari, "Bari (3 MeV)", "l" );
  cout << fBari->Eval(500) << endl;
  lgnd->Draw( "same" );

  c1.Print( "emission.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "emission.pdf]" ); // ] closes file
  cout << "evince emission.pdf" << endl;

}
