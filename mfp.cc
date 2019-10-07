
// calculate mean free path for ionization
// Hans Bichsel code

// make mfp
// mfp
// reads HEPS.TAB, MACOM.TAB, EMERC.TAB

// produces mfp.root

#include <cstdlib> // atoi
#include <iostream> // cout
#include <fstream> // files
#include <sstream> // stringstream
#include <cmath> // log
#include <random>
#include <ctime>
#include <stack>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>

using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
{
  unsigned npm = 4; // e

  double temp = 300; // [K]

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  double zi = 1;
  double elmm = 0.51099906; // e mass [MeV]

  double ptm = elmm; // e 0.51100 MeV
  if(      npm == 1 ) ptm = 938.2723; // proton
  else if( npm == 2 ) ptm = 139.578; // pion
  else if( npm == 3 ) ptm = 493.67; // K
  else if( npm == 5 ) ptm = 105.65932; // mu

  double elm = 1e6 * elmm; // me [eV]
  double pi = 3.1415926536;
  double Ry = 13.6056981;
  double fac = 8 * pi * Ry*Ry * pow( 0.529177e-8, 2 ) / elm;
  double log10 = log(10);

  // silicon:

  double Z = 14.0; // atomic number of absorber, Si
  double A = 28.086; // atomic weight of absorber
  double rho = 2.329; // rho= density of absorber material
  double radl = 9.36; // [cm]

  double atnu = 6.0221367e23 * rho / A; // atnu = # of atoms per cm**3

  cout
    << endl << "particle type " << npm << ", mass " << ptm << " MeV"
    << endl << "element " << Z << ", density " << rho << " g/cm3"
    << endl << "temperature " << temp << " K"
    << endl;

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // book histos

  TFile * histoFile = new TFile( "mfp.root", "RECREATE" );

  TProfile betavse( "betavse",
		    "beta;log_{10}(Ekin [MeV]);beta",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile lastvse( "lastvse",
		    "last dE/E;log_{10}(Ekin [MeV]);last de/E",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile zeffvse( "zeffvse",
		    "Zeff;log_{10}(Ekin [MeV]);Zeff",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile xmfpvse( "xmfpvse",
		    "inelastic mean free path;log_{10}(Ekin [MeV]);inelastic mean free path [#mum]",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile derdxvse( "derdxvse",
		     "restricted energy loss;log_{10}(Ekin [MeV]);restricted energy loss [eV/um]",
		     6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile dedxvse( "dedxvse",
		    "total energy loss;log_{10}(Ekin [MeV]);total energy loss [eV/um]",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile rm0vse( "rm0vse",
		   "tail interactions;log_{10}(Ekin [MeV]);tail interactions [1/#mum]",
		   6*50+1, -3.01, 3.01 ); // 50 steps/decade
  TProfile devse( "devse",
		  "energy loss per step;log_{10}(Ekin [MeV]);energy loss [eV/interaction]",
		  6*50+1, -3.01, 3.01 ); // 50 steps/decade

  TProfile emfpvse( "emfpvse",
		    "elastic mean free path;log_{10}(Ekin [MeV]);elastic mean feee path [#mum]",
		    6*50+1, -3.01, 3.01 ); // 50 steps/decade

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // initialize energy bins

  int n2    = 64;
  double u  = log(2.0) / n2;
  double um = exp(u);
  int ken  = log( 1839.0 / 1.5 ) / u; // intger
  double Emin = 1839.0 / pow( 2, ken/n2 ); // integer division intended

  // EMIN is chosen to give an E-value exactly at the K-shell edge, 1839 eV

  const unsigned lime = 1251; // see HEPS.TAB
  unsigned nume = lime-1;

  double E[lime], dE[lime];

  E[1] = Emin;

  for( unsigned j = 1; j < nume; ++j ) {
    E[j+1] = E[j] * um; // must agree with heps.tab
    dE[j]  = E[j+1] - E[j];
  }

  double Etop = E[nume]*sqrt(um); // [eV]

  cout
    << endl << "dE bins " << n2*log10 << " per decade, step factor " << um
    << endl << "Emin " << Emin << " eV"
    << endl << "Etop " << Etop*1e-6 << " MeV"
    << endl << "bins " << nume
    << endl;

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // read dielectric constants

  double ep[3][lime];
  double dfdE[lime];

  ifstream f14( "HEPS.TAB" );

  if( f14.bad() || ! f14.is_open() ) {
    cout << "Error opening HEPS.TAB" << endl;
    return 1;
  }
  else {

    string line;
    getline( f14, line );

    istringstream tokenizer( line );

    int n2t;
    unsigned numt;
    tokenizer >> n2t;
    tokenizer >> numt;

    cout << endl << "HEPS.TAB: n2t " << n2t << ", numt " << numt << endl;
    if( n2 != n2t ) cout << "  CAUTION: n2 & n2t differ" << endl;
    if( nume != numt ) cout << "  CAUTION: nume & numt differ" << endl;
    if( numt > nume ) numt = nume;

  // HEPS.TAB is the table of the dielectric constant for solid Si,
  // epsilon = ep(1,j) + i*ep(2,j), as a function of energy loss E(j),
  // section II.E in RMP, and rim is Im(-1/epsilon), Eq. (2.17), p.668.
  // Print statements are included to check that the file is read correctly.

    unsigned jt = 1;

    while( ! f14.eof() && jt < numt ) {

      getline( f14, line );

      istringstream tokenizer( line );

      double etbl, ep1, ep2, rimt;
      tokenizer >> jt >> etbl >> ep1 >> ep2 >> rimt;

      //cout << jt << "  " << etbl << "  " << ep1 << "  " << ep2 << "  " << rimt << endl;

      ep[1][jt] = ep1;
      ep[2][jt] = ep2;

     // The dipole oscillator strength df/dE is calculated,
     // essentially Eq. (2.20)

      dfdE[jt] = rimt * 0.0092456 * E[jt];

    }

    cout << "  read " << jt << " data lines from HEPS.TAB" << endl;

    // Mazziotta: 0.0 at 864
    // EP( 2, 864 ) = 0.5 * ( EP(2, 863) + EP(2, 865) )
    // RIM(864) = 0.5 * ( RIM(863) + RIM(865) )
    // DFDE(864) = RIM(864) * 0.0092456 * E(864)
    // DP: fixed in HEPS.TAB

  } // HEPS.TAB

  //= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // read integral over momentum transfer of the generalized oscillator strength

  double sig[7][lime];

  ifstream f15( "MACOM.TAB" );

  if( f15.bad() || ! f15.is_open() ) {
    cout << "Error opening MACOM.TAB" << endl;
    return 1;
  }
  else {

    string line;
    getline( f15, line );

    istringstream tokenizer( line );

    int n2t;
    unsigned numt;
    tokenizer >> n2t;
    tokenizer >> numt;

    cout << endl << "MACOM.TAB: n2t " << n2t << ", numt " << numt << endl;
    if( n2 != n2t ) cout << "  CAUTION: n2 & n2t differ" << endl;
    if( nume != numt ) cout << "  CAUTION: nume & numt differ" << endl;
    if( numt > nume ) numt = nume;

  // MACOM.TAB is the table of the integrals over momentum transfer K of the
  // generalized oscillator strength, summed for all shells, i.e. the A(E)
  // of Eq. (2.11), p. 667 of RMP

    unsigned jt = 1;

    while( ! f15.eof() && jt < numt ) {

      getline( f15, line );

      istringstream tokenizer( line );

      double etbl, sigt;
      tokenizer >> jt >> etbl >> sigt;

      //cout << jt << "  " << etbl << "  " << sigt << endl;

      sig[6][jt] = sigt;

    }

    cout << "  read " << jt << " data lines from MACOM.TAB" << endl;

  } // MACOM.TAB

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  double xkmn[200];

  ifstream f16( "EMERC.TAB" );

  if( f16.bad() || ! f16.is_open() ) {
    cout << "Error opening EMERC.TAB" << endl;
    return 1;
  }
  else {

    string line;
    getline( f16, line ); // header lines
    getline( f16, line );
    getline( f16, line );
    getline( f16, line );

    cout << endl << "EMERC.TAB" << endl;

    // EMERC.TAB is the table of the integral over k of generalized oscillator
    // strength for E < 11.9 eV with Im(-1/epsilon) from equations in the Appendix
    // of Emerson et al., Phys Rev B7, 1798 (1973) (also see CCS-63)

    unsigned jt = 1;

    while( ! f16.eof() && jt < 200 ) {

      getline( f16, line );

      istringstream tokenizer( line );

      double etbl, sigt, xk;
      tokenizer >> jt >> etbl >> sigt >> xk;

      //cout << jt << "  " << etbl << "  " << sigt << "  " << xk << endl;

      sig[6][jt] = sigt; // overwritten!
      xkmn[jt] = xk;

    }

    cout << "  read " << jt << " data lines from EMERC.TAB" << endl;

  } // EMERC.TAB

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  // incident energy:

  double dErdx = 0;
  double dEdx = 0;
  double xmfp = 0;

  double ekmev = 0.001; // [MeV]

  while( ekmev < 1001 ) {

    double l10ekmev = log(ekmev)/log10;
    double gam = ekmev / ptm + 1; // W = total energy / restmass
    double bg  = sqrt( gam*gam - 1 ); // bg = beta*gamma
    double betasq = bg*bg / ( 1 + bg*bg );
    betavse.Fill( l10ekmev, sqrt(betasq) );

    // maximum energy loss, see Uehling, also Sternheimer & Peierls Eq.(53)

    double Emax = ptm * ( gam*gam - 1 ) / ( 0.5*ptm/elmm + 0.5*elmm/ptm + gam ); // bug fixed
    if( npm == 4 ) // e
      Emax = 0.5*ekmev;

    Emax = 1e6 * Emax; // eV

    // Define parameters and calculate Inokuti"s sums,
    // Sect 3.3 in Rev Mod Phys 43, 297 (1971)

    double dec = zi*zi * atnu * fac / betasq;
    double bbemx = betasq / Emax; // [1/eV]
    double ekin  = ekmev * 1e6;
    double twombb = 2 * elm * betasq; // [eV]

    // Generate collision spectrum sigma(E) from df/dE, epsilon and AE.
    // sig(*,j) actually is dE**2 * sigma(E)

    double tsig[6];
    for( unsigned i = 1; i <= 5; ++i )
      tsig[i] = 0;

    double zeff = 0;
    double stpw = 0;
    unsigned nlast = nume;

    for( unsigned j = 1; j <= nume; ++j ) {

      if( E[j] > Emax ) break;

      nlast = j;

      zeff += dfdE[j] * dE[j];

      // Eq. (3.1) in RMP and red notebook CCS-33, 39 & 47

      double Q1 = Ry;
      if( E[j] < 100.0 ) Q1 = pow( 0.025, 2 ) * Ry;
      if( E[j] <  11.9 ) Q1 = pow( xkmn[j], 2 ) * Ry;

      double qmin = E[j]*E[j] / twombb; // twombb = 2 m beta**2 [eV]

      if( E[j] < 11.9 && Q1 < qmin )
	sig[1][j] = 0;
      else
	sig[1][j] = E[j] * dfdE[j] * log( Q1 / qmin );
      // longitudinal excitation, Eq. (46) in Fano; Eq. (2.9) in RMP

      double epbe = 1 - betasq * ep[1][j]; // Fano Eq. (47)
      if( epbe < 1e-20 ) epbe = 1E-20;

      double sgg = E[j] * dfdE[j] * (-0.5) *
	log( epbe*epbe + pow( betasq * ep[2][j], 2 ) );

      double thet = atan( ep[2][j] * betasq / epbe );
      if( thet < 0 ) thet = thet + pi; // plausible-otherwise I"d have a jump
      // Fano says [p 21]: "arctan approaches pi for betasq*eps1 > 1 "

      double sgh = 0.0092456 * E[j]*E[j] * thet *
	( betasq - ep[1][j] / ( pow( ep[1][j], 2 ) + pow( ep[2][j], 2 ) ) );

      sig[3][j] = sgg + sgh;

      double uef = 1 - E[j] * bbemx;
      if( npm == 4 )
	uef = 1 +
	  pow( E[j] / ( ekin - E[j] ), 2 ) +
	  pow( (gam-1) / gam * E[j]/ekin, 2 ) -
	  ( 2*gam-1 ) * E[j] / ( gam*gam * ( ekin - E[j] ) );

      // uef from  Eqs. 9 & 2 in Uehling, Ann Rev Nucl Sci 4, 315 (1954)
      // if( j == 1) PRINT*, " uef=",UEF

      sig[4][j] = 2 * sig[6][j] * uef;

      // there is a factor of 2 because the integral was over d(lnK) rather than d(lnQ)

      sig[2][j] = 0;
      sig[5][j] = 0;

      for( unsigned i = 1; i <= 4; ++i ) {

	// integrated total collision cross section:

	tsig[i]  = tsig[i]  + sig[i][j] * dE[j] / ( E[j]*E[j] );

	// sig(5,j] = total differential cross section, Eq. (3.8) in RMP

	sig[5][j] += sig[i][j];

      } // i

      tsig[5] += sig[5][j] / ( E[j]*E[j] ) * dE[j]; // total cross section

      stpw += sig[5][j] / E[j] * dE[j]; // dE/dx

    } // dE j

    double xm0 = tsig[5] * dec; // 1/path [1/cm]

    dErdx = stpw * dec; // [eV/cm]

    double Efin = E[nlast]*sqrt(um); // [eV]

    lastvse.Fill( l10ekmev, Efin/ekin );
    zeffvse.Fill( l10ekmev, zeff );
    xmfpvse.Fill( l10ekmev, 1e4/xm0 ); // [um]
    derdxvse.Fill( l10ekmev, dErdx*1e-4 ); // [eV/um]

    // add rest from Etop to Emax:

    double rst = 0;
    double rm0 = 0;

    if( Emax > Etop ) {

      if( npm == 4 ) // e, FSR-143 and Uehling Eq 9
	rst = // dEdx
	  log( Emax/Etop ) + log( Emax ) - log( ekin - Etop ) - 1 / ( 1 - Etop / ekin ) + 2 +
	  pow( ( gam-1 ) / gam, 2 ) * ( 0.125 - 0.5 * pow( Etop/ekin, 2 ) ) +
	  ( 2*gam - 1 ) / (gam*gam) * ( log( Emax ) - log( ekin-Etop ) );
      else
	rst =
	  ( 1 - 720*betasq / Emax ) * log( Emax/Etop ) +
	  720 * ( 1/Etop - 1/Emax ) - betasq * ( 1 - Etop/Emax );

      rm0 = // cross section
	( 1 - 720*betasq / Emax ) *
	( 1/Etop - 1/Emax + 2/(Etop*Etop) - 2/(Emax*Emax) ) -
	betasq * log( Emax/Etop ) / Emax;

    }

    dEdx = dErdx + Z * rst * dec;
    dedxvse.Fill( l10ekmev, dEdx*1e-4 ); // [eV/um]

    rm0vse.Fill( l10ekmev, Z*rm0*dec*1e-4 ); // [1/um] should be small, is small
    xmfp = xm0 + Z*rm0*dec; // [1/cm]
    devse.Fill( l10ekmev, dEdx/xmfp ); // [eV/step] 40..120

    // elastic:

    double xlel = 1;

    if( npm == 4 ) { // e

      double gn = 2*2.61 * pow( Z, 2.0/3.0 ) / ekin;
      double E2 = 14.4e-14; // [MeV*cm]
      double FF = 0.5*pi * E2*E2 * Z*Z / (ekmev*ekmev);
      double S0EL = 2*FF / ( gn * ( 2 + gn ) );
      // elastic total cross section  [cm2/atom]
      xlel = atnu*S0EL; // ATNU = N_A * rho / A = atoms/cm3

    }
    else { // other particles

      double getot = ekmev + ptm;
      double pmom = sqrt( ekmev * ( getot + ptm ) );
      xlel = min( 2232 * radl * pow( pmom*pmom / (getot*zi), 2 ), 10*radl );
      // units ?
    }

    emfpvse.Fill( l10ekmev, 1e4/xlel ); // [um]

    ekmev *= pow( 10, 0.02 );

  } // ekmev

  cout << endl << "final mfp  " << 1e4/xmfp << " um"
       << endl << "final dEdx " << dEdx*1e-4 << " eV/um"
       << endl << "restricted " << dErdx*1e-4 << " eV/um"
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  histoFile->Write();
  cout << endl;
  histoFile->ls();
  histoFile->Close();
  cout << endl
       << histoFile->GetName() << endl;
  cout << endl;

} // main
