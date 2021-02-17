
// simulation of ionization by charged particle tracks in silicon
// z = along track: columns
// x = transverse: rows, turn angle in z-x
// y = vertical = drift

// time ionizer -n 10100 -p 25 -d 285 -t 500 -a 9.5 -e 5000 -c 0.02

// -n events
// -p pixel width [mu]
// -d pixel Dicke [mu]
// -t pixel threshold [e]
// -c cross talk [fraction]
// -a angle of incidence [deg] default is ideal 2-pix
// -e kinetic energy [MeV]
// -f not fast

// History:

// COVPRI.f calculates the primary collision cross section table for Si
// Hans Bichsel, Seattle, deceased 2018
// Reference: Rev Mod Phys 60, 663 (1988)

// MCCOVPRI.f to make the energy loss distribution
// (straggling) in Si, using tBichsel's ionization cross section
// and to generate e-h pairs with shell transitions
// M. Nicola Mazziotta, INFN Sezione di Bari 2003
// Reference: NIM A533 (2004) 322

// MCCOVPRITMS.f to include the electron elastic cross section
// int the Si, using the Z.E.A. Chaoui prescriptions
// and temperature dependence of the band gap
// M. Nicola Mazziotta, INFN Sezione di Bari 2007
// mazziotta@ba.infn.it
// References:
// Nuclear Instruments and Methods in Physics Research A584 (2008) 436–439
// Surface and Interface Analysis 38 (2006) 664
// Applied Physics Letter 88 (2006) 024105
// Physics Letters A297 (2002) 432

// converted to C++ by Daniel Pitzl, DESY, Sep 2019

// Input files:
// HEPS.TAB
// MACOM.TAB
// EMERC.TAB

// Output file:
// ionizer.hist

#include <cstdlib> // atoi
#include <iostream> // cout
#include <fstream> // files
#include <sstream> // stringstream
#include <cmath> // log
#include <random>
#include <ctime>
#include <stack>

#include <TFile.h> // ROOT
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

using namespace std;

struct delta {
  double E; // [MeV]
  double x; // position
  double y;
  double z;
  double u; // direction
  double v;
  double w;
  unsigned npm; // particle type
};

struct cluster {
  int neh;
  double x; // position
  double y;
  double z;
  double E; // [eV] generating particle
};

// forward declarations:

double alph1( double x );
double alph2( double x );
double trian( double x );
double gena1();
double gena2();
double gentri();
void shells( double Eg, stack <double> &veh );
void TRL1MM( double Ev, stack <double> &veh );
void TRL23MM( double Ev, stack <double> &veh );

// global variables: (initialized in main, used in shells)

const unsigned lsh = 5;
double Eshel[lsh], augmi[lsh][10], augde[lsh][10];
int nvac[lsh];

const unsigned lep = 14;
double EPP[lep], PM[lep], PL23[lep], PL1[lep], PK[lep];
const unsigned nep = lep-1;

ranlux24 rgen; // C++11 random number engine
uniform_real_distribution <double> unirnd( 0, 1 );

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main( int argc, char* argv[] )
{
  // defaults:

  bool ldb = 0; // debug flag

  unsigned nev = 10*1000;

  double depth = 285; // [mu] pixel depth

  double pitch = 25; // [mu] pixels size

  double angle = 999; // flag

  double thr = 500; // threshold [e]

  double cx = 0; // cross talk

  double Ekin0 = 5000; // [MeV] kinetic energy

  bool fast = 1; // default is fast

  uint64_t seed = 0;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-p" ) )
      pitch = atof( argv[++i] ); // [mu]

    if( !strcmp( argv[i], "-d" ) )
      depth = atof( argv[++i] ); // [mu]

    if( !strcmp( argv[i], "-a" ) )
      angle = atof( argv[++i] ); // [deg]

    if( !strcmp( argv[i], "-t" ) )
      thr = atof( argv[++i] ); // pixel threshold [e]

    if( !strcmp( argv[i], "-c" ) )
      cx = atof( argv[++i] ); // cross talk fraction

    if( !strcmp( argv[i], "-e" ) )
      Ekin0 = atof( argv[++i] ); // [MeV]

    if( !strcmp( argv[i], "-n" ) )
      nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      fast = 0; // full ionization, not fast: simulate each e-h pair

    if( !strcmp( argv[i], "-s" ) )
      seed = atoi( argv[++i] ); // random seed

    if( !strcmp( argv[i], "-h" ) ) {
      cout
	<< "  simulate ionization by charged particles in silicon" << endl
	<< "  usage: ionizer [option] [option] [option]" << endl
	<< "         produces ionizer.hist" << endl
	<< "  options:" << endl
	<< "    -z pixel length [mu] (default 150)" << endl
	<< "    -p pixel width [mu] (default 25)" << endl
	<< "    -d pixel depth [mu] (default 285)" << endl
	<< "    -a angle of incidence [deg] (default atan(p/d))" << endl
	<< "    -t readout threshold [e] (default 500)" << endl
	<< "    -c cross talk fraction (default 0)" << endl
	<< "    -e incident kinetic energy [MeV] (default 5000)" << endl
	<< "    -f full ionization (slow, default neh=de/3.645)" << endl
	<< "    -n number of events (default 10000)" << endl
	;
      return 0;
    }

  } // argc

  double pi = 3.1415926536;
  double wt = 180/pi;
  double twopi = 2*pi;
  double w2 = sqrt(2);

  double turn = atan( pitch / depth ); // [rad] default
  if( fabs(angle) < 91 )
    turn = angle/wt;

  double width = depth*tan(turn); // [mu] projected track, default: pitch

  // [V/cm] mean electric field: Vbias-Vdepletion/2
  double Efield = (120-30)/depth*1e4; // UHH
  //double Efield = (70-60)/depth*1e4; // DepFET

  // delta ray range: 1 um at 10 keV (Mazziotta 2004)
  //double explicit_delta_energy_cut_keV = 2; Dec 2019
  double explicit_delta_energy_cut_keV = 9; // Apr 2020, faster, no effect on resolution
  //double explicit_delta_energy_cut_keV = 99; // no effect on resolution

  // p=1, pi=2, K=3, e=4, mu=5, He=6, Li=7, C=8, Fe=9

  unsigned npm0 = 4; // e

  double temp = 298; // [K]

  double elmm = 0.51099906; // e mass [MeV]
  double elm = 1e6 * elmm; // me [eV]
  double twome = 2*elm; // [eV]
  double Ry = 13.6056981;
  double fac = 8.0 * pi * Ry*Ry * pow( 0.529177e-8, 2 ) / elm;
  double log10 = log(10);

  cout << "  particle type     " << npm0 << endl;
  cout << "  kinetic energy    " << Ekin0 << " MeV" << endl;
  cout << "  number of events  " << nev << endl;
  cout << "  pixel pitch       " << pitch << " um" << endl;
  cout << "  pixel depth       " << depth << " um" << endl;
  cout << "  incident angle    " << turn*wt << " deg" << endl;
  cout << "  track width       " << width << " um" << endl;
  cout << "  temperature       " << temp << " K" << endl;
  cout << "  readout threshold " << thr << " e" << endl;
  cout << "  cross talk        " << cx*100 << "%" << endl;

  // mobility from pixelav:
  // 0 = e, 1 = h

  int j = 0; // e CMS, B2 pixel
  //int j = 1; // h for strips

  double cvm[2] = { 1.53e9, 1.62e8 }; // [cm/s] vmax at temp=1K
  double evm[2] = { -0.87, -0.52 };
  double vm = cvm[j] * pow( temp, evm[j] );

  double cec[2] = { 1.01, 1.24 }; // [V/cm] Ecrit
  double eec[2] = { 1.55, 1.68 };
  double Ec = cec[j] * pow( temp, eec[j] );
  double mu0 = vm / Ec;

  double cbeta[2] = { 0.0257, 0.46 };
  double ebeta[2] = { 0.66, 0.17 };
  double beta = cbeta[j] * pow( temp, ebeta[j] );
  double ibeta = 1 / beta;

  double d2 = Efield / Ec;
  double d3 = pow( d2, beta ) + 1.;
  double mu = mu0 / pow( d3, ibeta ); // mu0 / ( 1 + (E/Ec)^b )^(1/b)
  const double vd = Efield*mu; // [cm/s]

  // diffusion from mobility: D = kTmu/e
  // e = 1.602e-19 C
  // k = 1.38e-23 J/K
  // k/e = 8.6e-5 eV/K

  const double D = 8.61733e-5 * temp * mu; // diffuson constant

  cout
    << endl
    << "   mobility for " << Efield << " V/cm"
    << ": vm " << vm // cm/s = 100 um / ns
    << ", Ec " << Ec
    << ", mu0 " << mu0 << endl
    << "  beta " << beta
    << ", mu " << mu
    << ", v " << vd << " cm/s"
    << " = " << vd/1e5 << " mu/ns" << endl
    << "  D " << D
    << ", rms " << sqrt(2*D*4e-9)*1e4 << " mu" // for 4 ns drift
    << endl;

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // book histos

  TFile * histoFile = new
    TFile( Form( "ionizer_p%i_w%i_d%i_t%i_c%i.hist",
		 int(pitch+0.5), int(width+0.5), int(depth+0.5),
		 int(thr), int(100*cx+0.1) ),
	   "RECREATE" );

  // book histos:

  TProfile elvse( "elvse", "elastic mfp;log_{10}(E_{kin}[MeV]);elastic mfp [#mum]",
		  140, -3, 4 );
  TProfile invse( "invse", "inelastic mfp;log_{10}(E_{kin}[MeV]);inelastic mfp [#mum]",
		  140, -3, 4 );

  TH1I hstep5( "step5", "step length;step length [#mum];steps", 500, 0, 5 );
  TH1I hstep0( "step0", "step length;step length [#mum];steps", 500, 0, 0.05 );
  TH1I hzz( "zz", "z;depth z [#mum];steps", depth, 0, depth );

  TH1I hde0( "de0", "step E loss;step E loss [eV];steps", 200, 0, 200 );
  TH1I hde1( "de1", "step E loss;step E loss [eV];steps", 100, 0, 5000 );
  TH1I hde2( "de2", "step E loss;step E loss [keV];steps", 200, 0, 20 );
  TH1I hdel( "del", "log step E loss;log_{10}(step E loss [eV]);steps", 140, 0, 7 );
  TH1I htet( "tet", "delta emission angle;delta emission angle [deg];inelasic steps",
	     180, 0, 90 );
  TH1I hnprim( "nprim", "primary eh;primary e-h;scatters", 21, -0.5, 20.5 );
  TH1I hlogE( "logE", "log Eeh;log_{10}(E_{eh}) [eV]);eh", 140, 0, 7 );
  TH1I hlogn( "logn", "log neh;log_{10}(n_{eh});clusters",  80, 0, 4 );

  TH1I hscat( "scat", "elastic scattering angle;scattering angle [deg];elastic steps",
	      180, 0, 180 );

  TH1I hncl( "ncl", "clusters;e-h clusters;tracks", 4*depth*5, 0, 4*depth*5 );

  double lastbin = 5*0.35*depth; // 350 eV/micron
  if( Ekin0 < 1.1 )
    lastbin = 1.05*Ekin0*1e3; // [keV]
  TH1I htde( "tde", "sum E loss;sum E loss [keV];tracks / keV",
	     max(100,int(lastbin)), 0, int(lastbin) );
  TH1I htde0( "tde0", "sum E loss, no delta;sum E loss [keV];tracks, no delta",
	      max(100,int(lastbin)), 0, int(lastbin) );
  TH1I htde1( "tde1", "sum E loss, with delta;sum E loss [keV];tracks, with delta",
	      max(100,int(lastbin)), 0, int(lastbin) );

  TH1I hteh( "teh", "total e-h;total charge [ke];tracks",
	     max(100,int(50*0.1*depth)), 0, max(1,int(10*0.1*depth)) );
  TH1I hq0( "q0", "normal charge;normal charge [ke];tracks",
	     max(100,int(50*0.1*depth)), 0, max(1,int(10*0.1*depth)) );
  TH1I hrms( "rms", "RMS e-h;charge RMS [e];tracks",
	     100, 0, 50*depth );

  TH1I * h1zev[11];
  TH2I * h2zxev[11];
  for( unsigned i = 0; i < 11; ++i ) {
    h1zev[i] = new
      TH1I( Form( "z%02i", i ),
	    Form( "z event %i;z [#mum];clusters [eh-pairs]", i ),
	    4*depth, 0, depth );
    h2zxev[i] = new
      TH2I( Form( "zx%02i", i ),
	    Form( "z-x event %i;x [#mum];z [#mum];clusters [eh-pairs]", i ),
	    4*2*pitch, -pitch, pitch, 4*depth, 0, depth );
  }

  TH2I * h2xy = new
    TH2I( "xy","x-y clusters;x [#mum];y [#mum];clusters [eh-pairs]",
	  400, -200, 200, 400, -200, 200 );
  TH2I * h2zx = new
    TH2I( "zx","z-x clusters;x [#mum];z [#mum];clusters [eh-pairs]",
	  4*pitch, -2*pitch, 2*pitch, depth, 0, depth );

  TH1I hdtime( "dtime", "drift time;drift time [ns];clusters", 100, 0, 10+20*j );
  TH1I hdiff( "diff", "diffusion width;diffusion width [#mum];clusters", 100, 0, 10 );
  TH1I htf( "tf", "Gaussian tail fraction;Gausian tail fraction;clusters", 100, 0, 1 );
  TProfile tfvsx( "tfvsx", "Gaussian tail fraction vs x;x [#mum];Gausian tail fraction",
		  200, -0.5*pitch, 0.5*pitch );

  TH1I hcleh( "cleh", "cluster neh;log_{10}(cluster eh [pairs]);clusters", 80, 0, 4 );
  TProfile wvse( "wvse", "energy per eh pair;log_{10}(step E loss [eV]);<w> [eV/pair]",
		 80, 0, 4 );
  TH1I hreh( "reh", "eh/eV;eh/dE [pairs/eV];clusters", 160, 0, 0.8 );

  TH1I heta0( "eta0", "eta;eta;tracks", 201, -1.005, 1.005 );
  TProfile eta0vsxm( "eta0vsxm", "eta vs track;track x [#mum];<eta>",
		     200, -0.5*pitch, 0.5*pitch );
  TH1I hdx0( "dx0", "dx0;#Deltax [#mum];tracks",
	     501, -pitch*1.001, pitch*1.001 );
  TProfile madx0vsq( "madx0vsq",
		     "MAD(#Deltax) vs charge;charge [ke];MAD(#Deltax) [#mum]",
		     100, 0, 2*0.1*depth );
  TH1I hdx0q( "dx0q", "dx0 Landau peak;#Deltax [#mum];tracks",
	      501, -pitch*1.001, pitch*1.001 );
  TProfile dx0qvsxm( "dx0qvsxm", "#Deltax vs x;track x [#mum];<#Deltax> [#mum]",
		     pitch, -0.5*pitch, 0.5*pitch );
  TProfile
    madx0qvsxm( "madx0qvsxm", "MAD(#Deltax) vs x;track x [#mum];MAD(#Deltax) [#mum]",
		       pitch, -0.5*pitch, 0.5*pitch );
  TProfile
    madx0qvsxm0( "madx0qvsxm0",
		 "MAD(#Deltax) vs x no delta;track x [#mum];MAD(#Deltax) [#mum]",
		 pitch, -0.5*pitch, 0.5*pitch );
  TProfile
    madx0qvsxm1( "madx0qvsxm1",
		 "MAD(#Deltax) vs x delta;track x [#mum];MAD(#Deltax) [#mum]",
		 pitch, -0.5*pitch, 0.5*pitch );

  // threshold:

  TH1I hpxq1( "pxq1", "thresholded pixel charge;pixel charge [ke];pixels",
	      max(100,int(10*0.1*depth/1)), 0, max(1,int(5*0.1*depth/1)) );
  TH1I hq1( "q1", "thresholded charge;charge [ke];tracks",
	    max(100,int(50*0.1*depth)), 0, max(1,int(10*0.1*depth)) );
  TH1I hnpx1( "npx1", "npx after threshold;npx;tracks", 4, 0.5, 4.5 );
  TProfile npx1vsxm( "npx1vsxm", "npx threshold vs track;track x [#mum];<npx>",
		     200, -0.5*pitch, 0.5*pitch );
  TH1I heta1( "eta1", "eta threshold;eta;tracks", 201, -1.005, 1.005 );
  TProfile eta1vsxm( "eta1vsxm", "eta threshold vs track;track x [#mum];<eta>",
		     200, -0.5*pitch, 0.5*pitch );
  TProfile x1vsxm( "x1vsxm", "xcog threshold vs track;track x [#mum];<cog> [#mum]",
		   200, -0.5*pitch, 0.5*pitch );
  TH1I hdx1( "dx1",
	     "dx threshold;#Deltax [#mum];tracks", 501, -pitch*1.001, pitch*1.001 );
  TProfile madx1vsq( "madx1vsq",
		     "MAD(#Deltax) vs charge;charge [ke];MAD(#Deltax) [#mum]",
		     100, 0, 2*0.1*depth );
  TH1I hdx1q( "dx1q",
	      "dx threshold Landau peak;#Deltax [#mum];tracks",
	      501, -pitch*1.001, pitch*1.001 );
  TProfile dx1qvsxm( "dx1qvsxm", "#Deltax vs x;track x [#mum];<#Deltax> [#mum]",
		     pitch, -0.5*pitch, 0.5*pitch );
  TProfile madx1qvsxm( "madx1qvsxm",
		       "MAD(#Deltax) vs x;track x [#mum];MAD(#Deltax) [#mum]",
		       pitch, -0.5*pitch, 0.5*pitch );

  TH1I hda1( "da1", "da threshold;#Deltax [#mum];tracks",
	     501, -pitch*1.001, pitch*1.001 );
  TH1I hda1q( "da1q",
	      "da threshold Landau peak;#Deltax [#mum];tracks",
	      501, -pitch*1.001, pitch*1.001 );
  TProfile mada1qvsxm( "mada1qvsxm",
		       "MAD(#Deltaa) vs x;track x [#mum];MAD(#Deltaa) [#mum]",
		       pitch, -0.5*pitch, 0.5*pitch );

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // silicon:

  double ZA = 14.0; // ZA = atomic number of absorber, Si
  double AW = 28.086; // AW = atomic weight of absorber
  double rho = 2.329; // rho= density of absorber material
  double radl = 9.36; // [cm]

  double atnu = 6.0221367e23 * rho / AW; // atnu = # of atoms per cm**3

  if(seed != 0) {
      std::cout << "SEEDING with " << seed << std::endl;
      rgen.seed(seed); // seconds since 1.1.1970
  } else {
      rgen.seed( time(NULL) ); // seconds since 1.1.1970
  }

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // INITIALIZE ENERGY BINS

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

  cout
    << endl
    << "  n2 " << n2 << ", Emin " << Emin << ", um " << um
    << ", E[nume] " << E[nume] << endl;

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // READ DIELECTRIC CONSTANTS

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

    cout
      << endl
      << "  HEPS.TAB: n2t " << n2t << ", numt " << numt << endl;

    if( n2 != n2t ) cout << " CAUTION: n2 & n2t differ" << endl;
    if( nume != numt ) cout << " CAUTION: nume & numt differ" << endl;
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

    // MAZZIOTTA: 0.0 at 864
    // EP( 2, 864 ) = 0.5 * ( EP(2, 863) + EP(2, 865) )
    // RIM(864) = 0.5 * ( RIM(863) + RIM(865) )
    // DFDE(864) = RIM(864) * 0.0092456 * E(864)
    // DP: fixed in HEPS.TAB

  } // HEPS.TAB

  //= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  // READ INTEGRAL OVER MOMENTUM TRANSFER OF THE GENERALIZED OSCILLATOR STRENGTH

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

    cout
      << endl
      << "  MACOM.TAB: n2t " << n2t << ", numt " << numt << endl;

    if( n2 != n2t ) cout << " CAUTION: n2 & n2t differ" << endl;
    if( nume != numt ) cout << " CAUTION: nume & numt differ" << endl;
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

    istringstream tokenizer( line );

  // EMERC.TAB is the table of the integral over K of generalized oscillator
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
  // SHELL INITIALIZATION

  nvac[1] = 0;
  nvac[2] = 2;
  nvac[3] = 2;
  nvac[4] = 9; // possible transitions to it

  Eshel[1] =   12.0; // valence band upper edge (holes live below)
  Eshel[2] =   99.2; // M
  Eshel[3] =  148.7; // L
  Eshel[4] = 1839.0; // K

  for( unsigned n = 1; n <= 4; ++n )
    for( unsigned i = 1; i <= 9; ++i ) {
      augmi[n][i] = 0;
      augde[n][i] = 0;
    }

  // augmi(KSH,J) = PROBABILITA" INTEGRALI DEI VARI PROCESSI DI
  // EMISSIONE AUGR DALLA SHELL KSH

  // augmi(KSH,J) = 1 PER L"ULTIMO VALORE DI J

  // KSH = 4 --> SHELL K
  // KSH = 3 --> SHELL L1
  // KSH = 2 --> SHELL L23

  augmi[4][1] = 0.1920;
  augmi[4][2] = 0.3885 + augmi[4][1];
  augmi[4][3] = 0.2325 + augmi[4][2];
  augmi[4][4] = 0.0720 + augmi[4][3];
  augmi[4][5] = 0.0030 + augmi[4][4];
  augmi[4][6] = 0.1000 + augmi[4][5];
  augmi[4][7] = 0.0040 + augmi[4][6];
  augmi[4][8] = 0.0070 + augmi[4][7];
  augmi[4][9] = 0.0010 + augmi[4][8];
  augmi[3][1] = 0.0250;
  augmi[3][2] = 0.9750 + augmi[3][1];
  augmi[2][1] = 0.9990;
  augmi[2][2] = 0.0010 + augmi[2][1];

  // augde[KSH, J) = ENERGIA IN eV DELL"ELETTRONE AUGR
  // EMESSO DALLA SHELL KSH NEL PROCESSO J

  augde[4][1] = 1541.6;
  augde[4][2] = 1591.1;
  augde[4][3] = 1640.6;
  augde[4][4] = 1690.3;
  augde[4][5] = 1690.3;
  augde[4][6] = 1739.8;
  augde[4][7] = 1739.8;
  augde[4][8] = 1839.0;
  augde[4][9] = 1839.0;

  augde[3][1] = 148.7;
  augde[3][2] =  49.5;

  augde[2][1] = 99.2;
  augde[2][2] =  0.0;

  // EPP(I) = VALORI DI ENERGIA PER TABULARE LE PROBABILITA" DI
  // FOTOASSORBIMENTO NELLE VARIE SHELL
  // PM, PL23, PL1, PK = PROBABILITA" DI ASSORBIMENTO DA PARTE
  // DELLE SHELL M,L23,L1 E K

  // VALORI ESTRAPOLATI DA FRASER

  EPP[1]  = 40.0;
  PM[1]   = 1.0;
  PL23[1] = 0.0;
  PL1[1]  = 0.0;
  PK[1]   = 0.0;

  EPP[2]  = 50.0;
  PM[2]   = 1.0;
  PL23[2] = 0.0;
  PL1[2]  = 0.0;
  PK[2]   = 0.0;

  EPP[3]  = 99.2;
  PM[3]   = 1.0;
  PL23[3] = 0.0;
  PL1[3]  = 0.0;
  PK[3]   = 0.0;

  EPP[4]  = 99.2;
  PM[4]   = 0.03;
  PL23[4] = 0.97;
  PL1[4]  = 0.0;
  PK[4]   = 0.0;

  EPP[5]  = 148.7;
  PM[5]   = 0.03;
  PL23[5] = 0.92;
  PL1[5]  = 0.0;
  PK[5]   = 0.0;

  EPP[6]  = 148.7;
  PM[6]   = 0.02;
  PL23[6] = 0.88;
  PL1[6]  = 0.1;
  PK[6]   = 0.0;

  EPP[7]  = 150.0;
  PM[7]   = 0.02;
  PL23[7] = 0.88;
  PL1[7]  = 0.1;
  PK[7]   = 0.0;

  EPP[8]  = 300.0;
  PM[8]   = 0.02;
  PL23[8] = 0.83;
  PL1[8]  = 0.15;
  PK[8]   = 0.0;

  EPP[9]  = 500.0;
  PM[9]   = 0.02;
  PL23[9] = 0.70;
  PL1[9]  = 0.28;
  PK[9]   = 0.0;

  EPP[10]  = 1000.0;
  PM[10]   = 0.03;
  PL23[10] = 0.55;
  PL1[10]  = 0.42;
  PK[10]   = 0.0;

  EPP[11]  = 1839.0;
  PM[11]   = 0.05;
  PL23[11] = 0.39;
  PL1[11]  = 0.56;
  PK[11]   = 0.0;

  EPP[12]  = 1839.0;
  PM[12]   = 0.0;
  PL23[12] = 0.0;
  PL1[12]  = 0.08;
  PK[12]   = 0.92;

  EPP[13]  = 2000.0;
  PM[13]   = 0.0;
  PL23[13] = 0.0;
  PL1[13]  = 0.08;
  PK[13]   = 0.92;

  // EGAP = GAP ENERGY IN eV
  // EMIN = THRESHOLD ENERGY (ALIG ET AL., PRB22 (1980), 5565)

  double Egap = 1.17 - 4.73e-4 * temp*temp / (636+temp);

  double Ethr = 1.5*Egap; // energy conservation

  double eom0 = 0.063; // phonons
  double aaa = 5.2;    // Alig 1980

  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  // EVENT LOOP:

  for( unsigned iev = 0; iev < nev; ++ iev ) {

    cout << iev << endl;

    stack <delta> deltas;

    // put track on stack:

    double xm = pitch * ( unirnd(rgen) - 0.5 ); // [mu] -p/2..p/2 at track mid

    delta t;
    t.E = Ekin0; // [MeV]
    t.x = ( xm - 0.5*width ) * 1e-4; // entry point is left;
    t.y = 0; // [cm]
    t.z = 0; // pixel from 0 to depth [cm]
    t.u = sin(turn);
    t.v = 0;
    t.w = cos(turn); // along z
    t.npm = npm0;
    deltas.push(t); // beam particle is first "delta"

    unsigned it = 0;
    unsigned nscat = 0; // elastic
    unsigned nloss = 0; // ionization
    unsigned ndelta = 0;
    double tde = 0.0;
    unsigned meh = 0;
    unsigned sumeh2 = 0;
    vector <cluster> clusters;
    double Ekprev = 9e9; // update flag

    while( ! deltas.empty() ) {

      delta t = deltas.top();
      deltas.pop();

      double Ek = t.E; // [MeV] kinetic energy
      double xx = t.x;
      double yy = t.y;
      double zz = t.z;
      double vect[3];
      vect[0] = t.u; // direction cosines
      vect[1] = t.v;
      vect[2] = t.w;
      unsigned npm = t.npm;
      unsigned nlast = nume;

      double xm0 = 1;
      double xlel = 1;
      double gn = 1;
      double totsig[lime];
      double ptm = elmm; // e 0.51100 MeV

      cout << "  delta " << Ek*1e3 << " keV"
	   << ", cost " << t.w
	   << ", u " << t.u
	   << ", v " << t.v
	   << ", z " << zz*1e4;

      while(1) { // steps

	if( Ek < 0.9 * Ekprev ) { // update

	  double zi = 1.0;
	  ptm = elmm; // e 0.51100 MeV
	  if(      npm == 1 ) ptm = 938.2723; // proton
	  else if( npm == 2 ) ptm = 139.578; // pion
	  else if( npm == 3 ) ptm = 493.67; // K
	  else if( npm == 5 ) ptm = 105.65932; // mu

	  double gam = Ek / ptm + 1.0; // W = total energy / restmass
	  double bg  = sqrt( gam*gam - 1.0 ); // bg = beta*gamma = p/m
	  double pmom = ptm*bg; // [MeV/c]
	  double betasq = bg*bg / ( 1 + bg*bg );
	  double Emax = ptm * ( gam*gam - 1 ) / ( 0.5*ptm/elmm + 0.5*elmm/ptm + gam );
	  // Emax=maximum energy loss, see Uehling, also Sternheimer & Peierls Eq.(53)
	  if( npm == 4 ) Emax = 0.5*Ek;
	  // maximum energy loss for incident electrons
	  Emax = 1e6 * Emax; // eV

	  // Define parameters and calculate Inokuti"s sums,
	  // S ect 3.3 in Rev Mod Phys 43, 297 (1971)

	  double dec = zi*zi * atnu * fac / betasq;
	  double bemx = betasq / Emax;
	  double EkeV = Ek * 1e6; // [eV]
	  double twombb = 2 * elm * betasq; // [eV]

	  // Generate collision spectrum sigma(E) from df/dE, epsilon and AE.
	  // sig(*,j) actually is E**2 * sigma(E)

	  double tsig[6];

	  for( unsigned i = 1; i <= 5; ++i )
	    tsig[i] = 0;

	  double stpw = 0;
	  double H[lime];

	  for( unsigned j = 1; j <= nume; ++j ) {

	    if( E[j] > Emax ) break;

	    double Q1 = Ry;

	    // Eq. (3.1) in RMP and red notebook CCS-33, 39 & 47

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

	    sig[2][j] = sgg;
	    sig[3][j] = sgh; // small, negative

	    //sig[2][j] = 0; // TEST, worse resolution: more fluctuations
	    //sig[2][j] *= 2; // TEST, better resolution: less fluctuations
	    //if( E[j] > 1838 ) sig[2][j] = 0; // TEST, 7% better resolution

	    double uef = 1 - E[j] * bemx;
	    if( npm == 4 )
	      uef = 1 +
		pow( E[j] / ( EkeV - E[j] ), 2 ) +
		pow( (gam-1) / gam * E[j]/EkeV, 2 ) -
		( 2*gam-1 ) * E[j] / ( gam*gam * ( EkeV - E[j] ) );

	    // uef from  Eqs. 9 & 2 in Uehling, Ann Rev Nucl Sci 4, 315 (1954)
	    // if( j == 1) PRINT*, " uef=",UEF

	    sig[4][j] = 2 * sig[6][j] * uef;

	    // there is a factor of 2 because the integral was over d(lnK)
	    // rather than d(lnQ)

	    sig[5][j] = 0;

	    for( unsigned i = 1; i <= 4; ++i ) {

	      sig[5][j] += sig[i][j]; // sum

	      // divide by E**2 to get the differential collision cross section sigma
	      // Tsig = integrated total collision cross section

	      tsig[i]  = tsig[i]  + sig[i][j] * dE[j] / ( E[j]*E[j] );

	    } // i

	    tsig[5] += sig[5][j] * dE[j] / ( E[j]*E[j] ); // running sum

	    double HE2  = sig[5][j] * dec;
	    H[j] = HE2 / ( E[j]*E[j] );
	    stpw += H[j] * E[j] * dE[j]; // dE/dx
	    /*
	    if( ndelta == 0 )
	      cout << j
		   << "  " << E[j]
		   << "  " << dfdE[j]
		   << "  " << sgg
		   << "  " << sgh
		   << "  " << sig[1][j]
		   << "  " << sig[3][j]
		   << "  " << sig[4][j]
		   << "  " << sig[5][j]
		   << "  " << HE2 << endl; // compare Bichsel CONV-5000.OPA: agree
	    */
	    nlast = j;

	  } // j

	  xm0 = tsig[5] * dec; // 1/path

	  totsig[1] = H[1]*dE[1]; // running integral
	  double sst = H[1]*dE[1]; // total cross section (integral)
	  for( unsigned j = 2; j <= nlast; ++j ) {
	    totsig[j] = totsig[j-1] + H[j]*dE[j];
	    sst += H[j]*dE[j];
	  }

	  // NORMALIZE running integral:

	  for( unsigned j = 1; j <= nlast; ++j )
	    totsig[j] /= totsig[nlast]; // norm

	  // elastic:

	  if( npm == 4 ) { // ELECTRONS

	    //gn = 2*2.61 * pow( ZA, 2.0/3.0 ) / EkeV; // Mazziotta
	    gn = 2*2.61 * pow( ZA, 2.0/3.0 ) / (pmom*pmom)*1e-6; // Moliere
	    double E2 = 14.4e-14; // [MeV*cm]
	    double FF = 0.5*pi * E2*E2 * ZA*ZA / (Ek*Ek);
	    double S0EL = 2*FF / ( gn * ( 2 + gn ) );
	    // elastic total cross section  [cm2/atom]
	    xlel = atnu*S0EL; // ATNU = N_A * rho / A = atoms/cm3

	  }
	  else { //  OTHER PARTICLES

	    double getot = Ek + ptm;
	    xlel = min( 2232.0 * radl * pow( pmom*pmom / (getot*zi), 2 ), 10.0*radl );
	    // units ?
	  }

	  elvse.Fill( log(Ek)/log10, 1e4/xlel );
	  invse.Fill( log(Ek)/log10, 1e4/xm0 );

	  Ekprev = Ek;

	  if( ldb )
	    cout << "  ev " << iev << " type " << npm << ", Ekin " << Ek*1e3 << " keV"
		 << ", beta " << sqrt(betasq) << ", gam " << gam << endl
		 << "  Emax " << Emax << ", nlast " << nlast << ", Elast " << E[nlast]
		 << ", norm " << totsig[nlast] << endl
		 << "  inelastic " << 1e4/xm0 << "  " << 1e4/sst
		 << ", elastic " << 1e4/xlel << " um"
		 << ", mean dE " << stpw*depth*1e-4*1e-3 << " keV"
		 << endl << flush;

	} // update

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// step:

	double tlam = 1 / ( xm0 + xlel ); // [cm] TOTAL MEAN FREE PATH (MFP)

	double xr = -log( 1 - unirnd(rgen) ) * tlam; // exponential step length

	hstep5.Fill( xr*1e4 );
	hstep0.Fill( xr*1e4 );

	zz += xr*vect[2];

	if( ldb && Ek < 1 )
	  cout << "step " << xr*1e4 << ", z " << zz*1e4 << endl;

	hzz.Fill( zz*1e4 );

	if( zz < 0 || zz > depth*1e-4 ) break; // exit back or front

	xx += xr*vect[0];
	yy += xr*vect[1];

	if( fabs( yy ) > 0.0200 ) break; // save time

	++it;

	if( unirnd(rgen) > tlam*xlel ) { // INELASTIC (ionization) PROCESS

	  ++nloss;

	  // GENERATE VIRTUAL GAMMA:

	  double yr = unirnd(rgen); // inversion method
	  unsigned je = 2;
	  for( ; je <= nlast; ++je )
	    if( yr < totsig[je] ) break;

	  double Eg = E[je-1] + ( E[je] - E[je-1] ) * unirnd(rgen); // [eV]

	  hde0.Fill( Eg ); // M and L shells
	  hde1.Fill( Eg ); // K shell
	  hde2.Fill( Eg*1e-3 );
	  hdel.Fill( log(Eg)/log10 );

	  double resekin = Ek - Eg*1E-6; // [ MeV]

	  // cut off for further movement: [MeV]

	  if( resekin < explicit_delta_energy_cut_keV*1e-3 ) {

	    // cout << "@@@ NEG RESIDUAL ENERGY" << Ek*1e3 << Eg*1e-3 << resekin*1e-3
	    Eg = Ek*1E6; // [eV]
	    resekin = Ek - Eg; // zero
	    // cout << "LAST ENERGY LOSS" << Eg << resekin

	  }

	  //if( Eg < explicit_delta_energy_cut_keV*1e3 ) // avoid double counting
	  tde += Eg; // [eV]

	  // emission angle from delta:

	  // PRIMARY SCATTERING ANGLE:
	  // SINT = SQRT(ER*1.E-6/EK) ! Mazziotta = deflection angle
          // COST = SQRT(1.-SINT*SINT)
	  // STORE INFORMATION ABOUT DELTA-RAY:
	  // SINT = COST ! flip
	  // COST = SQRT(1.-SINT**2) ! sqrt( 1 - ER*1e-6 / Ek ) ! wrong

	  //double cost = sqrt( Eg / (twome + Eg) ); // M. Swartz
	  double cost = sqrt( Eg / (twome + Eg) * ( Ek + twome*1e-6 ) / Ek );
	  // Penelope, Geant4
	  double sint;
	  if( cost*cost <= 1 )
	    sint = sqrt( 1 - cost*cost ); // mostly 90 deg
	  else {
	    cout << " NAN 1-cost " << 1-cost << ", 1-cost^2 " << 1-cost*cost << endl;
	    sint = 0;
	  }
	  double phi = 2*pi*unirnd(rgen);

	  // G4PenelopeIonisationModel.cc

	  // rb = kineticEnergy + 2*electron_mass_c2;

	  // kineticEnergy1 = kineticEnergy - deltaE;
	  // cosThetaPrimary = sqrt( kineticEnergy1 * rb /
	  // ( kineticEnergy * ( rb - deltaE ) ) );

	  // cosThetaSecondary = sqrt( deltaE * rb /
	  // ( kineticEnergy * ( deltaE + 2*electron_mass_c2 ) ) );

	  // penelope.f90
	  // Energy and scattering angle ( primary electron ).
	  // EP = E - DE
	  // TME = 2 * ME
	  // RB = E + TME
	  // CDT = SQRT( EP * RB / ( E * ( RB - DE ) ) )

	  // emission angle of the delta ray:
	  // CDTS = SQRT( DE * RB / ( E * ( DE + TME ) ) ) // like Geant

	  vector <double> din(3);
	  din[0] = sint*cos(phi);
	  din[1] = sint*sin(phi);
	  din[2] = cost;

	  htet.Fill( wt*asin(sint) ); // peak at 90, tail to 45, elastic forward

	  // transform into detector system:

	  double cz = vect[2];     // delta direction
	  double sz = sqrt( 1 - cz*cz );
	  double phif = atan2( vect[1], vect[0] );
	  double sf = sin(phif);
	  double cf = cos(phif);
	  double uu =  cz*cf * din[0] - sf * din[1] + sz*cf * din[2];
	  double vv =  cz*sf * din[0] + cf * din[1] + sz*sf * din[2];
	  double ww = -sz    * din[0]               + cz    * din[2];

	  // GENERATE PRIMARY e-h:

	  stack <double> veh;

	  if( Eg > Ethr )
	    shells( Eg, veh );

	  hnprim.Fill( veh.size() );

	  // process e and h

	  double sumEeh{0};
	  unsigned neh{0};

	  while( ! veh.empty() ) {

	    double Eeh = veh.top();

	    if( Eeh > 1 )
	      hlogE.Fill( log(Eeh)/log10 );
	    else
	      hlogE.Fill( 0 );

	    //cout << "    eh "<< veh.size() << ", Eeh " << Eeh << ", neh " << neh << endl;

	    veh.pop();

	    if( Eeh > explicit_delta_energy_cut_keV*1e3 ) {

	      // put delta on stack:

	      delta t;
	      t.E = Eeh*1E-6; // Ekin [MeV]
	      t.x = xx;
	      t.y = yy;
	      t.z = zz;
	      t.u = uu;
	      t.v = vv;
	      t.w = ww;
	      t.npm = 4; // e
	      deltas.push(t);

	      ++ndelta;

	      tde -= Eeh; // [eV], avoid double counting

	      continue; // next ieh

	    } // new delta

	    sumEeh += Eeh;

	    // slow down low energy e and h: 95% of CPU time

	    while( fast == 0 && Eeh > Ethr ) {

	      double pion = 1 / ( 1 + aaa*105/twopi * sqrt(Eeh-eom0) /
				  pow( Eeh - Ethr, 3.5 ) );
	      // for e and h

	      if( unirnd(rgen) < pion ) { // ionization

		++neh;

		double E1 = gena1() * (Eeh-Ethr);
		double E2 = gena2() * (Eeh-Ethr-E1);

		//cout << "      ion " << Eeh << " => " << E1 << " + " << E2 << endl;

		if( E1 > Ethr )
		  veh.push( E1 );
		if( E2 > Ethr )
		  veh.push( E2 );

		Eeh = Eeh - E1 - E2 - Ethr;
	      }
	      else
		Eeh = Eeh - eom0; // phonon emission
	      // cout << "      fon " << ed

	    } // slow: while Eeh

	  } // while veh

	  if( fast ) {
	    poisson_distribution <int> poisson(sumEeh/3.645);
	    neh = poisson(rgen);
	  }

	  meh += neh;
	  sumeh2 += neh*neh;

	  //cout << "  dE " << Eg << " eV, neh " << neh << endl;

	  // store charge cluster:

	  if( neh > 0 ) {

	    hlogn.Fill( log(neh)/log10 );

	    cluster c;

	    c.neh = neh;
	    c.x = xx;
	    c.y = yy;
	    c.z = zz;
	    c.E = Eg; // [eV]
	    clusters.push_back(c);

	  } // neh

	  Ek -= Eg*1E-6; // [MeV]

	  if( ldb && Ek < 1 )
	    cout << "    Ek " << Ek*1e3
		 << " keV, z " << zz*1e4 << ", neh " << neh
		 << ", steps " << it << ", ion " << nloss << ", elas " << nscat
		 << ", cl " << clusters.size()
		 << endl;

	  if( Ek < 1E-6 || resekin < 1E-6 ) {
	    // cout << "  absorbed" << endl;
	    break;
	  }

	  if( npm == 4 ) { // electrons, update elastic cross section at new Ek

	    //gn = 2*2.61 * pow( ZA, 2.0/3.0 ) / (Ek*1E6); // Mazziotta
	    double pmom = sqrt( Ek * ( Ek + 2*ptm ) ); // [MeV/c] 2nd binomial
	    gn = 2*2.61 * pow( ZA, 2.0/3.0 ) / (pmom*pmom)*1e-6; // Moliere
	    double E2 = 14.4e-14; // [MeV*cm]
	    double FF = 0.5*pi * E2*E2 * ZA*ZA / (Ek*Ek);
	    double S0EL = 2*FF / ( gn * ( 2 + gn ) );
	    // elastic total cross section  [cm2/atom]
	    xlel = atnu*S0EL; // ATNU = N_A * rho / A = atoms/cm3

	  }

	}
	else { // ELASTIC SCATTERING: Chaoui 2006

	  ++nscat;

	  double r = unirnd(rgen);
	  double cost = 1 - 2*gn*r / ( 2 + gn - 2*r );
	  double sint = sqrt( 1 - cost*cost );

	  double phi = twopi*unirnd(rgen);

	  vector <double> din(3);
	  din[0] = sint*cos(phi);
	  din[1] = sint*sin(phi);
	  din[2] = cost;

	  hscat.Fill( wt*asin(sint) ); // forward peak, tail to 90

	  // change direction of delta VECT:

	  double cz = vect[2];       // delta direction
	  double sz = sqrt( 1.0 - cz*cz );
	  double phif = atan2( vect[1], vect[0] );
	  double sf = sin(phif);
	  double cf = cos(phif);
	  vect[0] =  cz*cf * din[0] - sf * din[1] + sz*cf * din[2];
	  vect[1] =  cz*sf * din[0] + cf * din[1] + sz*sf * din[2];
	  vect[2] = -sz    * din[0]               + cz    * din[2];

	} // elastic

      } // while steps

      cout << endl;

      Ekprev = 9e9; // update flag for next delta

    } // while deltas
    /*
    write( 69, * ) "ev " << iev << ncl << meh // event header

      do i = 1, ncl
	   write( 69, "( 3F7.1, x, f9.1, I6 )" )
	   vclu(i,1)*1e4 << vclu(i,2)*1e4 << vclu(i,3)*1e4 << // [mu]
	   vclu(i,4) << // [eV]
	   kclu(i)
	   enddo
    */

    cout
      << "  steps " << it << ", ion " << nloss << ", elas " << nscat
      << ", dE " << tde*1e-3 << " keV"
      << ", eh " << meh
      << ", cl " << clusters.size()
      << endl;

    hncl.Fill( clusters.size() );
    htde.Fill( tde*1e-3 ); // [keV] energy conservation - binding energy
    if( ndelta )
      htde1.Fill( tde*1e-3 ); // [keV]
    else
      htde0.Fill( tde*1e-3 ); // [keV]
    hteh.Fill( meh*1e-3 ); // [ke]
    hq0.Fill( meh*1e-3 ); // [ke]
    hrms.Fill( sqrt(sumeh2) );

    // 4 pixels along x:

    double q1[4];
    for( int ir = 0; ir < 4; ++ir ) {
      q1[ir] = 0;
    }

    for( unsigned i = 0; i < clusters.size(); ++i ) {

      double xx = clusters[i].x*1e4; // [mu]
      double yy = clusters[i].y*1e4; // [mu]
      double zz = clusters[i].z*1e4; // [mu]
      int neh = clusters[i].neh;

      if( iev < 11 ) {
	h1zev[iev]->Fill( zz, neh );
	h2zxev[iev]->Fill( xx, zz, neh );
      }
      h2xy->Fill( xx, yy, neh );
      h2zx->Fill( xx, zz, neh );

      double Eg = clusters[i].E;
      hcleh.Fill( log(neh)/log10 ); // per cluster
      hreh.Fill( neh/Eg );
      wvse.Fill( log(Eg)/log10, Eg/neh ); // dE per eh pair

      // 0 | 1 | 2 | 3, bins 0 and 3 are half-infinite
      // diffusion across x crack: |  |  |

      double xc = -pitch; // left
      int m = 0; // minus = left pixel
      int p = 1; // plus = right pixel
      if( xx > 0.5*pitch ) { // nearer x-crack is right
	xc = pitch;
	m = 2;
	p = 3;
      }
      else if( xx > -0.5*pitch ) { // mid crack
	xc = 0;
	m = 1;
	p = 2;
      }

      double dtime = zz*1e-4/vd; // [s] drift time along z (mean speed theorem)
      double diff = sqrt(2*D*dtime)*1e4; // [mu] rms diffusion (projected or 3D?)
      double uu = -(xx-xc)/w2/diff; // scaled diffusion distance for erfc
      double tf = 0.5*erfc(uu); // upper Gaussian tail fraction

      hdtime.Fill( dtime*1e9 );
      hdiff.Fill( diff );
      htf.Fill( tf );
      tfvsx.Fill( xx, tf ); // S-curve, x = 0 is a pixel boundary

      q1[p] += neh*tf;
      q1[m] += neh*(1-tf);

    } // clusters

    double q0 = q1[0]+q1[1]+q1[2]+q1[3];
    double eta = (q1[2]-q1[1])/(q1[2]+q1[1]); // central bins
    heta0.Fill( eta );
    eta0vsxm.Fill( xm, eta );

    double sumq0 = 0;
    double sumqx0 = 0;
    for( int ir = 0; ir < 4; ++ir ) {
      sumq0 += q1[ir];
      sumqx0 += q1[ir]*(ir-1.5); // -1.5, -0.5, 0.5, 1.5
    }
    double cog0 = sumqx0/sumq0;
    double dx0 = cog0*pitch - xm; // [mu]
    hdx0.Fill( dx0 );

    madx0vsq.Fill( q0*1e-3, fabs(dx0) ); // linear rise, too steep

    if( q0 < 95*depth ) { // keep 2/3 in 300 mu
      hdx0q.Fill( dx0 );
      dx0qvsxm.Fill( xm, dx0 );
      madx0qvsxm.Fill( xm, fabs(dx0) ); // inverted U-shape
      if( ndelta )
	madx0qvsxm1.Fill( xm, fabs(dx0) ); // inverted U-shape
      else
	madx0qvsxm0.Fill( xm, fabs(dx0) ); // inverted U-shape
    }

    // after threshold:

    int npx = 0;
    double sumq1 = 0;
    double sumqx1 = 0;
    for( int ir = 0; ir < 4; ++ir ) {
      if( q1[ir] > thr ) {
	++npx;
	hpxq1.Fill( q1[ir]*1e-3 ); // [ke]
	sumq1 += q1[ir];
	sumqx1 += q1[ir]*(ir-1.5); // -1.5, -0.5, 0.5, 1.5
      }
    }
    hq1.Fill( sumq1*1e-3 ); // [ke]
    hnpx1.Fill( npx );
    npx1vsxm.Fill( xm, npx );

    double eta1 = (q1[2]-q1[1])/(q1[2]+q1[1]); // central bins
    heta1.Fill( eta1 );
    eta1vsxm.Fill( xm, eta1 );

    double cog1 = sumqx1/sumq1;
    x1vsxm.Fill( xm, cog1*pitch );
    double dx1 = cog1*pitch - xm;
    hdx1.Fill( dx1 );

    madx1vsq.Fill( sumq1*1e-3, fabs(dx1) ); // linear rise, too steep

    if( sumq1 < 95*depth ) { // keep 2/3 in 300 mu
      hdx1q.Fill( dx1 );
      // pitch  thck thr  sigma
      // 25 mu   50  500  2.48 mu
      // 25 mu  100  500  1.69 mu
      // 25 mu  150  500  1.4  mu
      // 25 mu  285  500  1.08 mu
      // 25 mu  285  700  1.13 mu
      // 25 mu  450  500  0.93 mu
      // 17 mu  285  500  0.73 mu
      // 17 mu  285  700  0.77 mu
      // 10 mu  285  500  0.43 mu
      dx1qvsxm.Fill( xm, dx1 );
      madx1qvsxm.Fill( xm, fabs(dx1) ); // inverted U-shape
    }
    // cross talk:
    /*
    double am = qm*(1-cx) + qp*cx;
    double ap = qp*(1-cx) + qm*cx;
    cog = ( -0.5*pitch*am + 0.5*pitch*ap )/(ap+am); // = 0.5*pitch*eta
    double dat1 = cog - xm;
    hda1.Fill( dat1 );

    if( q0 < 90*depth ) { // keep 2/3
      hda1q.Fill( dat1 );
      da1qvsxm.Fill( xm, dat1 );
      mada1qvsxm.Fill( xm, fabs(dat1) ); // with cross talk: flat, bud bad
    }
    */

  } // events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  histoFile->Write();
  histoFile->ls();
  histoFile->Close();
  cout << endl;

  cout << "done: events " << nev << endl;

  cout << "  particle type     " << npm0 << endl;
  cout << "  kinetic energy    " << Ekin0 << " MeV" << endl;
  cout << "  number of events  " << nev << endl;
  cout << "  pixel pitch       " << pitch << " um" << endl;
  cout << "  thickness         " << depth << " um" << endl;
  cout << "  incident angle    " << turn*wt << " deg" << endl;
  cout << "  track width       " << width << " um" << endl;
  cout << "  temperature       " << temp << " K" << endl;
  cout << "  readout threshold " << thr << " e" << endl;
  cout << "  cross talk        " << cx*100 << "%" << endl;

  cout
    << endl
    << ( (j) ? "  holes" : "  electrons" ) << endl
    << "  mobility for " << Efield << " V/cm"
    << ": vm " << vm // cm/s = 100 um / ns
    << ", Ec " << Ec
    << ", mu0 " << mu0 << endl
    << "  beta " << beta
    << ", mu " << mu
    << ", v " << vd << " cm/s"
    << " = " << vd/1e5 << " mu/ns" << endl
    << "  D " << D
    << ", rms " << sqrt(2*D*4e-9)*1e4 << " mu" // for 4 ns drift
    << endl;

  cout << endl
       << "  " << histoFile->GetName() << endl;
  cout << endl;

} // main

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double alph1( double x ) // x = 0..1
{
  return 105./16. * (1.-x)*(1-x) * sqrt(x); // integral = 1, max = 1.8782971
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double alph2( double x ) // x = 0..1
{
  static const double pi = 3.1415926536;
  return 8/pi * sqrt( x*(1-x) );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double trian( double x ) // x = -1..1
{
  if( x < 0 )
    return x + 1;
  else
    return -x + 1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double gena1()
{
 l31:
  double r1 = unirnd(rgen);
  double r2 = unirnd(rgen);
  if( alph1( r1 ) > 1.8783*r2 ) goto l31; // rejection method

  return r1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double gena2()
{
 l32:
  double r1 = unirnd(rgen);
  double r2 = unirnd(rgen);
  if( alph2( r1 ) > 1.27324*r2 ) goto l32; // rejection method

 return r1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double gentri() // -1..1
{
 l33:
  double r1 = unirnd(rgen);
  double r2 = unirnd(rgen);
  double x = -1 + 2*r1;
  if( trian(x) > r2 ) goto l33; // rejection method

  return x;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void shells( double Eg, stack <double> &veh )
{
  // INPUT:
  // EG = VIRTUAL GAMMA ENERGY [eV]

  // OUTPUT:
  // veh ENERGIES OF PRIMARY e/h

  //Eg = Eg - Egap; // double counting?

  // EV = binding ENERGY OF THE TOP OF THE VALENCE BAND
  double Ev = Eshel[1]; // 12.0 eV

  double PV[5];

  int is = -1;
  if( Eg <= Eshel[1] )
    is = 0;
  else if( Eg <= EPP[3] )
    is = 1;
  else {
    if( Eg > EPP[nep] ) {
      PV[1] = PM[nep];
      PV[2] = PL23[nep];
      PV[3] = PL1[nep];
      PV[4] = PK[nep];
    }
    else {
      unsigned iep = 3;
      for( ; iep < nep; ++iep )
	if( Eg > EPP[iep] && Eg <= EPP[iep+1] ) break;

      // interpolate:

      PV[1] = PM[iep] + ( PM[iep+1] - PM[iep] ) / ( EPP[iep+1] - EPP[iep] ) * ( Eg - EPP[iep] );
      PV[2] = PL23[iep] + ( PL23[iep+1] - PL23[iep] ) / ( EPP[iep+1] - EPP[iep] ) * ( Eg - EPP[iep] );
      PV[3] = PL1[iep] + ( PL1[iep+1] - PL1[iep] ) / ( EPP[iep+1] - EPP[iep] ) * ( Eg - EPP[iep] );
      PV[4] = PK[iep] + ( PK[iep+1] - PK[iep] ) / ( EPP[iep+1] - EPP[iep] ) * ( Eg - EPP[iep] );
    }

    double PPV = PV[1] + PV[2] + PV[3] + PV[4];

    PV[2] = PV[1] + PV[2];
    PV[3] = PV[2] + PV[3];
    PV[4] = PV[3] + PV[4];

    double rs = unirnd(rgen);
    unsigned iv = 1;
    for( ; iv <= 4; ++iv ) {
      PV[iv] = PV[iv]/PPV; // normalize
      if( PV[iv] > rs ) break;
    }
    is = iv;
    if( is > 4 ) is = 4;
  }

  //cout << "  shells for " << Eg << " eV, Ev " << Ev << ", is " << is << endl;

  // PROCESSES:

  // PHOTOABSORPTION IN VALENCE BAND

  if( is <= 1 ) {

    double rv = unirnd(rgen);
    double Ee = Eg;
    if( Ee < 0.1 ) return;
    if( Ee < Ev ) {
      veh.push( rv*Ee );
      veh.push( (1-rv)*Ee );
    }
    else {
      veh.push( rv*Ev );
      veh.push( Ee - rv*Ev );
    }
    return;
  }

  // PHOTOABSORPTION IN AN INNER SHELL

  double Eth = Eshel[is];
  double Ephe = Eg - Eth;
  if( Ephe <= 0 ) {
    cout << "shells: photoelectron with negative energy "
	 << Eg << ", shell " << is << " at " << Eth << " eV" << endl;
    return;
  }

  // PRIMARY PHOTOELECTRON:

  veh.push( Ephe );

  // AUGR ELECTRONS:

  double raug = unirnd(rgen);

  int ks = 1;
  if( is <= 1 )
    ks = 1;
  else if( is <= 3 ) {
    ks = 1;
    if( raug > augmi[is][1] ) ks = 2;
  }
  else {
    if( raug < augmi[is][1] )
      ks = 1;
    else {
      for( int js = 2; js <= nvac[is]; ++js )
	if( raug >= augmi[is][js-1] && raug < augmi[is][js] )
	  ks = js;
    }
  }

  if( is == 2 ) {

    // L23-SHELL VACANCIES

     if( ks == 1 )
       TRL23MM( Ev, veh );

  }
  else if( is == 3 ) {

    // L1-SHELL VACANCIES

    if( ks == 2 ) {
      // TRANSITION L1 L23 M
      double rv = unirnd(rgen);
      veh.push( rv*Ev );
      veh.push( augde[is][ks] - rv*Ev );

      unsigned kks = 1;
      if( unirnd(rgen) > augmi[2][1] )
	kks = 2;
      if( kks == 1 )
	TRL23MM( Ev, veh );
    }
    else
      TRL1MM( Ev, veh );

  } // is 3

  if( is == 4 ) {

    // K-SHELL VACANCIES

    if( ks >= 8 ) {

      // TRANSITION K M M
      double r = gentri(); // -1..1
      // HF1( 112, r )
      double rEv = (1+r)*Ev; // 0..2*Ev
      veh.push( augde[is][ks] - rEv ); // adjust for energy conservation

    l654:
      double rv = unirnd(rgen);
      double Eh1 = rv * rEv;
      double Eh2 = (1-rv) * rEv;
      if( Eh1 > Ev || Eh2 > Ev ) goto l654; // holes stay below valence band edge (12 eV)

      veh.push( Eh1 );
      veh.push( Eh2 );

    }
    else if( ks == 6 || ks == 7 ) {

      // TRANSITION K L23 M

      double rEv = Ev*unirnd(rgen);
      veh.push( augde[is][ks] - rEv ); // adjust for energy conservation
      veh.push( rEv );

      unsigned kks = 1;
      if( unirnd(rgen) > augmi[2][1] ) kks = 2;
      if( kks == 1 )
	TRL23MM( Ev, veh );

    }
    else if( ks == 4 || ks == 5 ) {

      // TRANSITION K L1 M

      double rEv = Ev*unirnd(rgen);
      veh.push( augde[is][ks] - rEv ); // adjust for energy conservation
      veh.push( rEv );

      unsigned kks = 1;
      if( unirnd(rgen) > augmi[3][1] )
	kks = 2;

      if( kks == 1 )
	TRL1MM( Ev, veh );

      else {

	// TRANSITION L1 L23 M

	double rEv = Ev*unirnd(rgen);
	veh.push( rEv );
	veh.push( augde[3][kks] - rEv );

	unsigned kks = 1;
	if( unirnd(rgen) > augmi[2][1] )
	  kks = 2;
	if( kks == 1 )
	  TRL23MM( Ev, veh );

      }
    }
    else if( ks == 3 ) {

      // TRANSITION K L23 L23

      veh.push( augde[is][ks] ); // default

      unsigned kks = 1;
      if( unirnd(rgen) > augmi[2][1] ) kks = 2;
      if( kks == 1 )
	TRL23MM( Ev, veh );

      kks = 1;
      if( unirnd(rgen) > augmi[2][1] ) kks = 2;
      if( kks == 1 )
	TRL23MM( Ev, veh );

    }
    else if( ks == 2 ) {

      veh.push( augde[is][ks] ); // default

      // TRANSITION K L1 L23
      // L23-SHELL VACANCIES
      unsigned kks = 1;
      if( unirnd(rgen) > augmi[2][1] ) kks = 2;
      if( kks == 1 )
           TRL23MM( Ev, veh );

      // L1-SHELL VACANCIES
      kks = 1;
      if( unirnd(rgen) > augmi[3][1] ) kks = 2;
      if( kks == 2 ) {
	// TRANSITION L1 L23 M
	double rEv = Ev*unirnd(rgen);
	veh.push( rEv );
	veh.push( augde[3][kks] - rEv );

	kks = 1;
	if( unirnd(rgen) > augmi[2][1] ) kks = 2;
	if( kks == 1 )
	  TRL23MM( Ev, veh );
      }
      else
	TRL1MM( Ev, veh );
    }
    else if( ks == 1 ) {

      veh.push( augde[is][ks] ); // default
      // TRANSITION K L1 L1
      // L1-SHELL VACANCIES

      unsigned kks = 1;
      if( unirnd(rgen) > augmi[3][1] ) kks = 2;

      if( kks == 2 ) {
	// TRANSITION L1 L23 M
	double rEv = Ev*unirnd(rgen);
	veh.push( rEv );
	veh.push( augde[3][kks] - rEv );

	// L23-SHELL VACANCIES
	kks = 1;
	if( unirnd(rgen) > augmi[2][1] ) kks = 2;
	if( kks == 1 )
	  TRL23MM( Ev, veh );
      }
      else
	TRL1MM( Ev, veh );

      // L1-SHELL VACANCIES
      kks = 1;
      if( unirnd(rgen) > augmi[3][1] )
	kks = 2;

      if( kks == 2 ) {
	// TRANSITION L1 L23 M
	double rEv = Ev*unirnd(rgen);
	veh.push( rEv );
	veh.push( augde[3][kks] - rEv );

	// L23-SHELL
	kks = 1;
	if( unirnd(rgen) > augmi[2][1] ) kks = 2;
	if( kks == 1 )
	  TRL23MM( Ev, veh );
      }
      else
	TRL1MM( Ev, veh );

    } // ks

  } // is 4

} // SHELLS

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TRL1MM( double Ev, stack <double> &veh )
{
  // TRANSITION L1 M M

  // AUGER ELECTRON

  double rEv = ( 1 + gentri() ) * Ev; // 0..2*Ev

  veh.push( augde[3][1] - rEv );

  // ASSIGN ENERGIES TO THE HOLES

 l63:
  double rv = unirnd(rgen);
  double Eh1 = rv * rEv;
  double Eh2 = (1-rv) * rEv;
  if( Eh1 > Ev || Eh2 > Ev ) goto l63; // holes stay below valence band edge (12 eV)

  veh.push( Eh1 );
  veh.push( Eh2 );

} // TRL1MM

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TRL23MM( double Ev, stack <double> &veh )
{
  // TRANSITION L23 M M

  // AUGER ELECTRON

  double rEv = ( 1 + gentri() ) * Ev; // 0..2*Ev

  veh.push( augde[2][1] - rEv );

  // ASSIGN ENERGIES TO THE HOLES

 l62:
  double rv = unirnd(rgen);
  double Eh1 = rv * rEv;
  double Eh2 = (1-rv) * rEv;
  if( Eh1 > Ev || Eh2 > Ev ) goto l62; // holes stay below valence band edge (12 eV)

  veh.push( Eh1 );
  veh.push( Eh2 );

} // TRL23MM
