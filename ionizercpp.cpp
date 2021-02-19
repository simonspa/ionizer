
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

#include "DepositionBichsel.hpp"

#include <cmath>   // log
#include <cstdlib> // atoi
#include <ctime>
#include <iostream> // std::cout
#include <random>
#include <stack>

#include <Math/Vector3D.h>
#include <TFile.h> // ROOT
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

using namespace allpix;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char* argv[]) {
    // defaults:
    unsigned nev = 10 * 1000; // events
    double depth = DEPTH;     // [mu] pixel depth
    double pitch = 25;        // [mu] pixels size
    double angle = 999;       // flag
    double thr = 500;         // threshold [e]
    double cx = 0;            // cross talk
    double Ekin0 = EKIN;      // [MeV] kinetic energy

    uint64_t seed = 0;

    for(int i = 1; i < argc; ++i) {

        if(!strcmp(argv[i], "-p"))
            pitch = atof(argv[++i]); // [mu]

        if(!strcmp(argv[i], "-d"))
            depth = atof(argv[++i]); // [mu]

        if(!strcmp(argv[i], "-a"))
            angle = atof(argv[++i]); // [deg]

        if(!strcmp(argv[i], "-t"))
            thr = atof(argv[++i]); // pixel threshold [e]

        if(!strcmp(argv[i], "-c"))
            cx = atof(argv[++i]); // cross talk fraction

        if(!strcmp(argv[i], "-e"))
            Ekin0 = atof(argv[++i]); // [MeV]

        if(!strcmp(argv[i], "-n"))
            nev = atoi(argv[++i]);

        // if(!strcmp(argv[i], "-f"))
        // fast = 0; // full ionization, not fast: simulate each e-h pair

        if(!strcmp(argv[i], "-s"))
            seed = atoi(argv[++i]); // random seed

        if(!strcmp(argv[i], "-h")) {
            std::cout << "  simulate ionization by charged particles in silicon" << std::endl
                      << "  usage: ionizer [option] [option] [option]" << std::endl
                      << "         produces ionizer.hist" << std::endl
                      << "  options:" << std::endl
                      << "    -z pixel length [mu] (default 150)" << std::endl
                      << "    -p pixel width [mu] (default 25)" << std::endl
                      << "    -d pixel depth [mu] (default 285)" << std::endl
                      << "    -a angle of incidence [deg] (default atan(p/d))" << std::endl
                      << "    -t readout threshold [e] (default 500)" << std::endl
                      << "    -c cross talk fraction (default 0)" << std::endl
                      << "    -e incident kinetic energy [MeV] (default 5000)" << std::endl
                      << "    -f full ionization (slow, default neh=de/3.645)" << std::endl
                      << "    -n number of events (default 10000)" << std::endl;
            return 0;
        }
    } // argc

    double turn = atan(pitch / depth); // [rad] default
    if(fabs(angle) < 91)
        turn = angle / 180 * M_PI;

    double width = depth * tan(turn); // [mu] projected track, default: pitch

    // [V/cm] mean electric field: Vbias-Vdepletion/2
    double Efield = (120 - 30) / depth * 1e4; // UHH

    // from config:
    double temperature = TEMPERATURE;
    ParticleType default_particle_type = PARTICLE_TYPE;

    std::cout << "  particle type     " << default_particle_type << std::endl;
    std::cout << "  kinetic energy    " << Ekin0 << " MeV" << std::endl;
    std::cout << "  number of events  " << nev << std::endl;
    std::cout << "  pixel pitch       " << pitch << " um" << std::endl;
    std::cout << "  pixel depth       " << depth << " um" << std::endl;
    std::cout << "  incident angle    " << turn * 180 / M_PI << " deg" << std::endl;
    std::cout << "  track width       " << width << " um" << std::endl;
    std::cout << "  temperature       " << temperature << " K" << std::endl;
    std::cout << "  readout threshold " << thr << " e" << std::endl;
    std::cout << "  cross talk        " << cx * 100 << "%" << std::endl;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // mobility from pixelav:
    // 0 = e, 1 = h
    int j = 0;                        // e CMS, B2 pixel
    double cvm[2] = {1.53e9, 1.62e8}; // [cm/s] vmax at temperature=1K
    double evm[2] = {-0.87, -0.52};
    double vm = cvm[j] * pow(temperature, evm[j]);
    double cec[2] = {1.01, 1.24}; // [V/cm] Ecrit
    double eec[2] = {1.55, 1.68};
    double Ec = cec[j] * pow(temperature, eec[j]);
    double mu0 = vm / Ec;
    double cbeta[2] = {0.0257, 0.46};
    double ebeta[2] = {0.66, 0.17};
    double beta = cbeta[j] * pow(temperature, ebeta[j]);
    double ibeta = 1 / beta;
    double d2 = Efield / Ec;
    double d3 = pow(d2, beta) + 1.;
    double mu = mu0 / pow(d3, ibeta); // mu0 / ( 1 + (E/Ec)^b )^(1/b)
    const double vd = Efield * mu;    // [cm/s]
    // diffusion from mobility: D = kTmu/e
    // e = 1.602e-19 C
    // k = 1.38e-23 J/K
    // k/e = 8.6e-5 eV/K
    const double D = 8.61733e-5 * temperature * mu; // diffuson constant

    std::cout << std::endl
              << "   mobility for " << Efield << " V/cm"
              << ": vm " << vm // cm/s = 100 um / ns
              << ", Ec " << Ec << ", mu0 " << mu0 << std::endl
              << "  beta " << beta << ", mu " << mu << ", v " << vd << " cm/s"
              << " = " << vd / 1e5 << " mu/ns" << std::endl
              << "  D " << D << ", rms " << sqrt(2 * D * 4e-9) * 1e4 << " mu" // for 4 ns drift
              << std::endl;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // book histos

    TFile* histoFile = new TFile(Form("ionizercpp_p%i_w%i_d%i_t%i_c%i.hist",
                                      int(pitch + 0.5),
                                      int(width + 0.5),
                                      int(depth + 0.5),
                                      int(thr),
                                      int(100 * cx + 0.1)),
                                 "RECREATE");

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    TH1I* h1zev[11];
    TH2I* h2zxev[11];
    for(unsigned i = 0; i < 11; ++i) {
        h1zev[i] = new TH1I(Form("z%02i", i), Form("z event %i;z [#mum];clusters [eh-pairs]", i), 4 * depth, 0, depth);
        h2zxev[i] = new TH2I(Form("zx%02i", i),
                             Form("z-x event %i;x [#mum];z [#mum];clusters [eh-pairs]", i),
                             4 * 2 * pitch,
                             -pitch,
                             pitch,
                             4 * depth,
                             0,
                             depth);
    }

    TH2I* h2xy = new TH2I("xy", "x-y clusters;x [#mum];y [#mum];clusters [eh-pairs]", 400, -200, 200, 400, -200, 200);
    TH2I* h2zx = new TH2I(
        "zx", "z-x clusters;x [#mum];z [#mum];clusters [eh-pairs]", 4 * pitch, -2 * pitch, 2 * pitch, depth, 0, depth);

    TH1I hdtime("dtime", "drift time;drift time [ns];clusters", 100, 0, 10 + 20 * j);
    TH1I hdiff("diff", "diffusion width;diffusion width [#mum];clusters", 100, 0, 10);
    TH1I htf("tf", "Gaussian tail fraction;Gausian tail fraction;clusters", 100, 0, 1);
    TProfile tfvsx("tfvsx", "Gaussian tail fraction vs x;x [#mum];Gausian tail fraction", 200, -0.5 * pitch, 0.5 * pitch);

    TH1I hcleh("cleh", "cluster neh;log_{10}(cluster eh [pairs]);clusters", 80, 0, 4);
    TProfile wvse("wvse", "energy per eh pair;log_{10}(step E loss [eV]);<w> [eV/pair]", 80, 0, 4);
    TH1I hreh("reh", "eh/eV;eh/dE [pairs/eV];clusters", 160, 0, 0.8);

    TH1I heta0("eta0", "eta;eta;tracks", 201, -1.005, 1.005);
    TProfile eta0vsxm("eta0vsxm", "eta vs track;track x [#mum];<eta>", 200, -0.5 * pitch, 0.5 * pitch);
    TH1I hdx0("dx0", "dx0;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TProfile madx0vsq("madx0vsq", "MAD(#Deltax) vs charge;charge [ke];MAD(#Deltax) [#mum]", 100, 0, 2 * 0.1 * depth);
    TH1I hdx0q("dx0q", "dx0 Landau peak;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TProfile dx0qvsxm("dx0qvsxm", "#Deltax vs x;track x [#mum];<#Deltax> [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);
    TProfile madx0qvsxm(
        "madx0qvsxm", "MAD(#Deltax) vs x;track x [#mum];MAD(#Deltax) [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);
    TProfile madx0qvsxm0(
        "madx0qvsxm0", "MAD(#Deltax) vs x no delta;track x [#mum];MAD(#Deltax) [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);
    TProfile madx0qvsxm1(
        "madx0qvsxm1", "MAD(#Deltax) vs x delta;track x [#mum];MAD(#Deltax) [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);

    // threshold:
    TH1I hpxq1("pxq1",
               "thresholded pixel charge;pixel charge [ke];pixels",
               std::max(100, int(10 * 0.1 * depth / 1)),
               0,
               std::max(1, int(5 * 0.1 * depth / 1)));
    TH1I hq1("q1",
             "thresholded charge;charge [ke];tracks",
             std::max(100, int(50 * 0.1 * depth)),
             0,
             std::max(1, int(10 * 0.1 * depth)));
    TH1I hnpx1("npx1", "npx after threshold;npx;tracks", 4, 0.5, 4.5);
    TProfile npx1vsxm("npx1vsxm", "npx threshold vs track;track x [#mum];<npx>", 200, -0.5 * pitch, 0.5 * pitch);
    TH1I heta1("eta1", "eta threshold;eta;tracks", 201, -1.005, 1.005);
    TProfile eta1vsxm("eta1vsxm", "eta threshold vs track;track x [#mum];<eta>", 200, -0.5 * pitch, 0.5 * pitch);
    TProfile x1vsxm("x1vsxm", "xcog threshold vs track;track x [#mum];<cog> [#mum]", 200, -0.5 * pitch, 0.5 * pitch);
    TH1I hdx1("dx1", "dx threshold;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TProfile madx1vsq("madx1vsq", "MAD(#Deltax) vs charge;charge [ke];MAD(#Deltax) [#mum]", 100, 0, 2 * 0.1 * depth);
    TH1I hdx1q("dx1q", "dx threshold Landau peak;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TProfile dx1qvsxm("dx1qvsxm", "#Deltax vs x;track x [#mum];<#Deltax> [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);
    TProfile madx1qvsxm(
        "madx1qvsxm", "MAD(#Deltax) vs x;track x [#mum];MAD(#Deltax) [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);

    TH1I hda1("da1", "da threshold;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TH1I hda1q("da1q", "da threshold Landau peak;#Deltax [#mum];tracks", 501, -pitch * 1.001, pitch * 1.001);
    TProfile mada1qvsxm(
        "mada1qvsxm", "MAD(#Deltaa) vs x;track x [#mum];MAD(#Deltaa) [#mum]", pitch, -0.5 * pitch, 0.5 * pitch);

    std::ranlux24 rgen; // C++11 random number engine
    if(seed != 0) {
        std::cout << "SEEDING with " << seed << std::endl;
        rgen.seed(seed);
    } else {
        rgen.seed(time(NULL)); // seconds since 1.1.1970
    }

    DepositionBichsel module;
    module.setRandomEngine(&rgen);
    module.init();

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    std::uniform_real_distribution<double> unirnd(0, 1);
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // EVENT LOOP:

    for(unsigned iev = 0; iev < nev; ++iev) {

        std::cout << iev << std::endl;

        // put track on std::stack:
        double xm = pitch * (unirnd(rgen) - 0.5); // [mu] -p/2..p/2 at track mid
        ROOT::Math::XYZVector pos((xm - 0.5 * width) * 1e-4, 0, 0);
        ROOT::Math::XYZVector dir(sin(turn), 0, cos(turn));
        particle initial(Ekin0, pos, dir, default_particle_type); // beam particle is first "delta"
        // E : Ekin0; // [MeV]
        // x : entry point is left;
        // y :  [cm]
        // z :  pixel from 0 to depth [cm]

        unsigned ndelta = 0; // number of deltas generated

        auto clusters = module.stepping(initial, iev, depth, ndelta);

        // 4 pixels along x:

        double q1[4];
        for(int ir = 0; ir < 4; ++ir) {
            q1[ir] = 0;
        }

        for(unsigned i = 0; i < clusters.size(); ++i) {

            double xx = clusters[i].position.X() * 1e4; // [mu]
            double yy = clusters[i].position.Y() * 1e4; // [mu]
            double zz = clusters[i].position.Z() * 1e4; // [mu]
            int neh = clusters[i].neh;

            if(iev < 11) {
                h1zev[iev]->Fill(zz, neh);
                h2zxev[iev]->Fill(xx, zz, neh);
            }
            h2xy->Fill(xx, yy, neh);
            h2zx->Fill(xx, zz, neh);

            double Eg = clusters[i].E;
            hcleh.Fill(log(neh) / log(10)); // per cluster
            hreh.Fill(neh / Eg);
            wvse.Fill(log(Eg) / log(10), Eg / neh); // dE per eh pair

            // 0 | 1 | 2 | 3, bins 0 and 3 are half-infinite
            // diffusion across x crack: |  |  |

            double xc = -pitch;    // left
            int m = 0;             // minus = left pixel
            int p = 1;             // plus = right pixel
            if(xx > 0.5 * pitch) { // nearer x-crack is right
                xc = pitch;
                m = 2;
                p = 3;
            } else if(xx > -0.5 * pitch) { // mid crack
                xc = 0;
                m = 1;
                p = 2;
            }

            double dtime = zz * 1e-4 / vd;                // [s] drift time along z (mean speed theorem)
            double diff = sqrt(2 * D * dtime) * 1e4;      // [mu] rms diffusion (projected or 3D?)
            double uu = -(xx - xc) / std::sqrt(2) / diff; // scaled diffusion distance for erfc
            double tf = 0.5 * erfc(uu);                   // upper Gaussian tail fraction

            hdtime.Fill(dtime * 1e9);
            hdiff.Fill(diff);
            htf.Fill(tf);
            tfvsx.Fill(xx, tf); // S-curve, x = 0 is a pixel boundary

            q1[p] += neh * tf;
            q1[m] += neh * (1 - tf);

        } // clusters

        double q0 = q1[0] + q1[1] + q1[2] + q1[3];
        double eta = (q1[2] - q1[1]) / (q1[2] + q1[1]); // central bins
        heta0.Fill(eta);
        eta0vsxm.Fill(xm, eta);

        double sumq0 = 0;
        double sumqx0 = 0;
        for(int ir = 0; ir < 4; ++ir) {
            sumq0 += q1[ir];
            sumqx0 += q1[ir] * (ir - 1.5); // -1.5, -0.5, 0.5, 1.5
        }
        double cog0 = sumqx0 / sumq0;
        double dx0 = cog0 * pitch - xm; // [mu]
        hdx0.Fill(dx0);

        madx0vsq.Fill(q0 * 1e-3, fabs(dx0)); // linear rise, too steep

        if(q0 < 95 * depth) { // keep 2/3 in 300 mu
            hdx0q.Fill(dx0);
            dx0qvsxm.Fill(xm, dx0);
            madx0qvsxm.Fill(xm, fabs(dx0)); // inverted U-shape
            if(ndelta)
                madx0qvsxm1.Fill(xm, fabs(dx0)); // inverted U-shape
            else
                madx0qvsxm0.Fill(xm, fabs(dx0)); // inverted U-shape
        }

        // after threshold:

        int npx = 0;
        double sumq1 = 0;
        double sumqx1 = 0;
        for(int ir = 0; ir < 4; ++ir) {
            if(q1[ir] > thr) {
                ++npx;
                hpxq1.Fill(q1[ir] * 1e-3); // [ke]
                sumq1 += q1[ir];
                sumqx1 += q1[ir] * (ir - 1.5); // -1.5, -0.5, 0.5, 1.5
            }
        }
        hq1.Fill(sumq1 * 1e-3); // [ke]
        hnpx1.Fill(npx);
        npx1vsxm.Fill(xm, npx);

        double eta1 = (q1[2] - q1[1]) / (q1[2] + q1[1]); // central bins
        heta1.Fill(eta1);
        eta1vsxm.Fill(xm, eta1);

        double cog1 = sumqx1 / sumq1;
        x1vsxm.Fill(xm, cog1 * pitch);
        double dx1 = cog1 * pitch - xm;
        hdx1.Fill(dx1);

        madx1vsq.Fill(sumq1 * 1e-3, fabs(dx1)); // linear rise, too steep

        if(sumq1 < 95 * depth) { // keep 2/3 in 300 mu
            hdx1q.Fill(dx1);
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
            dx1qvsxm.Fill(xm, dx1);
            madx1qvsxm.Fill(xm, fabs(dx1)); // inverted U-shape
        }

    } // events

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    histoFile->Write();
    histoFile->ls();
    histoFile->Close();
    std::cout << std::endl;

    std::cout << "done: events " << nev << std::endl;

    std::cout << "  particle type     " << default_particle_type << std::endl;
    std::cout << "  kinetic energy    " << Ekin0 << " MeV" << std::endl;
    std::cout << "  number of events  " << nev << std::endl;
    std::cout << "  pixel pitch       " << pitch << " um" << std::endl;
    std::cout << "  thickness         " << depth << " um" << std::endl;
    std::cout << "  incident angle    " << turn * 180 / M_PI << " deg" << std::endl;
    std::cout << "  track width       " << width << " um" << std::endl;
    std::cout << "  temperature       " << temperature << " K" << std::endl;
    std::cout << "  readout threshold " << thr << " e" << std::endl;
    std::cout << "  cross talk        " << cx * 100 << "%" << std::endl;

    std::cout << std::endl
              << ((j) ? "  holes" : "  electrons") << std::endl
              << "  mobility for " << Efield << " V/cm"
              << ": vm " << vm // cm/s = 100 um / ns
              << ", Ec " << Ec << ", mu0 " << mu0 << std::endl
              << "  beta " << beta << ", mu " << mu << ", v " << vd << " cm/s"
              << " = " << vd / 1e5 << " mu/ns" << std::endl
              << "  D " << D << ", rms " << sqrt(2 * D * 4e-9) * 1e4 << " mu" // for 4 ns drift
              << std::endl;

    std::cout << std::endl << "  " << histoFile->GetName() << std::endl;
    std::cout << std::endl;

} // main
