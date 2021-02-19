#include "DepositionBichsel.hpp"
#include "MazziottaIonizer.hpp"

#include <fstream>
#include <sstream>

using namespace allpix;

void DepositionBichsel::init() {

    // Booking histograms:

    elvse = new TProfile("elvse", "elastic mfp;log_{10}(E_{kin}[MeV]);elastic mfp [#mum]", 140, -3, 4);
    invse = new TProfile("invse", "inelastic mfp;log_{10}(E_{kin}[MeV]);inelastic mfp [#mum]", 140, -3, 4);

    hstep5 = new TH1I("step5", "step length;step length [#mum];steps", 500, 0, 5);
    hstep0 = new TH1I("step0", "step length;step length [#mum];steps", 500, 0, 0.05);
    hzz = new TH1I("zz", "z;depth z [#mum];steps", DEPTH, 0, DEPTH);

    hde0 = new TH1I("de0", "step E loss;step E loss [eV];steps", 200, 0, 200);
    hde1 = new TH1I("de1", "step E loss;step E loss [eV];steps", 100, 0, 5000);
    hde2 = new TH1I("de2", "step E loss;step E loss [keV];steps", 200, 0, 20);
    hdel = new TH1I("del", "log step E loss;log_{10}(step E loss [eV]);steps", 140, 0, 7);
    htet = new TH1I("tet", "delta emission angle;delta emission angle [deg];inelasic steps", 180, 0, 90);
    hnprim = new TH1I("nprim", "primary eh;primary e-h;scatters", 21, -0.5, 20.5);
    hlogE = new TH1I("logE", "log Eeh;log_{10}(E_{eh}) [eV]);eh", 140, 0, 7);
    hlogn = new TH1I("logn", "log neh;log_{10}(n_{eh});clusters", 80, 0, 4);
    hscat = new TH1I("scat", "elastic scattering angle;scattering angle [deg];elastic steps", 180, 0, 180);
    hncl = new TH1I("ncl", "clusters;e-h clusters;tracks", 4 * DEPTH * 5, 0, 4 * DEPTH * 5);

    double lastbin = EKIN < 1.1 ? 1.05 * EKIN * 1e3 : 5 * 0.35 * DEPTH; // 350 eV/micron
    htde = new TH1I("tde", "sum E loss;sum E loss [keV];tracks / keV", std::max(100, int(lastbin)), 0, int(lastbin));
    htde0 = new TH1I(
        "tde0", "sum E loss, no delta;sum E loss [keV];tracks, no delta", std::max(100, int(lastbin)), 0, int(lastbin));
    htde1 = new TH1I(
        "tde1", "sum E loss, with delta;sum E loss [keV];tracks, with delta", std::max(100, int(lastbin)), 0, int(lastbin));

    hteh = new TH1I("teh",
                    "total e-h;total charge [ke];tracks",
                    std::max(100, int(50 * 0.1 * DEPTH)),
                    0,
                    std::max(1, int(10 * 0.1 * DEPTH)));
    hq0 = new TH1I("q0",
                   "normal charge;normal charge [ke];tracks",
                   std::max(100, int(50 * 0.1 * DEPTH)),
                   0,
                   std::max(1, int(10 * 0.1 * DEPTH)));
    hrms = new TH1I("rms", "RMS e-h;charge RMS [e];tracks", 100, 0, 50 * DEPTH);

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // INITIALIZE ENERGY BINS

    double u = log(2.0) / N2;
    double um = exp(u);
    int ken = log(1839.0 / 1.5) / u;         // intger
    double Emin = 1839.0 / pow(2, ken / N2); // integer division intended

    // EMIN is chosen to give an E-value exactly at the K-shell edge, 1839 eV
    E[1] = Emin;
    for(size_t j = 2; j < E.size(); ++j) {
        E[j] = E[j - 1] * um; // must agree with heps.tab
        dE[j - 1] = E[j] - E[j - 1];
    }

    std::cout << std::endl << "  n2 " << N2 << ", Emin " << Emin << ", um " << um << ", E[nume] " << E.back() << std::endl;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // READ DIELECTRIC CONSTANTS
    read_hepstab();

    //= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // READ INTEGRAL OVER MOMENTUM TRANSFER OF THE GENERALIZED OSCILLATOR STRENGTH
    read_macomtab();

    read_emerctab();

    // EGAP = GAP ENERGY IN eV
    // EMIN = THRESHOLD ENERGY (ALIG ET AL., PRB22 (1980), 5565)
    energy_gap = 1.17 - 4.73e-4 * temperature * temperature / (636 + temperature);
    energy_threshold = 1.5 * energy_gap; // energy conservation
    temperature = TEMPERATURE;
    default_particle_type = PARTICLE_TYPE;
    explicit_delta_energy_cut_keV = DELTA_ENERGY_CUT;
    fast = FAST;
}

std::vector<cluster> DepositionBichsel::stepping(const particle& init, unsigned iev, double depth, unsigned& ndelta) {

    MazziottaIonizer ionizer(random_engine_);
    std::uniform_real_distribution<double> unirnd(0, 1);

    std::vector<cluster> clusters;
    std::stack<particle> deltas;
    deltas.push(init);
    // Statistics:
    unsigned nsteps = 0; // number of steps for full event
    unsigned nscat = 0;  // elastic scattering
    unsigned nloss = 0;  // ionization
    double total_energy_loss = 0.0;
    unsigned nehpairs = 0;
    unsigned sumeh2 = 0;

    while(!deltas.empty()) {
        double Ekprev = 9e9; // update flag for next delta

        auto t = deltas.top();
        deltas.pop();

        unsigned nlast = E.size() - 1;

        double xm0 = 1;
        double xlel = 1;
        double gn = 1;
        table totsig;

        std::cout << "  delta " << t.E() * 1e3 << " keV"
                  << ", cost " << t.direction().Z() << ", u " << t.direction().X() << ", v " << t.direction().Y() << ", z "
                  << t.position().Z() * 1e4;

        while(1) { // steps

            if(t.E() < 0.9 * Ekprev) { // update

                // Emax = maximum energy loss, see Uehling, also Sternheimer & Peierls Eq.(53)
                double Emax = t.mass() * (t.gamma() * t.gamma() - 1) /
                              (0.5 * t.mass() / electron_mass + 0.5 * electron_mass / t.mass() + t.gamma());

                // maximum energy loss for incident electrons
                if(t.type() == ParticleType::ELECTRON) {
                    Emax = 0.5 * t.E();
                }
                Emax = 1e6 * Emax; // eV

                // Define parameters and calculate Inokuti"s sums,
                // S ect 3.3 in Rev Mod Phys 43, 297 (1971)

                double zi = 1.0;
                double dec = zi * zi * atnu * fac / t.betasquared();
                double EkeV = t.E() * 1e6; // [eV]

                // Generate collision spectrum sigma(E) from df/dE, epsilon and AE.
                // sig(*,j) actually is E**2 * sigma(E)

                std::array<double, 6> tsig;
                tsig.fill(0);
                table H;

                // Statistics:
                double stpw = 0;

                std::array<table, 6> sig;
                for(unsigned j = 1; j < E.size(); ++j) {

                    if(E[j] > Emax) {
                        break;
                    }

                    // Eq. (3.1) in RMP and red notebook CCS-33, 39 & 47
                    double Q1 = rydberg_constant;
                    if(E[j] < 11.9) {
                        Q1 = pow(xkmn[j], 2) * rydberg_constant;
                    } else if(E[j] < 100.0) {
                        Q1 = pow(0.025, 2) * rydberg_constant;
                    }

                    double qmin = E[j] * E[j] / (2 * electron_mass * 1e6 * t.betasquared()); // twombb = 2 m beta**2 [eV]
                    if(E[j] < 11.9 && Q1 < qmin) {
                        sig[1][j] = 0;
                    } else {
                        sig[1][j] = E[j] * dfdE[j] * log(Q1 / qmin);
                    }
                    // longitudinal excitation, Eq. (46) in Fano; Eq. (2.9) in RMP

                    double epbe = std::max(1 - t.betasquared() * dielectric_const_real[j], 1e-20); // Fano Eq. (47)
                    double sgg =
                        E[j] * dfdE[j] * (-0.5) * log(epbe * epbe + pow(t.betasquared() * dielectric_const_imag[j], 2));

                    double thet = atan(dielectric_const_imag[j] * t.betasquared() / epbe);
                    if(thet < 0) {
                        thet = thet + M_PI; // plausible-otherwise I"d have a jump
                    }
                    // Fano says [p 21]: "arctan approaches pi for betasq*eps1 > 1 "

                    double sgh = 0.0092456 * E[j] * E[j] * thet *
                                 (t.betasquared() - dielectric_const_real[j] / (pow(dielectric_const_real[j], 2) +
                                                                                pow(dielectric_const_imag[j], 2)));

                    sig[2][j] = sgg;
                    sig[3][j] = sgh; // small, negative

                    // uef from  Eqs. 9 & 2 in Uehling, Ann Rev Nucl Sci 4, 315 (1954)
                    double uef = 1 - E[j] * t.betasquared() / Emax;
                    if(t.type() == ParticleType::ELECTRON) {
                        uef = 1 + pow(E[j] / (EkeV - E[j]), 2) + pow((t.gamma() - 1) / t.gamma() * E[j] / EkeV, 2) -
                              (2 * t.gamma() - 1) * E[j] / (t.gamma() * t.gamma() * (EkeV - E[j]));
                    }
                    // there is a factor of 2 because the integral was over d(lnK) rather than d(lnQ)
                    sig[4][j] = 2 * oscillator_strength_ae[j] * uef;

                    sig[5][j] = 0;
                    for(unsigned i = 1; i <= 4; ++i) {
                        sig[5][j] += sig[i][j]; // sum

                        // divide by E**2 to get the differential collision cross section sigma
                        // Tsig = integrated total collision cross section
                        tsig[i] = tsig[i] + sig[i][j] * dE[j] / (E[j] * E[j]);
                    }                                             // i
                    tsig[5] += sig[5][j] * dE[j] / (E[j] * E[j]); // running sum

                    double HE2 = sig[5][j] * dec;
                    H[j] = HE2 / (E[j] * E[j]);
                    stpw += H[j] * E[j] * dE[j]; // dE/dx
                    nlast = j;
                }
                xm0 = tsig[5] * dec; // 1/path

                // Statistics:
                double sst = H[1] * dE[1]; // total cross section (integral)

                totsig[1] = H[1] * dE[1]; // running integral
                for(unsigned j = 2; j <= nlast; ++j) {
                    totsig[j] = totsig[j - 1] + H[j] * dE[j];
                    sst += H[j] * dE[j];
                }

                // NORMALIZE running integral:
                for(unsigned j = 1; j <= nlast; ++j) {
                    totsig[j] /= totsig[nlast]; // norm
                }

                // elastic:
                if(t.type() == ParticleType::ELECTRON) {
                    // gn = 2*2.61 * pow( atomic_number, 2.0/3.0 ) / EkeV; // Mazziotta
                    gn = 2 * 2.61 * pow(atomic_number, 2.0 / 3.0) / (t.momentum() * t.momentum()) * 1e-6; // Moliere
                    double E2 = 14.4e-14;                                                                 // [MeV*cm]
                    double FF = 0.5 * M_PI * E2 * E2 * atomic_number * atomic_number / (t.E() * t.E());
                    double S0EL = 2 * FF / (gn * (2 + gn));
                    // elastic total cross section  [cm2/atom]
                    xlel = atnu * S0EL; // ATNU = N_A * density / A = atoms/cm3
                } else {
                    double getot = t.E() + t.mass();
                    xlel = std::min(2232.0 * radiation_length * pow(t.momentum() * t.momentum() / (getot * zi), 2),
                                    10.0 * radiation_length);
                    // units ?
                }

                elvse->Fill(log(t.E()) / log(10), 1e4 / xlel);
                invse->Fill(log(t.E()) / log(10), 1e4 / xm0);

                Ekprev = t.E();

                if(ldb)
                    std::cout << "  ev " << iev << " type " << t.type() << ", Ekin " << t.E() * 1e3 << " keV"
                              << ", beta " << sqrt(t.betasquared()) << ", gam " << t.gamma() << std::endl
                              << "  Emax " << Emax << ", nlast " << nlast << ", Elast " << E[nlast] << ", norm "
                              << totsig[nlast] << std::endl
                              << "  inelastic " << 1e4 / xm0 << "  " << 1e4 / sst << ", elastic " << 1e4 / xlel << " um"
                              << ", mean dE " << stpw * depth * 1e-4 * 1e-3 << " keV" << std::endl
                              << std::flush;
            } // update

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // step:

            double tlam = 1 / (xm0 + xlel);                           // [cm] TOTAL MEAN FREE PATH (MFP)
            double step = -log(1 - unirnd(getRandomEngine())) * tlam; // exponential step length
            double pos_z = t.position().Z() + step * t.direction().Z();

            if(ldb && t.E() < 1) {
                std::cout << "step " << step * 1e4 << ", z " << pos_z * 1e4 << std::endl;
            }
            hstep5->Fill(step * 1e4);
            hstep0->Fill(step * 1e4);
            hzz->Fill(pos_z * 1e4);

            // Outside the sensor
            if(pos_z < 0 || pos_z > depth * 1e-4) {
                break; // exit back or front
            }

            // Update position after step
            t.setPosition(t.position() + step * t.direction());

            if(fabs(t.position().Y()) > 0.0200) {
                break; // save time
            }

            ++nsteps;

            // INELASTIC (ionization) PROCESS
            if(unirnd(getRandomEngine()) > tlam * xlel) {
                ++nloss;

                // GENERATE VIRTUAL GAMMA:
                double yr = unirnd(getRandomEngine()); // inversion method
                unsigned je = 2;
                for(; je <= nlast; ++je) {
                    if(yr < totsig[je]) {
                        break;
                    }
                }

                double energy_gamma = E[je - 1] + (E[je] - E[je - 1]) * unirnd(getRandomEngine()); // [eV]

                hde0->Fill(energy_gamma); // M and L shells
                hde1->Fill(energy_gamma); // K shell
                hde2->Fill(energy_gamma * 1e-3);
                hdel->Fill(log(energy_gamma) / log(10));

                double residual_kin_energy = t.E() - energy_gamma * 1E-6; // [ MeV]

                // cut off for further movement: [MeV]
                if(residual_kin_energy < explicit_delta_energy_cut_keV * 1e-3) {
                    energy_gamma = t.E() * 1E6;                 // [eV]
                    residual_kin_energy = t.E() - energy_gamma; // zero
                    // std::cout << "LAST ENERGY LOSS" << energy_gamma << residual_kin_energy
                }

                total_energy_loss += energy_gamma; // [eV]

                // emission angle from delta:

                // PRIMARY SCATTERING ANGLE:
                // SINT = SQRT(ER*1.E-6/EK) ! Mazziotta = deflection angle
                // COST = SQRT(1.-SINT*SINT)
                // STORE INFORMATION ABOUT DELTA-RAY:
                // SINT = COST ! flip
                // COST = SQRT(1.-SINT**2) ! sqrt( 1 - ER*1e-6 / t.E() ) ! wrong

                // double cost = sqrt( energy_gamma / (2*electron_mass_ev + energy_gamma) ); // M. Swartz
                double cost =
                    sqrt(energy_gamma / (2 * electron_mass * 1e6 + energy_gamma) * (t.E() + 2 * electron_mass) / t.E());
                // Penelope, Geant4
                double sint = 0;
                if(cost * cost <= 1) {
                    sint = sqrt(1 - cost * cost); // mostly 90 deg
                }
                double phi = 2 * M_PI * unirnd(getRandomEngine());

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

                std::vector<double> din{sint * cos(phi), sint * sin(phi), cost};
                htet->Fill(180 / M_PI * asin(sint)); // peak at 90, tail to 45, elastic forward

                // transform into detector system:
                double cz = t.direction().Z(); // delta direction
                double sz = sqrt(1 - cz * cz);
                double phif = atan2(t.direction().Y(), t.direction().X());
                ROOT::Math::XYZVector delta_direction(cz * cos(phif) * din[0] - sin(phif) * din[1] + sz * cos(phif) * din[2],
                                                      cz * sin(phif) * din[0] + cos(phif) * din[1] + sz * sin(phif) * din[2],
                                                      -sz * din[0] + cz * din[2]);

                // GENERATE PRIMARY e-h:
                std::stack<double> veh;
                if(energy_gamma > energy_threshold) {
                    veh = ionizer.getIonization(energy_gamma);
                }

                hnprim->Fill(veh.size());

                double sumEeh{0};
                unsigned neh{0};

                // PROCESS e and h
                while(!veh.empty()) {
                    double Eeh = veh.top();
                    veh.pop();

                    hlogE->Fill(Eeh > 1 ? log(Eeh) / log(10) : 0);

                    if(Eeh > explicit_delta_energy_cut_keV * 1e3) {
                        // Put new delta on std::stack:
                        deltas.emplace(Eeh * 1E-6, t.position(), delta_direction, ParticleType::ELECTRON);

                        ++ndelta;
                        total_energy_loss -= Eeh; // [eV], avoid double counting

                        continue; // next ieh
                    }             // new delta

                    sumEeh += Eeh;

                    // slow down low energy e and h: 95% of CPU time
                    while(!fast && Eeh > energy_threshold) {
                        // for e and h
                        double p_ionization =
                            1 / (1 + aaa * 105 / 2 / M_PI * sqrt(Eeh - eom0) / pow(Eeh - energy_threshold, 3.5));

                        if(unirnd(getRandomEngine()) < p_ionization) { // ionization
                            ++neh;
                            double E1 = gena1() * (Eeh - energy_threshold);
                            double E2 = gena2() * (Eeh - energy_threshold - E1);
                            // cout << "      ion " << Eeh << " => " << E1 << " + " << E2 << std::endl;

                            if(E1 > energy_threshold) {
                                veh.push(E1);
                            }
                            if(E2 > energy_threshold) {
                                veh.push(E2);
                            }

                            Eeh = Eeh - E1 - E2 - energy_threshold;
                        } else {
                            Eeh -= eom0; // phonon emission
                        }
                    } // slow: while Eeh
                }     // while veh

                if(fast) {
                    std::poisson_distribution<int> poisson(sumEeh / 3.645);
                    neh = poisson(getRandomEngine());
                }

                nehpairs += neh;
                sumeh2 += neh * neh;

                // cout << "  dE " << energy_gamma << " eV, neh " << neh << std::endl;

                // Store charge cluster:
                if(neh > 0) {
                    clusters.emplace_back(neh, t.position(), energy_gamma);

                    hlogn->Fill(log(neh) / log(10));
                }

                // Update particle energy
                t.setE(t.E() - energy_gamma * 1E-6); // [MeV]

                if(ldb && t.E() < 1) {
                    std::cout << "    Ek " << t.E() * 1e3 << " keV, z " << t.position().Z() * 1e4 << ", neh " << neh
                              << ", steps " << nsteps << ", ion " << nloss << ", elas " << nscat << ", cl "
                              << clusters.size() << std::endl;
                }

                if(t.E() < 1E-6 || residual_kin_energy < 1E-6) {
                    // std::cout << "  absorbed" << std::endl;
                    break;
                }

                // For electrons, update elastic cross section at new energy
                if(t.type() == ParticleType::ELECTRON) {
                    // gn = 2*2.61 * pow( atomic_number, 2.0/3.0 ) / (t.E()*1E6); // Mazziotta
                    double pmom = sqrt(t.E() * (t.E() + 2 * t.mass()));                   // [MeV/c] 2nd binomial
                    gn = 2 * 2.61 * pow(atomic_number, 2.0 / 3.0) / (pmom * pmom) * 1e-6; // Moliere
                    double E2 = 14.4e-14;                                                 // [MeV*cm]
                    double FF = 0.5 * M_PI * E2 * E2 * atomic_number * atomic_number / (t.E() * t.E());
                    double S0EL = 2 * FF / (gn * (2 + gn));
                    // elastic total cross section  [cm2/atom]
                    xlel = atnu * S0EL; // ATNU = N_A * density / A = atoms/cm3
                }

            } else { // ELASTIC SCATTERING: Chaoui 2006
                ++nscat;

                double r = unirnd(getRandomEngine());
                double cost = 1 - 2 * gn * r / (2 + gn - 2 * r);
                double sint = sqrt(1 - cost * cost);
                double phi = 2 * M_PI * unirnd(getRandomEngine());
                std::vector<double> din{sint * cos(phi), sint * sin(phi), cost};

                hscat->Fill(180 / M_PI * asin(sint)); // forward peak, tail to 90

                // Change direction of particle:
                double cz = t.direction().Z(); // delta direction
                double sz = sqrt(1.0 - cz * cz);
                double phif = atan2(t.direction().Y(), t.direction().X());
                t.setDirection(ROOT::Math::XYZVector(cz * cos(phif) * din[0] - sin(phif) * din[1] + sz * cos(phif) * din[2],
                                                     cz * sin(phif) * din[0] + cos(phif) * din[1] + sz * sin(phif) * din[2],
                                                     -sz * din[0] + cz * din[2]));
            } // elastic
        }     // while steps

        std::cout << std::endl;
    } // while deltas

    std::cout << "  steps " << nsteps << ", ion " << nloss << ", elas " << nscat << ", dE " << total_energy_loss * 1e-3
              << " keV"
              << ", eh " << nehpairs << ", cl " << clusters.size() << std::endl;

    hncl->Fill(clusters.size());
    htde->Fill(total_energy_loss * 1e-3); // [keV] energy conservation - binding energy
    if(ndelta) {
        htde1->Fill(total_energy_loss * 1e-3); // [keV]
    } else {
        htde0->Fill(total_energy_loss * 1e-3); // [keV]
    }
    hteh->Fill(nehpairs * 1e-3); // [ke]
    hq0->Fill(nehpairs * 1e-3);  // [ke]
    hrms->Fill(sqrt(sumeh2));

    return clusters;
}

void DepositionBichsel::read_hepstab() {
    std::ifstream heps("HEPS.TAB");
    if(heps.bad() || !heps.is_open()) {
        std::cout << "Error opening HEPS.TAB" << std::endl;
        exit(1);
    }

    std::string line;
    getline(heps, line);
    std::istringstream tokenizer(line);

    int n2t;
    unsigned numt;
    tokenizer >> n2t >> numt;

    std::cout << " HEPS.TAB: n2t " << n2t << ", numt " << numt << std::endl;
    if(N2 != n2t) {
        std::cout << " CAUTION: n2 & n2t differ" << std::endl;
    }
    if(E.size() - 1 != numt) {
        std::cout << " CAUTION: nume & numt differ" << std::endl;
    }
    if(numt > E.size() - 1) {
        numt = E.size() - 1;
    }

    unsigned jt = 1;
    while(!heps.eof() && jt < numt) {
        getline(heps, line);
        std::istringstream tokenizer(line);

        double etbl, rimt;
        tokenizer >> jt >> etbl >> dielectric_const_real[jt] >> dielectric_const_imag[jt] >> rimt;

        // The dipole oscillator strength df/dE is calculated, essentially Eq. (2.20)
        dfdE[jt] = rimt * 0.0092456 * E[jt];
    }

    std::cout << "read " << jt << " data lines from HEPS.TAB" << std::endl;
    // MAZZIOTTA: 0.0 at 864
    // EP( 2, 864 ) = 0.5 * ( EP(2, 863) + EP(2, 865) )
    // RIM(864) = 0.5 * ( RIM(863) + RIM(865) )
    // DFDE(864) = RIM(864) * 0.0092456 * E(864)
    // DP: fixed in HEPS.TAB
}

void DepositionBichsel::read_macomtab() {
    std::ifstream macom("MACOM.TAB");
    if(macom.bad() || !macom.is_open()) {
        std::cout << "Error opening MACOM.TAB" << std::endl;
        exit(1);
    }

    std::string line;
    getline(macom, line);
    std::istringstream tokenizer(line);

    int n2t;
    unsigned numt;
    tokenizer >> n2t >> numt;

    unsigned int nume = E.size() - 1;
    std::cout << " MACOM.TAB: n2t " << n2t << ", numt " << numt << std::endl;
    if(N2 != n2t)
        std::cout << " CAUTION: n2 & n2t differ" << std::endl;
    if(nume != numt)
        std::cout << " CAUTION: nume & numt differ" << std::endl;
    if(numt > nume)
        numt = nume;

    unsigned jt = 1;
    while(!macom.eof() && jt < numt) {
        getline(macom, line);
        std::istringstream tokenizer(line);

        double etbl;
        tokenizer >> jt >> etbl >> oscillator_strength_ae[jt];
    }
    std::cout << "read " << jt << " data lines from MACOM.TAB" << std::endl;
}

void DepositionBichsel::read_emerctab() {
    std::ifstream emerc("EMERC.TAB");
    if(emerc.bad() || !emerc.is_open()) {
        std::cout << "Error opening EMERC.TAB" << std::endl;
        exit(1);
    }

    std::string line;
    getline(emerc, line); // header lines
    getline(emerc, line);
    getline(emerc, line);
    getline(emerc, line);

    std::istringstream tokenizer(line);

    unsigned jt = 1;
    while(!emerc.eof() && jt < 200) {
        getline(emerc, line);
        std::istringstream tokenizer(line);

        double etbl;
        tokenizer >> jt >> etbl >> oscillator_strength_ae[jt] >> xkmn[jt];
    }
    std::cout << "  read " << jt << " data lines from EMERC.TAB" << std::endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double DepositionBichsel::gena1() {
    std::uniform_real_distribution<double> uniform_dist(0, 1);

    double r1 = 0, r2 = 0, alph1 = 0;
    do {
        r1 = uniform_dist(getRandomEngine());
        r2 = uniform_dist(getRandomEngine());
        alph1 = 105. / 16. * (1. - r1) * (1 - r1) * sqrt(r1); // integral = 1, max = 1.8782971
    } while(alph1 > 1.8783 * r2);                             // rejection method

    return r1;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double DepositionBichsel::gena2() {
    std::uniform_real_distribution<double> uniform_dist(0, 1);

    double r1 = 0, r2 = 0, alph2 = 0;
    do {
        r1 = uniform_dist(getRandomEngine());
        r2 = uniform_dist(getRandomEngine());
        alph2 = 8 / M_PI * sqrt(r1 * (1 - r1));
    } while(alph2 > 1.27324 * r2); // rejection method

    return r1;
}
