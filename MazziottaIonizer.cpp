#include "MazziottaIonizer.hpp"

using namespace ionizer;

MazziottaIonizer::MazziottaIonizer(std::ranlux24* random_engine) : random_engine_(random_engine) {
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // SHELL INITIALIZATION

    for(unsigned n = 1; n <= 4; ++n) {
        for(unsigned i = 1; i <= 9; ++i) {
            auger_prob_integral[n][i] = 0;
            auger_energy[n][i] = 0;
        }
    }

    // auger_prob_integral(KSH,J) = PROBABILITA" INTEGRALI DEI VARI PROCESSI DI
    // EMISSIONE AUGR DALLA SHELL KSH

    // auger_prob_integral(KSH,J) = 1 PER L"ULTIMO VALORE DI J

    // KSH = 4 --> SHELL K
    // KSH = 3 --> SHELL L1
    // KSH = 2 --> SHELL L23

    auger_prob_integral[4][1] = 0.1920;
    auger_prob_integral[4][2] = 0.3885 + auger_prob_integral[4][1];
    auger_prob_integral[4][3] = 0.2325 + auger_prob_integral[4][2];
    auger_prob_integral[4][4] = 0.0720 + auger_prob_integral[4][3];
    auger_prob_integral[4][5] = 0.0030 + auger_prob_integral[4][4];
    auger_prob_integral[4][6] = 0.1000 + auger_prob_integral[4][5];
    auger_prob_integral[4][7] = 0.0040 + auger_prob_integral[4][6];
    auger_prob_integral[4][8] = 0.0070 + auger_prob_integral[4][7];
    auger_prob_integral[4][9] = 0.0010 + auger_prob_integral[4][8];
    auger_prob_integral[3][1] = 0.0250;
    auger_prob_integral[3][2] = 0.9750 + auger_prob_integral[3][1];
    auger_prob_integral[2][1] = 0.9990;
    auger_prob_integral[2][2] = 0.0010 + auger_prob_integral[2][1];

    // auger_energy[KSH, J) = ENERGIA IN eV DELL"ELETTRONE AUGR
    // EMESSO DALLA SHELL KSH NEL PROCESSO J

    auger_energy[4][1] = 1541.6;
    auger_energy[4][2] = 1591.1;
    auger_energy[4][3] = 1640.6;
    auger_energy[4][4] = 1690.3;
    auger_energy[4][5] = 1690.3;
    auger_energy[4][6] = 1739.8;
    auger_energy[4][7] = 1739.8;
    auger_energy[4][8] = 1839.0;
    auger_energy[4][9] = 1839.0;

    auger_energy[3][1] = 148.7;
    auger_energy[3][2] = 49.5;

    auger_energy[2][1] = 99.2;
    auger_energy[2][2] = 0.0;
}

double MazziottaIonizer::get_uniform_prn() {
    if(random_engine_ == nullptr) {
        exit(1);
    }

    return uniform_dist(*random_engine_);
}

std::stack<double> MazziottaIonizer::getIonization(double energy_gamma) {

    // INPUT:
    // EG = VIRTUAL GAMMA ENERGY [eV]

    // OUTPUT:
    // veh ENERGIES OF PRIMARY e/h
    std::stack<double> veh;

    // energy_gamma = energy_gamma - energy_gap; // double counting?

    // EV = binding ENERGY OF THE TOP OF THE VALENCE BAND
    const double energy_valence = energy_shell[1]; // 12.0 eV

    int is = -1;
    if(energy_gamma <= energy_shell[1]) {
        is = 0;
    } else if(energy_gamma <= EPP[3]) {
        is = 1;
    } else {
        double PV[5];
        if(energy_gamma > EPP[13]) {
            PV[1] = PM[13];
            PV[2] = PL23[13];
            PV[3] = PL1[13];
            PV[4] = PK[13];
        } else {
            unsigned iep = 3;
            for(; iep < 13; ++iep) {
                if(energy_gamma > EPP[iep] && energy_gamma <= EPP[iep + 1]) {
                    break;
                }
            }

            // interpolate:
            PV[1] = PM[iep] + (PM[iep + 1] - PM[iep]) / (EPP[iep + 1] - EPP[iep]) * (energy_gamma - EPP[iep]);
            PV[2] = PL23[iep] + (PL23[iep + 1] - PL23[iep]) / (EPP[iep + 1] - EPP[iep]) * (energy_gamma - EPP[iep]);
            PV[3] = PL1[iep] + (PL1[iep + 1] - PL1[iep]) / (EPP[iep + 1] - EPP[iep]) * (energy_gamma - EPP[iep]);
            PV[4] = PK[iep] + (PK[iep + 1] - PK[iep]) / (EPP[iep + 1] - EPP[iep]) * (energy_gamma - EPP[iep]);
        }

        double PPV = PV[1] + PV[2] + PV[3] + PV[4];

        PV[2] = PV[1] + PV[2];
        PV[3] = PV[2] + PV[3];
        PV[4] = PV[3] + PV[4];

        double rs = get_uniform_prn();
        unsigned iv = 1;
        for(; iv <= 4; ++iv) {
            PV[iv] = PV[iv] / PPV; // normalize
            if(PV[iv] > rs) {
                break;
            }
        }
        // Restrict to 4 maximum
        is = std::min(iv, 4u);
    }

    // cout << "  shells for " << energy_gamma << " eV, energy_valence " << energy_valence << ", is " << is << std::endl;

    // PROCESSES:

    // PHOTOABSORPTION IN VALENCE BAND
    if(is <= 1) {
        if(energy_gamma < 0.1) {
            return veh;
        }

        double rv = get_uniform_prn();
        if(energy_gamma < energy_valence) {
            veh.push(rv * energy_gamma);
            veh.push((1 - rv) * energy_gamma);
        } else {
            veh.push(rv * energy_valence);
            veh.push(energy_gamma - rv * energy_valence);
        }
        return veh;
    }

    // PHOTOABSORPTION IN AN INNER SHELL
    double Ephe = energy_gamma - energy_shell[is];
    if(Ephe <= 0) {
        std::cout << "shells: photoelectron with negative energy " << energy_gamma << ", shell " << is << " at "
                  << energy_shell[is] << " eV" << std::endl;
        return veh;
    }

    // PRIMARY PHOTOELECTRON:
    veh.push(Ephe);

    // AUGER ELECTRONS:
    double raug = get_uniform_prn();

    int ks = 1;
    if(is <= 1) {
        ks = 1;
    } else if(is <= 3) {
        if(raug > auger_prob_integral[is][1]) {
            ks = 2;
        }
    } else {
        if(raug >= auger_prob_integral[is][1]) {
            for(int js = 2; js <= nvac[is]; ++js) {
                if(raug >= auger_prob_integral[is][js - 1] && raug < auger_prob_integral[is][js]) {
                    ks = js;
                }
            }
        }
    }

    if(is == 2) {
        // L23-SHELL VACANCIES
        if(ks == 1) {
            // TRANSITION L23 M M
            transition(energy_valence, auger_energy[2][1], veh);
        }
    } else if(is == 3) {
        // L1-SHELL VACANCIES
        if(ks == 2) {
            // TRANSITION L1 L23 M
            double energy = energy_valence * get_uniform_prn();
            veh.push(energy);
            veh.push(auger_energy[is][ks] - energy);

            if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                // TRANSITION L23 M M
                transition(energy_valence, auger_energy[2][1], veh);
            }
        } else {
            // TRANSITION L1 M M
            transition(energy_valence, auger_energy[3][1], veh);
        }
    } else if(is == 4) {
        // K-SHELL VACANCIES
        if(ks >= 8) {
            // TRANSITION K M M
            transition(energy_valence, auger_energy[is][ks], veh);
        } else if(ks == 6 || ks == 7) {
            // TRANSITION K L23 M
            double energy = energy_valence * get_uniform_prn();
            veh.push(energy);
            veh.push(auger_energy[is][ks] - energy); // adjust for energy conservation

            if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                // TRANSITION L23 M M
                transition(energy_valence, auger_energy[2][1], veh);
            }
        } else if(ks == 4 || ks == 5) {
            // TRANSITION K L1 M
            double energy = energy_valence * get_uniform_prn();
            veh.push(energy);
            veh.push(auger_energy[is][ks] - energy); // adjust for energy conservation

            if(get_uniform_prn() <= auger_prob_integral[3][1]) {
                // TRANSITION L1 M M
                transition(energy_valence, auger_energy[3][1], veh);
            } else {
                // TRANSITION L1 L23 M
                double energy = energy_valence * get_uniform_prn();
                veh.push(energy);
                veh.push(auger_energy[3][2] - energy);

                if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                    // TRANSITION L23 M M
                    transition(energy_valence, auger_energy[2][1], veh);
                }
            }
        } else if(ks == 3) {
            // TRANSITION K L23 L23
            veh.push(auger_energy[is][ks]); // default

            if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                // TRANSITION L23 M M
                transition(energy_valence, auger_energy[2][1], veh);
            }

            if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                // TRANSITION L23 M M
                transition(energy_valence, auger_energy[2][1], veh);
            }
        } else if(ks == 2) {
            // TRANSITION K L1 L23
            veh.push(auger_energy[is][ks]); // default

            // L23-SHELL VACANCIES
            if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                // TRANSITION L23 M M
                transition(energy_valence, auger_energy[2][1], veh);
            }

            // L1-SHELL VACANCIES
            if(get_uniform_prn() > auger_prob_integral[3][1]) {
                // TRANSITION L1 L23 M
                double energy = energy_valence * get_uniform_prn();
                veh.push(energy);
                veh.push(auger_energy[3][2] - energy);

                if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                    // TRANSITION L23 M M
                    transition(energy_valence, auger_energy[2][1], veh);
                }
            } else {
                // TRANSITION L1 M M
                transition(energy_valence, auger_energy[3][1], veh);
            }
        } else if(ks == 1) {
            // TRANSITION K L1 L1
            veh.push(auger_energy[is][ks]); // default

            // L1-SHELL VACANCIES
            if(get_uniform_prn() > auger_prob_integral[3][1]) {
                // TRANSITION L1 L23 M
                double energy = energy_valence * get_uniform_prn();
                veh.push(energy);
                veh.push(auger_energy[3][2] - energy);

                // L23-SHELL VACANCIES
                if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                    // TRANSITION L23 M M
                    transition(energy_valence, auger_energy[2][1], veh);
                }
            } else {
                // TRANSITION L1 M M
                transition(energy_valence, auger_energy[3][1], veh);
            }

            // L1-SHELL VACANCIES
            if(get_uniform_prn() > auger_prob_integral[3][1]) {
                // TRANSITION L1 L23 M
                double energy = energy_valence * get_uniform_prn();
                veh.push(energy);
                veh.push(auger_energy[3][2] - energy);

                // L23-SHELL
                if(get_uniform_prn() <= auger_prob_integral[2][1]) {
                    // TRANSITION L23 M M
                    transition(energy_valence, auger_energy[2][1], veh);
                }
            } else {
                // TRANSITION L1 M M
                transition(energy_valence, auger_energy[3][1], veh);
            }
        } // ks
    }     // is 4

    return veh;
} // SHELLS

void MazziottaIonizer::transition(double energy_valence, double energy_auger, std::stack<double>& veh) {

    auto gentri = [&]() { // -1..1
        double x = 0, r2 = 0, trian = 0;
        do {
            x = -1 + 2 * get_uniform_prn();
            r2 = get_uniform_prn();
            trian = x < 0 ? x + 1 : -x + 1;
        } while(trian > r2); // rejection method

        return x;
    };

    // AUGER ELECTRON
    double rEv = (1 + gentri()) * energy_valence; // 0..2*Ev
    veh.push(energy_auger - rEv);

    // ASSIGN ENERGIES TO THE HOLES
    double energy_hole1 = 0, energy_hole2 = 0;
    do {
        double rv = get_uniform_prn();
        energy_hole1 = rv * rEv;
        energy_hole2 = (1 - rv) * rEv;
    } while(energy_hole1 > energy_valence || energy_hole2 > energy_valence); // holes stay below valence band edge (12 eV)

    veh.push(energy_hole1);
    veh.push(energy_hole2);
}
