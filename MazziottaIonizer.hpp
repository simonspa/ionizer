#include <random>
#include <stack>

#include <Math/Vector3D.h>

namespace ionizer {

    class MazziottaIonizer {
    public:
        MazziottaIonizer(std::ranlux24* random_engine);
        std::stack <double> getIonization(double energy_gamma);

    private:
        void transition(double energy_valence, double energy_auger, std::stack <double> &veh);

        std::ranlux24* random_engine_{nullptr};
        std::uniform_real_distribution <double> uniform_dist{0, 1};
        double get_uniform_prn();

        // Shells
        // Possible transitions to this shell:
        const std::array<int, 5> nvac{{0, 0, 2, 2, 9}};
        // [1]: valence band upper edge (holes live below)
        // [2]: M, [3]: L, [4]: K shells
        const std::array<double, 5> energy_shell{{0, 12.0, 99.2, 148.7, 1839.0,}};

        double auger_prob_integral[5][10];
        double auger_energy[5][10];

        // EPP(I) = VALORI DI ENERGIA PER TABULARE LE PROBABILITA" DI FOTOASSORBIMENTO NELLE VARIE SHELL
        // PM, PL23, PL1, PK = PROBABILITA" DI ASSORBIMENTO DA PARTE DELLE SHELL M,L23,L1 E K

        // VALORI ESTRAPOLATI DA FRASER
        const std::array<double, 14> EPP{{0.0, 40.0, 50.0, 99.2, 99.2, 148.7, 148.7, 150.0, 300.0, 500.0, 1000.0, 1839.0, 1839.0, 2000.0}};
        const std::array<double, 14> PM{{0, 1.0, 1.0, 1.0, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.03, 0.05, 0.0, 0.0}};
        const std::array<double, 14> PL23{{0, 0.0, 0.0, 0.0, 0.97, 0.92, 0.88, 0.88, 0.83, 0.70, 0.55, 0.39, 0.0, 0.0}};
        const std::array<double, 14> PL1{{0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.15, 0.28, 0.42, 0.56, 0.08, 0.08}};
        const std::array<double, 14> PK{{0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.92, 0.92}};
    };
}
