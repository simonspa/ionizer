#include <random>
#include <stack>

#include <Math/Vector3D.h>

namespace ionizer {

    /**
     * @brief Type of particles
     */
    enum class ParticleType {
        NONE = 0, ///< No particle
        PROTON,
        PION,
        KAON,
        ELECTRON,
        MUON,
        HELIUM,
        LITHIUM,
        CARBON,
        IRON,
    };

    inline std::ostream& operator<<(std::ostream& os, const ParticleType type) {
        os << static_cast<std::underlying_type<ParticleType>::type>(type);
        return os;
    }

    class particle {
    public:
        particle(double energy, ROOT::Math::XYZVector pos, ROOT::Math::XYZVector dir, ParticleType particle_type) : E(energy), position(std::move(pos)), direction(std::move(dir)), type(particle_type)
        {};
        particle() = default;
        double E; // [MeV]
        ROOT::Math::XYZVector position;
        ROOT::Math::XYZVector direction;
        ParticleType type; // particle type
        double mass() {
            return mass_.at(static_cast<std::underlying_type<ParticleType>::type>(type));
        };

    private:
        std::vector<double> mass_{0,
            938.2723, // proton
            139.578, // pion
            493.67, // K
            0.51099906, // e
            105.65932 // mu
        };
    };

    class cluster {
    public:
        cluster() = default;
        cluster(int eh_pairs, ROOT::Math::XYZVector pos, double energy) : neh(eh_pairs), position(pos), E(energy) {};
        int neh;
        ROOT::Math::XYZVector position;
        double E; // [eV] generating particle
    };

    class DepositionBichsel {
    public:
        DepositionBichsel();
        std::stack <double> shells(double energy_gamma);

    private:
        void transition(double energy_valence, double energy_auger, std::stack <double> &veh);

        int nvac[5];
        double energy_shell[5];
        double auger_prob_integral[5][10];
        double auger_energy[5][10];

        // EPP(I) = VALORI DI ENERGIA PER TABULARE LE PROBABILITA" DI FOTOASSORBIMENTO NELLE VARIE SHELL
        // PM, PL23, PL1, PK = PROBABILITA" DI ASSORBIMENTO DA PARTE DELLE SHELL M,L23,L1 E K

        // VALORI ESTRAPOLATI DA FRASER
        std::vector<double> EPP{0.0, 40.0, 50.0, 99.2, 99.2, 148.7, 148.7, 150.0, 300.0, 500.0, 1000.0, 1839.0, 1839.0, 2000.0};;
        std::vector<double> PM{0, 1.0, 1.0, 1.0, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.03, 0.05, 0.0, 0.0};
        std::vector<double> PL23{0, 0.0, 0.0, 0.0, 0.97, 0.92, 0.88, 0.88, 0.83, 0.70, 0.55, 0.39, 0.0, 0.0};
        std::vector<double> PL1{0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.15, 0.28, 0.42, 0.56, 0.08, 0.08};
        std::vector<double> PK{0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.92, 0.92};
    };
}
