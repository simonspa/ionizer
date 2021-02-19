#include <random>
#include <stack>

#include <Math/Vector3D.h>

namespace ionizer {

    #define HEPS_ENTRIES 1251
    #define N2 64
    using table = std::array<double, HEPS_ENTRIES>;

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
}
