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

    /**
     * @brief Particle
     */
    class particle {
    public:
        /**
         * Constructor for new particle
         * @param energy        Kinetic energy of the particle
         * @param pos           Position of generation
         * @param dir           Direction of motion
         * @param particle_type Type of particle
         */
        particle(double energy, ROOT::Math::XYZVector pos, ROOT::Math::XYZVector dir, ParticleType particle_type)
            : E(energy), position(std::move(pos)), direction(std::move(dir)), type(particle_type){};

        /**
         * Default constructor
         */
        particle() = default;
        double E; // [MeV]
        ROOT::Math::XYZVector position;
        ROOT::Math::XYZVector direction;
        ParticleType type; // particle type

        /**
         * Helper to obtain particle rest mass in units of MeV
         * @return Particle rest mass in MeV
         */
        double mass() { return mass_.at(static_cast<std::underlying_type<ParticleType>::type>(type)); };

    private:
        std::vector<double> mass_{
            0,
            938.2723,   // proton
            139.578,    // pion
            493.67,     // K
            0.51099906, // e
            105.65932   // mu
        };
    };

    /**
     * @brief Deposited clusters of electron-hole pairs generated via ionization
     */
    class cluster {
    public:
        /**
         * @ brief default constructor
         */
        cluster() = default;

        /**
         * Constructor for e/h pair cluster
         * @param eh_pairs Number of electron-hole pairs
         * @param pos      Position of the cluster in local coordinates
         * @param energy   Energy of the generating particle
         */
        cluster(int eh_pairs, ROOT::Math::XYZVector pos, double energy) : neh(eh_pairs), position(pos), E(energy){};
        int neh;
        ROOT::Math::XYZVector position;
        double E; // [eV] generating particle
    };

    /**
     * Reading HEPS.TAB data file
     *
     * HEPS.TAB is the table of the dielectric constant for solid Si, epsilon = ep(1,j) + i*ep(2,j), as a function of energy
     * loss E(j), section II.E in RMP, and rim is Im(-1/epsilon), Eq. (2.17), p.668. Print statements are included to check
     * that the file is read correctly.
     *
     * @param n2   [description]
     * @param E    Energy levels
     * @param ep_1 Array of real parts of complex dielectric constant, split by energy levels
     * @param ep_2 Array of imaginary parts of complex dielectric constant, split by energy levels
     * @param dfdE [description]
     */
    void read_hepstab(int n2, ionizer::table E, ionizer::table& ep_1, ionizer::table& ep_2, ionizer::table& dfdE);

    /**
     * Reading MACOM.TAB data file
     *
     * MACOM.TAB is the table of the integrals over momentum transfer K of the generalized oscillator strength, summed for
     * all shells, i.e. the A(E) of Eq. (2.11), p. 667 of RMP
     *
     * @param n2   [description]
     * @param nume Number of energy levels
     * @param sig  Array of integrals of generalized oscillator strength, split by energy level
     */
    void read_macomtab(int n2, unsigned int nume, ionizer::table& sig);

    /**
     * Reading EMERC.TAB data file
     *
     * EMERC.TAB is the table of the integral over K of generalized oscillator strength for E < 11.9 eV with Im(-1/epsilon)
     * from equations in the Appendix of Emerson et al., Phys Rev B7, 1798 (1973) (also see CCS-63)
     *
     * @param sig  Array of integrals of generalized oscillator strength, split by energy level
     * @param xkmn [description]
     */
    void read_emerctab(ionizer::table& sig, ionizer::table& xkmn);

} // namespace ionizer
