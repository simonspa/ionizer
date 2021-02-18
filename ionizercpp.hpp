#include <random>
#include <stack>

namespace ionizer {

    class DepositionBichsel {
    public:
        DepositionBichsel();
        std::stack <double> shells(double energy_gamma);

    private:
        void transition(double energy_valence, double energy_auger, std::stack <double> &veh);

        double energy_shell[5];
        double auger_prob_integral[5][10];
        double auger_energy[5][10];
    };
}
