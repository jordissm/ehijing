#include "event_output.hpp"
#include "spacetime_history.hpp"
#include "HepMC3/ReaderAscii.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

Pythia8::Particle particle(Pythia8::ParticleDataEntry& data, int sign = 1,
                          int status = 83) {
    Pythia8::Particle p(sign * data.id(), status, 0, 0, 0, 0, 0, 0,
                       0., 0., 0., data.m0(), data.m0());
    p.setPDEPtr(&data);
    return p;
}
} // namespace

int main() {
    try {
        using Pythia8::ParticleDataEntry;
        ParticleDataEntry proton(2212, "p+", "pbar", 2, 3, 0, 0.938);
        ParticleDataEntry neutron(2112, "n", "nbar", 2, 0, 0, 0.939);
        ParticleDataEntry delta_plus(2214, "Delta+", "anti-Delta+", 4, 3, 0, 1.232);
        ParticleDataEntry delta_plus_plus(2224, "Delta++", "anti-Delta++", 4, 6, 0, 1.232);
        ParticleDataEntry delta_minus(1114, "Delta-", "anti-Delta-", 4, -3, 0, 1.232);
        ParticleDataEntry pion(211, "pi+", "pi-", 1, 3, 0, 0.138);
        ParticleDataEntry muon(13, "mu-", "mu+", 2, -3, 0, 0.106);
        ParticleDataEntry quark(2, "u", "ubar", 2, 2, 1, 0.33);

        auto balance = calculate_spectator_balance(6, 12, StruckNucleon::Proton,
                                                   {particle(proton)});
        require(balance.admissible && balance.protons == 5 && balance.neutrons == 6,
                "Primary proton bookkeeping changed");
        require(balance.baryon_number_conserved && balance.electric_charge_conserved,
                "Primary proton conservation failed");
        balance = calculate_spectator_balance(6, 12, StruckNucleon::Neutron,
                                             {particle(neutron)});
        require(balance.admissible && balance.protons == 6 && balance.neutrons == 5,
                "Primary neutron bookkeeping changed");

        // A baryon-antibaryon pair contributes no net B or Q; a final charged
        // meson does contribute Q. Intermediate particles must not be counted.
        balance = calculate_spectator_balance(6, 12, StruckNucleon::Neutron,
            {particle(neutron), particle(delta_plus), particle(delta_plus, -1),
             particle(pion), particle(proton, 1, -83)});
        require(balance.final_baryon_number_times3 == 3 && balance.final_charge_times3 == 3,
                "Antibaryon, resonance, meson or intermediate counting failed");
        require(balance.admissible && balance.protons == 5 && balance.neutrons == 6,
                "Charged meson was omitted from spectator balance");

        balance = calculate_spectator_balance(1, 1, StruckNucleon::Proton,
                                             {particle(neutron), particle(muon, -1)});
        require(balance.admissible && balance.protons == 0 && balance.neutrons == 0,
                "Charge of final decay lepton was omitted");
        balance = calculate_spectator_balance(1, 2, StruckNucleon::Proton,
                                             {particle(proton), particle(neutron)});
        require(balance.admissible && balance.protons == 0 && balance.neutrons == 0,
                "Fully depleted target needs zero spectators");
        balance = calculate_spectator_balance(1, 4, StruckNucleon::Proton,
                                             {particle(delta_plus_plus)});
        require(!balance.admissible && balance.required_protons_times3 == -3
                && balance.required_neutrons_times3 == 12,
                "Excess charge must be non-admissible even when baryon count fits");
        balance = calculate_spectator_balance(1, 1, StruckNucleon::Proton,
                                             {particle(neutron)});
        require(!balance.admissible && balance.required_protons_times3 == 3
                && balance.required_neutrons_times3 == -3,
                "Negative neutron requirement was missed");

        // A fractional required count must not be rounded into apparent conservation.
        balance = calculate_spectator_balance(6, 12, StruckNucleon::Proton,
                                             {particle(quark)});
        require(!balance.admissible && balance.protons == 5 && balance.neutrons == 6,
                "Fractional balance was accepted or legacy fallback changed");

        DISKinematics kin{};
        std::vector<Pythia8::Particle> products = {
            particle(proton), particle(proton), particle(neutron)};
        balance = calculate_spectator_balance(6, 12, StruckNucleon::Proton, products);
        require(balance.admissible && balance.protons == 4 && balance.neutrons == 5,
                "Medium baryons were not removed from spectators");
        std::ostringstream oscar, metadata;
        SpacetimeHistory history(42, 1);
        history.capture_particles(products, "hadronization_and_decays", "hadronizer",
                                  "absolute_mm_saved_vertex");
        write_event_output(42, 6, 12, StruckNucleon::Proton, kin, products,
                           oscar, metadata, &history);
        // Parse the actual emitted rows: count, charges, baryons and IDs must
        // agree with the event header and the independently specified target.
        std::istringstream lines(oscar.str());
        std::string line;
        int rows = 0, baryons = 0, charge = 0, spectator_protons = 0, spectator_neutrons = 0;
        while (std::getline(lines, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            double t, x, y, z, m, e, px, py, pz, production, formation, factor;
            int pdg, id, q;
            require(bool(row >> t >> x >> y >> z >> m >> e >> px >> py >> pz
                             >> pdg >> id >> q >> production >> formation >> factor),
                    "Invalid OSCAR row");
            require(id == rows, "Output IDs are not contiguous");
            charge += q;
            baryons += pdg == 2212 || pdg == 2112;
            if (factor == 1.) {
                spectator_protons += pdg == 2212;
                spectator_neutrons += pdg == 2112;
            }
            ++rows;
        }
        require(rows == 12 && baryons == 12 && charge == 6, "Emitted event violates B or Q");
        require(spectator_protons == 4 && spectator_neutrons == 5, "Wrong emitted spectators");
        require(oscar.str().find("# event 42 out 12\n") == 0, "OSCAR header omits spectators");
        require(metadata.str().find("\"admissible\": true") != std::string::npos,
                "Metadata does not report admissibility");
        std::ostringstream ascii;
        history.write(ascii);
        std::istringstream input(ascii.str());
        HepMC3::ReaderAscii reader(input);
        HepMC3::GenEvent snapshot;
        reader.read_event(snapshot);
        reader.read_event(snapshot);
        require(!reader.failed(), "Updated history is not HepMC3 ASCII");
        require(snapshot.attribute_as_string("target_A") == "12", "History omits target A");
        require(snapshot.attribute_as_string("nuclear_bookkeeping").find(
                    "\"electric_charge_conserved\": true") != std::string::npos,
                "History omits conservation status");

        // Explicit spectators retain their identities/geometry and never use
        // count repair to remove or invent nucleons. Repeated hits count once.
        EHIJING::NuclearTarget target(2,1,.5);
        target.sample([](){return .25;});
        require(target.primary()==0,"Fixture primary must be proton");
        target.nucleons()[1].position={1.,2.,3.};
        products={particle(proton)};
        balance=calculate_spectator_balance(1,2,StruckNucleon::Proton,products);
        require(target_balance_matches(target,balance),"Untouched identity rejected");
        kin.pProton=Pythia8::Vec4(0.,0.,0.,.938);
        std::ostringstream target_oscar,target_meta;
        SpacetimeHistory target_history(44,1);
        write_event_output(44,1,2,StruckNucleon::Proton,kin,products,target_oscar,target_meta,&target_history,&target);
        require(target_oscar.str().find("0.19733 0.39465 0.59198")!=std::string::npos,
                "Original spectator position discarded");
        require(target_meta.str().find("untouched_target_identities")!=std::string::npos,
                "Explicit spectator policy missing");
        target.nucleons()[1].hits=5;
        require(!target_balance_matches(target,balance),"Extra baryon silently repaired with spectator deletion");
        products={particle(proton),particle(neutron)};
        balance=calculate_spectator_balance(1,2,StruckNucleon::Proton,products);
        require(target_balance_matches(target,balance),"Repeated hits counted as extra baryons");
        // Match the user's deuteron example: B_H=5, Q_H=2 cannot be repaired by spectators.
        products = {particle(proton), particle(proton), particle(delta_plus),
                    particle(delta_minus), particle(delta_minus), particle(pion)};
        balance = calculate_spectator_balance(1, 2, StruckNucleon::Proton, products);
        require(!balance.admissible && balance.required_protons_times3 == -3
                && balance.required_neutrons_times3 == -6, "Impossible deuteron balance missed");
        require(balance.protons == 0 && balance.neutrons == 1,
                "Non-admissible event handling was changed");
        require(!balance.baryon_number_conserved && !balance.electric_charge_conserved,
                "Impossible balance falsely reported as conserved");
        std::ostringstream invalid_oscar, invalid_meta;
        write_event_output(43, 1, 2, StruckNucleon::Proton, kin, products,
                           invalid_oscar, invalid_meta);
        require(invalid_meta.str().find("\"admissible\": false") != std::string::npos
                && invalid_meta.str().find("\"spectator_policy\": \"legacy_unmodified\"")
                    != std::string::npos, "Impossible event fallback not disclosed");
        require(invalid_oscar.str().find("# event 43 out 7\n") == 0,
                "Non-admissible event was rejected or its spectators changed");
        std::cout << "Spectator baryon and charge balance checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
