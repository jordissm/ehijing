#include "event_output.hpp"
#include "spacetime_history.hpp"
#include "ehijing_constants.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>

using namespace Pythia8;

namespace {

namespace constants = ehijing::constants;

double sample_point_in_sphere(double R, std::mt19937& gen) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double u = uniform(gen);
    return R * std::cbrt(u);
}

int count_final_particles(const std::vector<Particle>& particles) {
    int count = 0;
    for (const auto& p : particles) {
        if (p.isFinal()) {
            ++count;
        }
    }
    return count;
}

std::string bookkeeping_json(const SpectatorBalance& balance, bool explicit_target = false) {
    std::ostringstream out;
    out << std::setprecision(16)
        << "{\"admissible\": " << (balance.admissible ? "true" : "false")
        << ", \"final_baryon_number\": " << balance.final_baryon_number_times3 / 3.0
        << ", \"final_charge\": " << balance.final_charge_times3 / 3.0
        << ", \"required_spectator_protons\": " << balance.required_protons_times3 / 3.0
        << ", \"required_spectator_neutrons\": " << balance.required_neutrons_times3 / 3.0
        << ", \"spectator_protons\": " << balance.protons
        << ", \"spectator_neutrons\": " << balance.neutrons
        << ", \"output_baryon_number\": "
        << balance.final_baryon_number_times3 / 3.0 + balance.protons + balance.neutrons
        << ", \"output_charge\": " << balance.final_charge_times3 / 3.0 + balance.protons
        << ", \"baryon_number_conserved\": " << (balance.baryon_number_conserved ? "true" : "false")
        << ", \"electric_charge_conserved\": " << (balance.electric_charge_conserved ? "true" : "false")
        << ", \"spectator_policy\": \"" << (explicit_target ? "untouched_target_identities" : balance.admissible ? "conservation" : "legacy_unmodified")
        << "\"}";
    return out.str();
}

void write_metadata_json(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    const DISKinematics& kin,
    const std::string& bookkeeping,
    const std::string& nuclear_target,
    std::ostream& out) {

    out << std::setprecision(16);
    out << "{\n";
    out << "  \"event\": " << event_number << ",\n";
    out << "  \"Z\": " << atomic_number << ", \"A\": " << mass_number << ",\n";
    out << "  \"xB\": " << kin.bjorken_x << ",\n";
    out << "  \"Q2\": " << kin.Q2 << ",\n";
    out << "  \"y\": "  << kin.y  << ",\n";
    out << "  \"nu\": " << kin.nu << ",\n";
    out << "  \"P4\": [" << kin.pProton.px() << ", "
                        << kin.pProton.py() << ", "
                        << kin.pProton.pz() << ", "
                        << kin.pProton.e()  << "],\n";
    out << "  \"q4\": [" << kin.pGamma.px() << ", "
                        << kin.pGamma.py() << ", "
                        << kin.pGamma.pz() << ", "
                        << kin.pGamma.e()  << "],\n";
    out << "  \"nuclear_bookkeeping\": " << bookkeeping;
    if (!nuclear_target.empty()) out << ",\n  \"nuclear_target\": " << nuclear_target;
    out << '\n';
    out << "}\n";
}

void write_final_hadrons(
    const std::vector<Particle>& particles,
    std::ostream& out,
    int32_t& particle_index,
    SpacetimeHistory* history) {

    for (const auto& p : particles) {
        if (!p.isFinal()) {
            continue;
        }

        // Convert the complete formation vertex from PYTHIA's mm to output fm
        // before freestreaming to t = 0.
        Vec4 local_pos = p.vProd() * constants::units::mm_to_fm;
        // Vec4 p_com = p.p();
        // Vec4 phadron = p.p();

        // phadron.bstback(p_com);
        // local_pos.bstback(p_com);

        const double form_time = local_pos.e();
        // Production uses the stored references, including origin-only collapses.
        const double prod_time = p.hasProductionTime()
            ? p.tProduction() * constants::units::mm_to_fm
            : 0.0;
        const int prod_time_status = p.hasProductionTime() ? 1 : 0;
        if (p.hasProductionTime()
            && (!std::isfinite(p.vProductionFirst().e())
                || !std::isfinite(p.vProductionSecond().e())
                || !std::isfinite(prod_time))) {
            throw std::invalid_argument("Non-finite production time for PDG "
                                        + std::to_string(p.id()));
        }
        double x_hadron = local_pos.px();
        double y_hadron = local_pos.py();
        double z_hadron = local_pos.pz();

        const double mass = p.m();
        const double e = p.e();
        const int pdgid = p.id();
        const int pid = particle_index;
        const int charge = p.charge();

        // Freestream from the saved formation point to t = 0.
        // For form_time > 0 this moves it backward in time; for form_time < 0
        // the same expression moves it forward in time.
        if (std::isfinite(form_time) && std::isfinite(e) && e > 0.0) {
            const double dt_to_zero = -form_time;
            x_hadron += p.px() / e * dt_to_zero;
            y_hadron += p.py() / e * dt_to_zero;
            z_hadron += p.pz() / e * dt_to_zero;
        }

        out << std::fixed << std::setprecision(2) << 0.00 << " "
            << std::fixed << std::setprecision(5)
            << x_hadron << " "
            << y_hadron << " "
            << z_hadron << " "
            << mass << " "
            << e << " "
            << p.px() << " "
            << p.py() << " "
            << p.pz() << " "
            << pdgid << " "
            << pid << " "
            << charge << " "
            << prod_time << " "
            << form_time << " "
            << 0.0 << " "
            << prod_time_status << '\n';

        if (history) history->output_particle(p, particle_index,
            Vec4(x_hadron, y_hadron, z_hadron, 0.0) * constants::units::fm_to_mm);
        ++particle_index;
    }
}

void write_spectator_nucleons(
    int mass_number,
    int64_t total_protons,
    int64_t total_neutrons,
    std::ostream& out,
    int32_t& particle_index,
    SpacetimeHistory* history) {

    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const double RA = constants::nuclear::radius_coefficient_fm *
                      std::pow(static_cast<double>(mass_number), 1.0 / 3.0);

    for (int64_t i = 0; i < total_protons; ++i) {
        const double rr = sample_point_in_sphere(RA, gen);
        const double phi = constants::math::two_pi * uniform(gen);
        const double costheta = 1.0 - 2.0 * uniform(gen);
        const double sintheta = std::sqrt(1.0 - costheta * costheta);

        const double rx = rr * sintheta * std::cos(phi);
        const double ry = rr * sintheta * std::sin(phi);
        const double rz = rr * costheta;

        const double t = 0.0;
        const double px = 0.0;
        const double py = 0.0;
        const double pz = 0.0;
        const double p0 = std::sqrt(
            constants::mass::proton_gev *
            constants::mass::proton_gev + px * px + py * py + pz * pz);

        out << std::fixed << std::setprecision(2) << t << " "
            << std::fixed << std::setprecision(5)
            << rx << " "
            << ry << " "
            << rz << " "
            << constants::mass::proton_gev << " "
            << p0 << " "
            << px << " "
            << py << " "
            << pz << " "
            << constants::pdg::proton_id << " "
            << particle_index << " "
            << 1 << " "
            << 0.0 << " "
            << 0.0 << " "
            << 1.0 << " 2\n";

        if (history) history->output_particle(
            Particle(constants::pdg::proton_id, 1, 0, 0, 0, 0, 0, 0,
                     px, py, pz, p0, constants::mass::proton_gev), particle_index,
            Vec4(rx, ry, rz, t) * constants::units::fm_to_mm, true);
        ++particle_index;
    }

    for (int64_t i = 0; i < total_neutrons; ++i) {
        const double rr = sample_point_in_sphere(RA, gen);
        const double phi = constants::math::two_pi * uniform(gen);
        const double costheta = 1.0 - 2.0 * uniform(gen);
        const double sintheta = std::sqrt(1.0 - costheta * costheta);

        const double rx = rr * sintheta * std::cos(phi);
        const double ry = rr * sintheta * std::sin(phi);
        const double rz = rr * costheta;

        const double t = 0.0;
        const double px = 0.0;
        const double py = 0.0;
        const double pz = 0.0;
        const double p0 = std::sqrt(
            constants::mass::neutron_gev *
            constants::mass::neutron_gev + px * px + py * py + pz * pz);

        out << std::fixed << std::setprecision(2) << t << " "
            << std::fixed << std::setprecision(5)
            << rx << " "
            << ry << " "
            << rz << " "
            << constants::mass::neutron_gev << " "
            << p0 << " "
            << px << " "
            << py << " "
            << pz << " "
            << constants::pdg::neutron_id << " "
            << particle_index << " "
            << 0 << " "
            << 0.0 << " "
            << 0.0 << " "
            << 1.0 << " 2\n";

        if (history) history->output_particle(
            Particle(constants::pdg::neutron_id, 1, 0, 0, 0, 0, 0, 0,
                     px, py, pz, p0, constants::mass::neutron_gev), particle_index,
            Vec4(rx, ry, rz, t) * constants::units::fm_to_mm, true);
        ++particle_index;
    }
}

void write_target_spectators(const EHIJING::NuclearTarget& target,
    const DISKinematics& kin, std::ostream& out, int32_t& index, SpacetimeHistory* history) {
    for (size_t id = 0; id < target.nucleons().size(); ++id) {
        const auto& n = target.nucleons()[id]; if (n.hits != 0) continue;
        const double mass = n.proton ? constants::mass::proton_gev : constants::mass::neutron_gev;
        Vec4 p(0.,0.,0.,mass); p.bst(kin.pProton);
        Vec4 v(n.position[0],n.position[1],n.position[2],0.);
        v *= constants::units::hbarc_gev_fm; v.bst(kin.pProton);
        // Common t=0 output hypersurface in the generator/output frame.
        const double x = v.px()-p.px()/p.e()*v.e();
        const double y = v.py()-p.py()/p.e()*v.e();
        const double z = v.pz()-p.pz()/p.e()*v.e();
        const int pid = n.proton ? 2212 : 2112;
        out << std::fixed << std::setprecision(5) << 0. << ' ' << x << ' ' << y << ' ' << z
            << ' ' << mass << ' ' << p.e() << ' ' << p.px() << ' ' << p.py() << ' ' << p.pz()
            << ' ' << pid << ' ' << index << ' ' << (n.proton ? 1 : 0) << " 0 0 1 2\n";
        if (history) {
            Particle spectator(pid,1,0,0,0,0,0,0,p,mass);
            spectator.nuclearTargetId(id);
            history->output_particle(spectator,index,Vec4(x,y,z,0.)*constants::units::fm_to_mm,true);
        }
        ++index;
    }
}
} // namespace

bool target_balance_matches(const EHIJING::NuclearTarget& target, const SpectatorBalance& balance) {
    int64_t protons = 0, neutrons = 0;
    for (const auto& n : target.nucleons()) if (n.hits == 0) { if (n.proton) ++protons; else ++neutrons; }
    return balance.admissible && balance.protons == protons && balance.neutrons == neutrons;
}

SpectatorBalance calculate_spectator_balance(
    int atomic_number, int mass_number, StruckNucleon struck_nucleon,
    const std::vector<Particle>& particles) {
    SpectatorBalance balance;
    for (const auto& p : particles) {
        if (!p.isFinal()) continue;
        // isHadron/isParton also guard the ParticleDataEntry pointer. PYTHIA's
        // baryonNumberType counts signed net quarks, including antibaryons.
        if (p.isHadron() || p.isParton())
            balance.final_baryon_number_times3 += p.particleDataEntry().baryonNumberType(p.id());
        balance.final_charge_times3 += p.chargeType();
    }
    balance.required_protons_times3 = 3LL * atomic_number - balance.final_charge_times3;
    balance.required_neutrons_times3 = 3LL * mass_number - balance.final_baryon_number_times3
                                      - balance.required_protons_times3;
    balance.admissible = balance.required_protons_times3 >= 0
        && balance.required_neutrons_times3 >= 0
        && balance.required_protons_times3 % 3 == 0
        && balance.required_neutrons_times3 % 3 == 0;
    if (balance.admissible) {
        balance.protons = balance.required_protons_times3 / 3;
        balance.neutrons = balance.required_neutrons_times3 / 3;
    } else {
        // No rejection, resampling or clipping policy is introduced here.
        balance.protons = atomic_number - (struck_nucleon == StruckNucleon::Proton ? 1 : 0);
        balance.neutrons = mass_number - atomic_number
                          - (struck_nucleon == StruckNucleon::Neutron ? 1 : 0);
    }
    balance.baryon_number_conserved = balance.final_baryon_number_times3
        + 3 * (balance.protons + balance.neutrons) == 3LL * mass_number;
    balance.electric_charge_conserved = balance.final_charge_times3
        + 3 * balance.protons == 3LL * atomic_number;
    return balance;
}

void write_event_headers(std::ostream& out) {
    out << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge prod_time form_time xsecfac prod_time_status\n";
    out << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none fm fm none none\n";
}

void write_event_output(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    StruckNucleon struck_nucleon,
    const DISKinematics& kin,
    const std::vector<Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out,
    SpacetimeHistory* history,
    const EHIJING::NuclearTarget* target) {

    const auto balance = calculate_spectator_balance(atomic_number, mass_number,
                                                     struck_nucleon, particles);
    if (target && !target_balance_matches(*target,balance))
        throw std::logic_error("Explicit target B/Q mismatch must be rejected before output");
    const auto bookkeeping = bookkeeping_json(balance,target != nullptr);
    std::ostringstream ledger;
    if (target) {
        ledger << std::setprecision(17) << "{\"primary_id\":" << target->primary()
            << ",\"width_fm\":" << target->width()*constants::units::hbarc_gev_fm
            << ",\"radius_fm\":" << target->radius()*constants::units::hbarc_gev_fm
            << ",\"surface_fm\":" << target->surface()*constants::units::hbarc_gev_fm
            << ",\"frame\":\"nuclear_rest\",\"nucleons\":[";
        for (size_t id=0;id<target->nucleons().size();++id) {
            const auto& n=target->nucleons()[id];
            if(id) ledger << ',';
            ledger << "{\"id\":" << id << ",\"pdg\":" << (n.proton?2212:2112)
                << ",\"hits\":" << n.hits << ",\"position_GeV_inv\":["
                << n.position[0] << ',' << n.position[1] << ',' << n.position[2] << "]}";
        }
        ledger << "]}";
    }
    write_metadata_json(event_number, atomic_number, mass_number, kin, bookkeeping,ledger.str(), meta_out);

    // The OSCAR header counts every emitted final product and spectator.
    const int64_t count = count_final_particles(particles) + balance.protons + balance.neutrons;
    // Write the event header (OSCAR-style)
    event_out << "# event " << event_number << " out " << count << "\n";

    int32_t particle_index = 0;
    write_final_hadrons(particles, event_out, particle_index, history);
    if (target) write_target_spectators(*target,kin,event_out,particle_index,history);
    else write_spectator_nucleons(mass_number, balance.protons, balance.neutrons,
                            event_out, particle_index, history);
    if (history) {
        history->finish_output();
        history->output_attribute("nuclear_bookkeeping", bookkeeping);
        if (target) history->output_attribute("nuclear_target",ledger.str());
        history->output_attribute("target_A", std::to_string(mass_number));
        history->output_attribute("target_Z", std::to_string(atomic_number));
    }

    // Write event footer
    event_out << "# event " << event_number << " end 0\n";
}
