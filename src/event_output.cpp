#include "event_output.hpp"
#include "ehijing_constants.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <random>
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

void write_metadata_json(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    const DISKinematics& kin,
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
                        << kin.pGamma.e()  << "]\n";
    out << "}\n";
}

void write_final_hadrons(
    const std::vector<Particle>& particles,
    std::ostream& out,
    int32_t& particle_index) {

    for (const auto& p : particles) {
        if (!p.isFinal()) {
            continue;
        }

        Vec4 local_pos{p.xProd(), p.yProd(), p.zProd(), p.tProd()};
        // Vec4 p_com = p.p();
        // Vec4 phadron = p.p();

        // phadron.bstback(p_com);
        // local_pos.bstback(p_com);

        const double t_hadron = local_pos.e();
        double x_hadron = local_pos.px();
        double y_hadron = local_pos.py();
        double z_hadron = local_pos.pz();

        const double mass = p.m();
        const double e = p.e();
        const int pdgid = p.id();
        const int pid = particle_index;
        const int charge = p.charge();

        // Freestream each hadron along its lab-frame velocity to t = 0.
        // For t_hadron > 0 this moves it backward in time; for t_hadron < 0
        // the same expression moves it forward in time.
        if (std::isfinite(t_hadron) && std::isfinite(e) && e > 0.0) {
            const double dt_to_zero = -t_hadron;
            x_hadron += p.px() / e * dt_to_zero;
            y_hadron += p.py() / e * dt_to_zero;
            z_hadron += p.pz() / e * dt_to_zero;
        }

        // Smear transverse output positions to avoid exact spatial overlaps.
        static thread_local std::mt19937 smear_gen(std::random_device{}());
        std::normal_distribution<double> transverse_smear(0.0, 0.5);
        x_hadron += transverse_smear(smear_gen);
        y_hadron += transverse_smear(smear_gen);

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
            << t_hadron << " "
            << 0.0 << '\n';

        ++particle_index;
    }
}

void write_spectator_nucleons(
    int atomic_number,
    int mass_number,
    StruckNucleon struck_nucleon,
    std::ostream& out,
    int32_t& particle_index) {

    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const double RA = constants::nuclear::radius_coefficient_fm *
                      std::pow(static_cast<double>(mass_number), 1.0 / 3.0);

    const int total_protons =
        atomic_number - (struck_nucleon == StruckNucleon::Proton ? 1 : 0);
    for (int i = 0; i < total_protons; ++i) {
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
            << 1.0 << '\n';

        ++particle_index;
    }

    const int total_neutrons =
        mass_number - atomic_number -
        (struck_nucleon == StruckNucleon::Neutron ? 1 : 0);
    for (int i = 0; i < total_neutrons; ++i) {
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
            << 1.0 << '\n';

        ++particle_index;
    }
}

} // namespace

void write_event_headers(std::ostream& out) {
    out << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge begin_form_time xsecfac\n";
    out << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none fm none\n";
}

void write_event_output(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    StruckNucleon struck_nucleon,
    const DISKinematics& kin,
    const std::vector<Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out) {

    write_metadata_json(event_number, atomic_number, mass_number, kin, meta_out);

    // Count number of hadrons in the shower
    const int count = count_final_particles(particles);
    // Write the event header (OSCAR-style)
    event_out << "# event " << event_number << " out " << count << "\n";

    int32_t particle_index = 0;
    write_final_hadrons(particles, event_out, particle_index);
    write_spectator_nucleons(atomic_number, mass_number, struck_nucleon, event_out, particle_index);

    // Write event footer
    event_out << "# event " << event_number << " end 0\n";
}
