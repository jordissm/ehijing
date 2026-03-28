#include "event_output.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <random>
#include <stdexcept>

using namespace Pythia8;

namespace {

constexpr double kFormationEpsilon = 1.0e-5;
constexpr double kProtonMass = 0.938272;
constexpr double kNeutronMass = 0.939565;

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
    int Z,
    int A,
    const DISKinematics& kin,
    std::ostream& out) {

    out << std::setprecision(16);
    out << "{\n";
    out << "  \"event\": " << event_number << ",\n";
    out << "  \"Z\": " << Z << ", \"A\": " << A << ",\n";
    out << "  \"xB\": " << kin.xBj << ",\n";
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
        Vec4 p_com = p.p();
        Vec4 phadron = p.p();

        // phadron.bstback(p_com);
        // local_pos.bstback(p_com);

        // Vec4 formation_4{
        //     kFormationEpsilon * phadron.px() / phadron.e() + local_pos.px(),
        //     kFormationEpsilon * phadron.py() / phadron.e() + local_pos.py(),
        //     kFormationEpsilon * phadron.pz() / phadron.e() + local_pos.pz(),
        //     kFormationEpsilon + local_pos.e()
        // };

        // formation_4.bst(p_com);

        // const double t_hadron = formation_4.e();
        // const double x_hadron = formation_4.px();
        // const double y_hadron = formation_4.py();
        // const double z_hadron = formation_4.pz();

        const double t_hadron = local_pos.e();
        const double x_hadron = local_pos.px();
        const double y_hadron = local_pos.py();
        const double z_hadron = local_pos.pz();

        const double mass = p.m();
        const double e = p.e();
        const int pdgid = p.id();
        const int pid = particle_index;
        const int charge = p.charge();

        out << std::fixed << std::setprecision(2) << t_hadron << " "
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
            << charge << '\n';

        ++particle_index;
    }
}

void write_spectator_nucleons(
    int Z,
    int A,
    std::ostream& out,
    int32_t& particle_index) {

    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_int_distribution<int> binary(0, 1);

    const double RA = 1.2 * std::pow(static_cast<double>(A), 1.0 / 3.0);
    const int random_binary = binary(gen);

    const int total_protons = Z - random_binary;
    for (int i = 0; i < total_protons; ++i) {
        const double rr = sample_point_in_sphere(RA, gen);
        const double phi = 2.0 * M_PI * uniform(gen);
        const double costheta = 1.0 - 2.0 * uniform(gen);
        const double sintheta = std::sqrt(1.0 - costheta * costheta);

        const double rx = rr * sintheta * std::cos(phi);
        const double ry = rr * sintheta * std::sin(phi);
        const double rz = rr * costheta;

        const double t = 0.0;
        const double px = 0.0;
        const double py = 0.0;
        const double pz = 0.0;
        const double p0 = std::sqrt(kProtonMass * kProtonMass + px * px + py * py + pz * pz);

        out << std::fixed << std::setprecision(2) << t << " "
            << std::fixed << std::setprecision(5)
            << rx << " "
            << ry << " "
            << rz << " "
            << kProtonMass << " "
            << p0 << " "
            << px << " "
            << py << " "
            << pz << " "
            << 2212 << " "
            << particle_index << " "
            << 1 << '\n';

        ++particle_index;
    }

    const int total_neutrons = A - Z + random_binary;
    for (int i = 0; i < total_neutrons; ++i) {
        const double rr = sample_point_in_sphere(RA, gen);
        const double phi = 2.0 * M_PI * uniform(gen);
        const double costheta = 1.0 - 2.0 * uniform(gen);
        const double sintheta = std::sqrt(1.0 - costheta * costheta);

        const double rx = rr * sintheta * std::cos(phi);
        const double ry = rr * sintheta * std::sin(phi);
        const double rz = rr * costheta;

        const double t = 0.0;
        const double px = 0.0;
        const double py = 0.0;
        const double pz = 0.0;
        const double p0 = std::sqrt(kNeutronMass * kNeutronMass + px * px + py * py + pz * pz);

        out << std::fixed << std::setprecision(2) << t << " "
            << std::fixed << std::setprecision(5)
            << rx << " "
            << ry << " "
            << rz << " "
            << kNeutronMass << " "
            << p0 << " "
            << px << " "
            << py << " "
            << pz << " "
            << 2112 << " "
            << particle_index << " "
            << 0 << '\n';

        ++particle_index;
    }
}

} // namespace

void write_event_headers(std::ostream& out) {
    out << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n";
    out << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none\n";
}

void write_event_output(
    int32_t event_number,
    int Z,
    int A,
    const DISKinematics& kin,
    const std::vector<Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out) {

    write_metadata_json(event_number, Z, A, kin, meta_out);

    // Count number of hadrons in the shower
    const int count = count_final_particles(particles);
    // Write the event header (OSCAR-style)
    event_out << "# event " << event_number << " out " << count << "\n";

    int32_t particle_index = 0;
    write_final_hadrons(particles, event_out, particle_index);
    write_spectator_nucleons(Z, A, event_out, particle_index);

    // Write event footer
    event_out << "# event " << event_number << " end 0\n";
}