#include "spacetime_history.hpp"
#include "ehijing_constants.hpp"

#include "HepMC3/Attribute.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
namespace units = ehijing::constants::units;

template <typename Object>
void attribute(Object& object, const std::string& name, const std::string& value) {
    object.add_attribute(name, std::make_shared<HepMC3::StringAttribute>(value));
}

std::string vector_text(const Pythia8::Vec4& value) {
    std::ostringstream out;
    out << std::setprecision(17) << value.px() << ' ' << value.py() << ' '
        << value.pz() << ' ' << value.e();
    return out.str();
}

std::string scalar_text(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::string collisions(const Pythia8::Particle& particle) {
    std::ostringstream out;
    out << std::setprecision(17);
    const auto& times = particle.coll_ts();
    const auto& kicks = particle.coll_qt2s();
    const auto& angles = particle.coll_phis();
    // Preserve sampler order, including any inconsistent vector sizes.
    out << "counts(t,qt2,phi)=" << times.size() << ',' << kicks.size() << ','
        << angles.size() << "; sampler_order: index t_GeV^-1 qt2_GeV^2 phi_rad";
    const auto count = std::max(times.size(), std::max(kicks.size(), angles.size()));
    for (size_t i = 0; i < count; ++i) {
        out << '\n' << i;
        if (i < times.size()) out << ' ' << times[i]; else out << " unavailable";
        if (i < kicks.size()) out << ' ' << kicks[i]; else out << " unavailable";
        if (i < angles.size()) out << ' ' << angles[i]; else out << " unavailable";
    }
    return out.str();
}
} // namespace

SpacetimeHistory::SpacetimeHistory(int64_t event_id, int attempt)
    : event_id_(event_id), attempt_(attempt) {}

void SpacetimeHistory::set_context(const Pythia8::Vec4& origin,
                                  const Pythia8::Vec4& target,
                                  bool source_vertices_in_fm) {
    origin_ = origin;
    target_ = target;
    source_vertices_in_fm_ = source_vertices_in_fm;
}

void SpacetimeHistory::capture(const Pythia8::Event& event,
                               const std::string& stage, const std::string& record,
                               const std::string& convention,
                               const std::vector<int>& sources) {
    std::vector<Pythia8::Particle> particles;
    particles.reserve(event.size());
    for (int i = 0; i < event.size(); ++i) particles.push_back(event[i]);
    capture_particles(particles, stage, record, convention, {}, sources);
}

void SpacetimeHistory::capture_particles(
    const std::vector<Pythia8::Particle>& particles, const std::string& stage,
    const std::string& record, const std::string& convention,
    const std::vector<int>& indices, const std::vector<int>& sources,
    const std::string& comment) {
    if ((!indices.empty() && indices.size() != particles.size()) ||
        (!sources.empty() && sources.size() != particles.size()))
        throw std::invalid_argument("History index map does not match snapshot");

    HepMC3::GenEvent snapshot(HepMC3::Units::GEV, HepMC3::Units::MM);
    snapshot.set_event_number(static_cast<int>(snapshots_.size()));
    attribute(snapshot, "event_id", std::to_string(event_id_));
    attribute(snapshot, "attempt", std::to_string(attempt_));
    attribute(snapshot, "history_schema", "electra-spacetime-v1");
    attribute(snapshot, "stage", stage);
    attribute(snapshot, "record", record);
    attribute(snapshot, "vertex_convention", convention);
    attribute(snapshot, "momentum_frame", "generator_event_frame; no diagnostic boosts applied");
    attribute(snapshot, "hard_origin_xyzt_GeV_inv", vector_text(origin_));
    attribute(snapshot, "target_px_py_pz_E_GeV", vector_text(target_));
    attribute(snapshot, "source_parton_vertex_setting", source_vertices_in_fm_ ? "on" : "off");
    if (source_vertices_in_fm_)
        attribute(snapshot, "comment_source_units",
            "Source PartonVertex is enabled. Original source parton vertices use fm; "
            "medium-appended displacements use mm. The operative hadronizer "
            "import assumes mm and requires review for this unsupported configuration.");
    attribute(snapshot, "comment",
        "Diagnostic observation graph, not a physical interaction graph. E numbers "
        "order snapshots, not physical time. Identity is (record,record_index). "
        "PYTHIA mother/daughter fields retain native range/copy semantics. "
        "Vertex status 900 denotes an observation. A zero position with "
        "has_vertex=0 is unavailable, not a measured origin. No shower delays, "
        "collision positions or Lorentz transformations are reconstructed. " + comment);
    attribute(snapshot, "comment_anchors",
        "Particles with diagnostic_anchor=1 (PDG 0, status 4) are serialization "
        "anchors, not simulated particles. They give observation vertices an "
        "incoming edge for compatibility with HepMC3 readers rejecting empty "
        "incoming lists. Filter them when counting simulated particles.");

    if (stage == "hadronization_and_decays") final_hadronizer_indices_.clear();
    if (stage == "high_scale_shower") {
        source_shower_size_ = 0;
        for (size_t i = 0; i < particles.size(); ++i)
            source_shower_size_ = std::max(source_shower_size_,
                (indices.empty() ? static_cast<int>(i) : indices[i]) + 1);
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const int index = indices.empty() ? static_cast<int>(i) : indices[i];
        auto particle = std::make_shared<HepMC3::GenParticle>(
            HepMC3::FourVector(p.px(), p.py(), p.pz(), p.e()), p.id(),
            p.isFinal() ? 1 : 2);
        particle->set_generated_mass(p.m());
        // Each stored position is preserved independently. Joining siblings
        // at one HepMC vertex would erase conflicting per-particle positions.
        double factor = 1.0;
        if (source_vertices_in_fm_ &&
            (record == "process" || (record == "shower" && index < source_shower_size_)))
            factor = units::fm_to_mm;
        const auto v = p.vProd() * factor;
        auto vertex = std::make_shared<HepMC3::GenVertex>(
            HepMC3::FourVector(v.px(), v.py(), v.pz(), v.e()));
        vertex->set_status(900);
        auto anchor = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(), 0, 4);
        snapshot.add_particle(anchor);
        attribute(*anchor, "diagnostic_anchor", "1");
        vertex->add_particle_in(anchor);
        vertex->add_particle_out(particle);
        snapshot.add_vertex(vertex);
        attribute(*particle, "record_index", std::to_string(index));
        attribute(*particle, "pythia_status", std::to_string(p.status()));
        attribute(*particle, "mother1", std::to_string(p.mother1()));
        attribute(*particle, "mother2", std::to_string(p.mother2()));
        attribute(*particle, "daughter1", std::to_string(p.daughter1()));
        attribute(*particle, "daughter2", std::to_string(p.daughter2()));
        attribute(*particle, "colour", std::to_string(p.col()));
        attribute(*particle, "anticolour", std::to_string(p.acol()));
        attribute(*particle, "has_vertex", p.hasVertex() ? "1" : "0");
        attribute(*particle, "stored_vertex_xyzt", vector_text(p.vProd()));
        attribute(*particle, "stored_vertex_units", factor == 1.0 ? "mm" : "fm");
        attribute(*particle, "vertex_convention", record == "shower" &&
            index < source_shower_size_ ? "source_parton_vertices" : convention);
        attribute(*particle, "tau_mm", scalar_text(p.tau()));
        attribute(*particle, "scale_GeV", scalar_text(p.scale()));
        attribute(*particle, "bookkeeping_system_entry", index == 0 &&
                  (record == "process" || record == "shower" || record == "hadronizer") ? "1" : "0");
        attribute(*particle, "geometry_origin_xyzt_GeV_inv",
                  vector_text(Pythia8::Vec4(p.Rx(), p.Ry(), p.Rz(), 0.0)));
        attribute(*particle, "collisions", collisions(p));
        attribute(*particle, "collisions_sampled", p.collset() ? "1" : "0");
        attribute(*particle, "has_production_references", p.hasProductionTime() ? "1" : "0");
        attribute(*particle, "has_breakup_pair", p.hasBreakupPair() ? "1" : "0");
        if (record == "hadronizer" && p.hasVertex() && p.tau() > 0. && p.m() > 0.)
            attribute(*particle, "native_decay_xyzt_mm", vector_text(p.vDec()));
        if (!sources.empty()) {
            attribute(*particle, "source_record", "shower");
            attribute(*particle, "source_index", std::to_string(sources[i]));
        }
        if (p.hasProductionTime()) {
            attribute(*particle, "production_first_xyzt_mm", vector_text(p.vProductionFirst()));
            attribute(*particle, "production_second_xyzt_mm", vector_text(p.vProductionSecond()));
            attribute(*particle, "production_origin_mode", std::to_string(p.productionOriginMode()));
        }
        if (p.hasBreakupPair()) {
            attribute(*particle, "breakup_first_xyzt_mm", vector_text(p.vBreakupFirst()));
            attribute(*particle, "breakup_second_xyzt_mm", vector_text(p.vBreakupSecond()));
        }
        if (stage == "hadronization_and_decays" && p.isFinal()) {
            attribute(*particle, "oscar_id", std::to_string(final_hadronizer_indices_.size()));
            final_hadronizer_indices_.push_back(index);
        }
    }
    snapshots_.push_back(std::move(snapshot));
}

void SpacetimeHistory::fragmentation(const std::string& trace) {
    fragmentation_trace_ = trace;
}

void SpacetimeHistory::output_particle(const Pythia8::Particle& particle,
                                      int output_id, const Pythia8::Vec4& position,
                                      bool spectator) {
    if (output_id != static_cast<int>(output_particles_.size()))
        throw std::invalid_argument("Non-contiguous OSCAR history ID");
    auto copy = particle;
    copy.vProd(position);
    output_particles_.push_back(copy);
    output_spectators_.push_back(spectator);
    output_formation_times_mm_.push_back(spectator ? 0.0 : particle.tProd());
}

void SpacetimeHistory::finish_output() {
    capture_particles(output_particles_, "oscar_output", "oscar",
                      "absolute_mm_t0_back_extrapolated",
                      {}, {}, "Positions precede OSCAR rounding to 5 decimals in fm. "
                      "Spectators are independently sampled by the existing output model.");
    auto& snapshot = snapshots_.back();
    for (const auto& particle : snapshot.particles()) {
        auto& p = *particle;
        const auto index = p.attribute_as_string("record_index");
        if (index.empty()) continue; // Serialization anchor.
        const auto i = static_cast<size_t>(std::stoul(index));
        attribute(p, "oscar_id", std::to_string(i));
        attribute(p, "spectator", output_spectators_[i] ? "1" : "0");
        attribute(p, "oscar_form_time_fm", scalar_text(output_formation_times_mm_[i] * units::mm_to_fm));
        attribute(p, "oscar_prod_time_fm", scalar_text(output_particles_[i].tProduction() * units::mm_to_fm));
        attribute(p, "oscar_xsecfac", output_spectators_[i] ? "1" : "0");
        if (!output_spectators_[i] && i < final_hadronizer_indices_.size()) {
            attribute(p, "source_record", "hadronizer");
            attribute(p, "source_index", std::to_string(final_hadronizer_indices_[i]));
        }
    }
}

void SpacetimeHistory::write(std::ostream& out) const {
    if (!out) throw std::runtime_error("Cannot write spacetime history");
    // Some HepMC3 versions close a supplied ofstream in WriterAscii::close().
    // Give the writer a plain ostream over the borrowed buffer so it can
    // finish the HepMC footer without closing the caller-owned file.
    std::ostream stream(out.rdbuf());
    HepMC3::WriterAscii writer(stream);
    writer.set_precision(17);
    for (size_t i = 0; i < snapshots_.size(); ++i) {
        // Serialize a diagnostic copy: neither the source generator nor the
        // captured state is exposed to the writer.
        auto snapshot = snapshots_[i];
        if (i + 1 == snapshots_.size() && !fragmentation_trace_.empty())
            attribute(snapshot, "comment_fragmentation", fragmentation_trace_);
        writer.write_event(snapshot);
        if (writer.failed()) throw std::runtime_error("Cannot write spacetime history");
    }
    writer.close();
    out.flush();
    if (writer.failed() || !stream || !out)
        throw std::runtime_error("Cannot finish spacetime history");
}
