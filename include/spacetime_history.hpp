#pragma once

#include "Pythia8/Event.h"
#include "HepMC3/GenEvent.h"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

/** Read-only diagnostic observations, serialized as HepMC3 Asciiv3 snapshots.
 * Observation vertices are not physical interactions. Exact PYTHIA ancestry
 * remains in attributes, so absent vertices and conflicting daughter positions
 * are never silently replaced by a reconstructed spacetime model.
 */
class SpacetimeHistory {
public:
    SpacetimeHistory(int64_t event_id, int64_t attempt);
    void capture(const Pythia8::Event& event, const std::string& stage,
                 const std::string& record, const std::string& vertex_convention,
                 const std::vector<int>& source_indices = {});
    void capture_particles(const std::vector<Pythia8::Particle>& particles,
                           const std::string& stage, const std::string& record,
                           const std::string& vertex_convention,
                           const std::vector<int>& indices = {},
                           const std::vector<int>& source_indices = {},
                           const std::string& comment = "");
    void set_context(const Pythia8::Vec4& hard_origin_gev_inv,
                     const Pythia8::Vec4& target_momentum,
                     bool source_parton_vertices_in_fm);
    void fragmentation(const std::string& trace);
    void output_particle(const Pythia8::Particle& particle, int output_id,
                         const Pythia8::Vec4& position_mm, bool spectator = false);
    void finish_output();
    /// Attach bookkeeping metadata to the completed OSCAR observation.
    void output_attribute(const std::string& name, const std::string& value);
    void write(std::ostream& out) const;

private:
    int64_t event_id_;
    int64_t attempt_;
    Pythia8::Vec4 origin_, target_;
    bool source_vertices_in_fm_ = false;
    int source_shower_size_ = 0;
    std::vector<HepMC3::GenEvent> snapshots_;
    std::vector<int> final_hadronizer_indices_;
    std::vector<Pythia8::Particle> output_particles_;
    std::vector<bool> output_spectators_;
    std::vector<double> output_formation_times_mm_;
    std::string fragmentation_trace_;
};
