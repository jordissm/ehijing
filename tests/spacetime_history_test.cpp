#include "spacetime_history.hpp"
#include "ehijing_constants.hpp"
#include "filesystem_utils.hpp"
#include "event_output.hpp"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void near(double actual, double expected) {
    require(std::abs(actual - expected) < 1e-25, "Vertex conversion mismatch");
}
HepMC3::GenEvent read(HepMC3::ReaderAscii& reader) {
    HepMC3::GenEvent result;
    reader.read_event(result);
    require(!reader.failed(), "HepMC3 reader rejected history");
    return result;
}
std::vector<HepMC3::ConstGenParticlePtr> observations(const HepMC3::GenEvent& event) {
    std::vector<HepMC3::ConstGenParticlePtr> result;
    for (const auto& p : event.particles())
        if (p->attribute_as_string("diagnostic_anchor") != "1") result.push_back(p);
    return result;
}
} // namespace

int main() {
    try {
        using Pythia8::Particle;
        using Pythia8::Vec4;
        // These fixtures exercise serialization without generating events or
        // running the physics generator.
        Particle missing(21, -51, 6, 0, 9, 10, 101, 102, 1., 2., 3., 4., 0.);
        missing.AddMultipleCollision({7., 3.}, {0.4, 0.2}, {1., 2.});
        Particle zero(1, 51, 6, 0, 0, 0, 103, 0, 1., 0., 0., 1., 0.);
        zero.vProd(Vec4());
        Particle hadron(211, 83, 1, 2, 0, 0, 0, 0, 0.1, 0.2, 0.3, 1., 0.139);
        hadron.vProd(Vec4(2e-12, 3e-12, 4e-12, 5e-12));
        hadron.productionVertices(Vec4(1e-12, 2e-12, 3e-12, 2e-12),
                                  Vec4(3e-12, 4e-12, 5e-12, 4e-12));
        hadron.tau(1.234567890123456e-12);
        SpacetimeHistory history(1234567890123LL, 2);
        history.set_context(Vec4(1., 2., 3., 0.), Vec4(0., 0., 0., 0.938), false);
        history.capture_particles({missing, zero}, "high_scale_shower", "shower",
                                  "source_parton_vertices", {6, 9});
        history.capture_particles({hadron}, "hadronizer_import", "hadronizer",
                                  "relative_mm_hard_interaction", {1}, {9});
        history.capture_particles({hadron}, "hadronization_and_decays", "hadronizer",
                                  "absolute_mm_saved_vertex", {3});
        history.fragmentation("# diagnostic comment\nB0=(0, 1, 2, 3)\nbackslash \\ preserved");
        history.output_particle(hadron, 0, Vec4(1e-12, 2e-12, 3e-12, 0.));
        history.output_particle(Particle(2112, 1, 0, 0, 0, 0, 0, 0,
                                        0., 0., 0., 0.939, 0.939), 1,
                                Vec4(-1e-12, 0., 0., 0.), true);
        history.finish_output();
        std::ostringstream output;
        history.write(output);
        require(output.str().find("HepMC::Asciiv3-START_EVENT_LISTING") != std::string::npos,
                "Not native HepMC3 ASCII");
        require(output.str().find("HepMC::Asciiv3-END_EVENT_LISTING") != std::string::npos,
                "Missing HepMC footer");
        std::istringstream input(output.str());
        HepMC3::ReaderAscii reader(input);
        auto shower = read(reader);
        require(shower.event_number() == 0, "Snapshot order lost");
        require(shower.attribute_as_string("event_id") == "1234567890123", "64-bit ID truncated");
        require(observations(shower).size() == 2, "Particle lost");
        auto p = observations(shower)[0];
        require(p->attribute_as_string("record_index") == "6", "Source identity lost");
        require(p->attribute_as_string("pythia_status") == "-51", "Native status lost");
        require(p->attribute_as_string("mother1") == "6", "Ancestry lost");
        require(p->attribute_as_string("has_vertex") == "0", "Missing vertex misrepresented");
        require(observations(shower)[1]->attribute_as_string("has_vertex") == "1",
                "Known zero vertex misrepresented");
        require(p->attribute_as_string("collisions").find("0 7 ") != std::string::npos,
                "Sampled collision order lost");
        require(p->production_vertex()->status() == 900, "Observation vertex marker lost");
        auto imported = read(reader);
        require(observations(imported)[0]->attribute_as_string("source_index") == "9",
                "Import mapping lost");
        auto fragmented = read(reader);
        near(observations(fragmented)[0]->production_vertex()->position().t(), 5e-12);
        require(observations(fragmented)[0]->attribute_as_string("oscar_id") == "0",
                "Final hadron mapping lost");
        require(!observations(fragmented)[0]->attribute_as_string("breakup_first_xyzt_mm").empty(),
                "Breakup references lost");
        require(!observations(fragmented)[0]->attribute_as_string("native_decay_xyzt_mm").empty(),
                "Native decay endpoint lost");
        auto oscar = read(reader);
        require(observations(oscar).size() == 2, "Spectator omitted");
        require(observations(oscar)[0]->attribute_as_string("source_index") == "3",
                "OSCAR to hadronizer mapping lost");
        require(observations(oscar)[1]->attribute_as_string("spectator") == "1",
                "Spectator provenance lost");
        require(observations(oscar)[0]->attribute_as_string("oscar_form_time_fm") == "5",
                "OSCAR formation time lost during extrapolation");
        near(observations(oscar)[0]->production_vertex()->position().t(), 0.);
        require(oscar.attribute_as_string("comment_fragmentation").find("\nB0=") != std::string::npos,
                "Multiline comments did not round trip");
        require(missing.coll_ts()[0] == 7. && !missing.hasVertex(), "Observer mutated source");
        near(hadron.tProd(), 5e-12);

        SpacetimeHistory fm_history(0, 1);
        fm_history.set_context(Vec4(), Vec4(), true);
        Particle fm_parton = zero;
        fm_parton.vProd(Vec4(1., 2., 3., 4.));
        fm_history.capture_particles({fm_parton}, "high_scale_shower", "shower",
                                     "source_parton_vertices");
        Particle mm_parton = hadron;
        fm_history.capture_particles({fm_parton, mm_parton}, "medium_modified", "shower",
                                     "relative_mm_hard_interaction");
        std::ostringstream fm_output;
        fm_history.write(fm_output);
        std::istringstream fm_input(fm_output.str());
        HepMC3::ReaderAscii fm_reader(fm_input);
        auto converted = read(fm_reader);
        near(observations(converted)[0]->production_vertex()->position().t(), 4e-12);
        require(observations(converted)[0]->attribute_as_string("stored_vertex_units") == "fm",
                "Source unit exception lost");
        auto mixed = read(fm_reader);
        near(observations(mixed)[0]->production_vertex()->position().t(), 4e-12);
        near(observations(mixed)[1]->production_vertex()->position().t(), 5e-12);
        require(observations(mixed)[1]->attribute_as_string("stored_vertex_units") == "mm",
                "Appended medium units were mistaken for source units");
        bool rejected = false;
        try { fm_history.capture_particles({zero}, "test", "shower", "mm", {1, 2}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Invalid index map accepted");
        require(zero_pad_int(42) == "00000042", "Output numbering changed");
        const auto paths = make_event_paths(std::filesystem::temp_directory_path() /
                                           "electra-history-path-test", 42, 10);
        require(paths.historyPath.filename() == "event_00000042.history", "Wrong history filename");
        require(paths.historyPath.parent_path() == paths.eventPath.parent_path(), "Wrong history shard");

        // A=Z=1 with a struck proton has no spectators, so both calls consume
        // no output RNG draws. Verify that diagnostics leave existing files exact.
        DISKinematics kin{};
        std::ostringstream plain_oscar, plain_meta, traced_oscar, traced_meta;
        write_event_output(42, 1, 1, StruckNucleon::Proton, kin, {hadron},
                           plain_oscar, plain_meta);
        SpacetimeHistory output_history(42, 1);
        output_history.capture_particles({hadron}, "hadronization_and_decays",
                                         "hadronizer", "absolute_mm_saved_vertex", {3});
        write_event_output(42, 1, 1, StruckNucleon::Proton, kin, {hadron},
                           traced_oscar, traced_meta, &output_history);
        require(plain_oscar.str() == traced_oscar.str(), "Diagnostics changed OSCAR output");
        require(plain_meta.str() == traced_meta.str(), "Diagnostics changed metadata output");
        near(hadron.tProd(), 5e-12);
        std::cout << "Spacetime history round-trip checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
