// Optional hadronization spacetime diagnostics for the eHIJING PYTHIA fork.
// Coordinates are printed as (t, x, y, z), in fm/c and fm respectively.
// Callers supply an output stream; diagnostics never fall back to stderr.
#ifndef Pythia8_HadronVertexTrace_H
#define Pythia8_HadronVertexTrace_H

#include "Pythia8/Event.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/ColourConnectionHistory.h"
#include <algorithm>
#include <iomanip>
#include <ostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Pythia8 {
namespace HadronVertexTrace {

struct Tag {
  int pass, string;
  std::string branch;
};

struct Snapshot {
  Vec4 vertex;
  bool hasPair;
  Vec4 first, second;
};

struct Context {
  bool active = false;
  std::ostream* output = nullptr;
  long long event = -1;
  int attempt = 0, pass = 0, string = -1;
  std::map<int, Tag> tags;
  std::map<int, Snapshot> snapshots;
  const ColourConnectionHistory* colourHistory = nullptr;
};

// One context per thread, shared by the wrapper and fragmentation routines.
inline Context& context() {
  static thread_local Context value;
  return value;
}

inline void prefix(std::ostream& out) {
  const Context& c = context();
  out << std::defaultfloat << std::setprecision(12) << std::boolalpha
      << "[pass=" << c.pass << " string=" << c.string << "] ";
}

inline void write(const std::string& text) {
  if (context().active && context().output) *context().output << text;
}

inline void vertex(std::ostream& out, const char* label, const Vec4& v,
  double factor = MM2FM) {
  out << label << "=(" << v.e() * factor << ", " << v.px() * factor
      << ", " << v.py() * factor << ", " << v.pz() * factor << ")";
}

inline void momentum(std::ostream& out, const Vec4& p) {
  out << "p(E,px,py,pz)_GeV=(" << p.e() << ", " << p.px() << ", "
      << p.py() << ", " << p.pz() << ")";
}

inline bool sameVertex(const Vec4& a, const Vec4& b) {
  return a.e() == b.e() && a.px() == b.px() && a.py() == b.py()
      && a.pz() == b.pz();
}

inline void beginEvent(long long event, int attempt, std::ostream& output) {
  Context& c = context();
  c = Context();
  c.active = true;
  c.output = &output;
  c.event = event;
  c.attempt = attempt;
  std::ostringstream out;
  out << "# eHIJING fragmentation trace v5\n"
      << "Event: " << event << "   Attempt: " << attempt << '\n'
      << "Coordinates: (t, x, y, z); time in fm/c; position in fm; momentum in GeV.\n"
      << "Breakup B labels are local to each string leg.\n"
      << "For an open reference string, B0 and the last B are endpoint turns, not pair creations.\n"
      << "Absolute breakup points include the inherited origin.\n"
      << "saved_vertex is the formation point in HadronVertex:mode = 2.\n"
      << "output_ID matches the OSCAR hadron row; spectators are not traced.\n"
      << "An unavailable production time is not a measured zero.\n\n";
  out << "GEOMETRY records describe the longitudinal massless-reference Lund worldsheet.\n"
      << "They exclude transverse tunnelling/smearing and heavy-endpoint corrections.\n"
      << "Reference constituent meeting points can differ from saved_vertex.\n"
      << "REFERENCE_BREAKUP points include spacetime conversion corrections, before heavy offsets and smearing.\n"
      << "STRING partons give input endpoint states, not sampled constituent histories.\n\n";
  out << "Colour connection times describe the recorded medium-stage topology, not hadron production/formation.\n"
      << "Connection time dates the start of the current uninterrupted complete source component; first_time also records its first appearance.\n"
      << "Untimed pre-medium/shower connections remain unavailable. Scalar connection clocks cannot be boosted without an interaction four-vertex.\n\n";
  write(out.str());
}

class EventScope {
public:
  EventScope(bool enabledIn, long long event, int attempt,
    std::ostream* output = nullptr)
    : enabled(enabledIn && output != nullptr) {
    if (enabled) beginEvent(event, attempt, *output);
  }
  ~EventScope() { finish(); }
  EventScope(const EventScope&) = delete;
  EventScope& operator=(const EventScope&) = delete;
  void finish() {
    if (!enabled) return;
    write("\nEnd of fragmentation trace.\n");
    context() = Context();
    enabled = false;
  }
private:
  bool enabled;
};

// Direct PYTHIA callers can supply their own stream through EventScope.
class PassScope {
public:
  PassScope(bool enabledIn, const Event& event)
    : enabled(enabledIn && context().active) {
    if (!enabled) return;
    // A failed attempt may have been rolled back by PYTHIA. Drop tags for
    // entries no longer present before their event indices can be reused.
    context().tags.erase(context().tags.lower_bound(event.size()),
      context().tags.end());
    context().snapshots.erase(context().snapshots.lower_bound(event.size()),
      context().snapshots.end());
    ++context().pass;
    context().string = -1;
  }
private:
  bool enabled;
};

inline void colourConnection(std::ostream& out, const Event& event,
  const std::vector<int>& partons) {
  if (!context().colourHistory) {
    out << "color_connection_status=unavailable_no_history color_connection_time=unavailable\n";
    return;
  }
  const auto connection = context().colourHistory->connection(event, partons);
  out << "color_connection_status=" << connection.status << " color_connection_time=";
  if (connection.timed) out << connection.time; else out << "unavailable";
  out << " color_connection_first_status=" << connection.firstStatus
      << " color_connection_first_time=";
  if (connection.firstTimed) out << connection.firstTime; else out << "unavailable";
  out << " color_latest_known_link_time=";
  if (connection.hasLatestLink) out << connection.latestLinkTime; else out << "unavailable";
  out << " color_source_ids=(";
  bool first = true;
  for (int id : connection.sources) { if (!first) out << ','; out << id; first = false; }
  out << ") color_link_ids=(";
  first = true;
  for (int id : connection.links) { if (!first) out << ','; out << id; first = false; }
  out << ")\n";
}

inline void beginString(int index, const ColSinglet& system,
  const Event& event, const char* branch) {
  if (!context().active) return;
  context().string = index;
  std::ostringstream out;
  out << "\n--- String " << index << " / pass " << context().pass << " ---\n";
  prefix(out);
  out << "STRING branch=" << branch << " topology="
      << (system.hasJunction ? "junction" : system.isClosed ? "closed" : "open")
      << " mass_GeV=" << system.mass << " mass_excess_GeV="
      << system.massExcess << ' ';
  momentum(out, system.pSum);
  out << "\n  ordered_partons:";
  for (int i : system.iParton) out << ' ' << i;
  out << '\n';
  for (int i : system.iParton) {
    if (i < 0) {
      out << "  junction_marker=" << i << '\n';
      continue;
    }
    const Particle& p = event[i];
    out << "  parton index=" << i << " pdg=" << p.id()
        << " name=\"" << p.name() << '"'
        << " mass_GeV=" << p.m()
        << " role=" << (system.hasJunction ? "junction_parton"
          : system.isClosed ? "gluon_loop"
          : i == system.iParton.front() ? "positive_endpoint"
          : i == system.iParton.back() ? "negative_endpoint" : "gluon_kink")
        << " status=" << p.status() << " colour_source_id=" << p.colourHistoryId()
        << " colors=" << p.col() << ','
        << p.acol() << " mothers=" << p.mother1() << ',' << p.mother2()
        << " hasVertex=" << p.hasVertex() << ' ';
    momentum(out, p.p());
    out << ' ';
    vertex(out, "absolute_vertex", p.vProd());
    if (context().colourHistory) {
      std::set<int> sources;
      const bool complete = context().colourHistory->sourceIds(event, i, sources);
      out << " color_provenance_complete=" << complete << " color_source_ids=(";
      bool first = true;
      for (int id : sources) { if (!first) out << ','; out << id; first = false; }
      out << ')';
    }
    out << '\n';
  }
  prefix(out);
  out << "STRING_CONNECTION ";
  colourConnection(out, event, system.iParton);
  write(out.str());
}

inline void recordHadron(const Event& event, int index, const char* branch) {
  if (!context().active) return;
  Context& c = context();
  c.tags[index] = Tag{c.pass, c.string, branch};
  const Particle& p = event[index];
  c.snapshots[index] = Snapshot{p.vProd(), p.hasProductionTime(),
    p.vBreakupFirst(), p.vBreakupSecond()};
}

inline void breakups(const char* leg, const std::vector<Vec4>& points) {
  if (!context().active) return;
  std::ostringstream out;
  out << "\n  Breakup chain (relative to the inherited origin):\n";
  for (int i = 0; i < int(points.size()); ++i) {
    prefix(out);
    out << "BREAKUP leg=" << leg << " B" << i << ' ';
    vertex(out, "relative_vertex", points[i], 1.);
    out << '\n';
  }
  write(out.str());
}

// Copy the parameter-plane geometry used by PYTHIA. Filling unused regions in
// this diagnostic copy does not alter fragmentation, particle data or RNG state.
// Null characteristics of this worldsheet give the classical longitudinal
// constituent paths, including changes of direction at gluon-region boundaries.
inline void geometry(const char* leg, const StringSystem& system,
  std::vector<int> partons, Event& event,
  const std::vector<StringVertex>& breaks, double kappa,
  const RotBstMatrix* toJRF = nullptr, const RotBstMatrix* fromJRF = nullptr,
  const std::vector<Vec4>* effectivePoints = nullptr,
  const Vec4* origin = nullptr) {
  if (!context().active) return;
  StringSystem reference = system;
  std::ostringstream out;
  prefix(out);
  out << "GEOMETRY leg=" << leg << " max_region=" << reference.iMax
      << " kappa_GeV_per_fm=" << kappa
      << " model=longitudinal_massless_reference";
  if (origin) { out << ' '; vertex(out, "origin_fm", *origin); }
  out << '\n';
  for (int pos = 0; pos <= reference.iMax; ++pos) {
    for (int neg = 0; neg <= reference.iMax - pos; ++neg) {
      StringRegion region = reference.region(pos, neg);
      if (!region.isSetUp) region.setUp(reference.regionLowPos(pos).pPos,
        reference.regionLowNeg(neg).pNeg, true);
      prefix(out);
      out << "REGION leg=" << leg << " i_pos=" << pos << " i_neg=" << neg
          << " available=" << !region.isEmpty;
      if (!region.isEmpty) {
        // Both offset and basis are recorded in the original event frame.
        Vec4 offset = toJRF
          ? region.gluonOffsetJRF(partons, event, pos, neg, *toJRF)
          : region.gluonOffset(partons, event, pos, neg);
        offset /= kappa;
        Vec4 plus = region.pPos, minus = region.pNeg;
        if (fromJRF) {
          offset.rotbst(*fromJRF);
          plus.rotbst(*fromJRF);
          minus.rotbst(*fromJRF);
        }
        out << ' '; vertex(out, "offset_fm", offset, 1.);
        out << ' '; vertex(out, "plus_GeV", plus, 1.);
        out << ' '; vertex(out, "minus_GeV", minus, 1.);
      }
      out << '\n';
    }
  }
  for (int i = 0; i < int(breaks.size()); ++i) {
    const StringVertex& b = breaks[i];
    prefix(out);
    out << "PARAMETER_BREAKUP leg=" << leg << " index=" << i
        << " i_pos=" << b.iRegPos << " i_neg=" << b.iRegNeg
        << " x_pos=" << b.xRegPos << " x_neg=" << b.xRegNeg
        << " sampled_from=" << (b.fromPos ? "positive" : "negative")
        << " role=" << ((i == 0 || i == int(breaks.size()) - 1)
          ? "boundary" : "pair_creation") << '\n';
    if (effectivePoints && i < int(effectivePoints->size())) {
      Vec4 position = (*effectivePoints)[i];
      if (fromJRF) position.rotbst(*fromJRF);
      prefix(out);
      out << "REFERENCE_BREAKUP leg=" << leg << " index=" << i << ' ';
      vertex(out, "position_fm", position, 1.);
      out << '\n';
    }
  }
  write(out.str());
}

inline void hadron(const Event& event, int index, const char* branch,
  int firstIndex, int secondIndex, const Vec4& origin,
  const Vec4& first, const Vec4& second, const Vec4& displacement,
  double redOsc, double kappa) {
  if (!context().active) return;
  recordHadron(event, index, branch);
  const Particle& p = event[index];
  Vec4 a = origin * MM2FM + first;
  Vec4 b = origin * MM2FM + second;
  double earliest = std::min(a.e(), b.e());
  std::ostringstream out;
  out << '\n'; prefix(out);
  out << "HADRON index=" << index << " pdg=" << p.id()
      << " name=\"" << p.name() << '"'
      << " status=" << p.status() << " branch=" << branch
      << " breaks=B" << firstIndex << ",B" << secondIndex
      << " redOsc=" << redOsc << " kappa_GeV_per_fm=" << kappa << '\n';
  out << "  fragmentation_side="
      << ((p.statusAbs() == 83 || p.statusAbs() == 87) ? "positive"
        : (p.statusAbs() == 84 || p.statusAbs() == 88) ? "negative" : "junction/mini")
      << ' '; momentum(out, p.p());
  out << '\n';
  out << "  "; vertex(out, "inherited_origin", origin);
  out << "\n  "; vertex(out, "absolute_breakup_first", a, 1.);
  out << "\n  "; vertex(out, "absolute_breakup_second", b, 1.);
  out << "\n  "; vertex(out, "prescription_displacement", displacement, 1.);
  out << "\n  "; vertex(out, "saved_vertex", p.vProd());
  out << "\n  computed_earliest_breakup=" << earliest
      << " saved_minus_earliest=" << p.tProd() * MM2FM - earliest
      << " hasProductionTime=" << p.hasProductionTime();
  if (p.hasProductionTime()) out << " stored_prod_time=" << p.tProduction() * MM2FM;
  else out << " stored_prod_time=unavailable(no_stored_breakup_pair)";
  out << '\n';
  write(out.str());
}

inline void collapse(const Event& event, int index, int recoilString,
  const Vec4& origin, const std::vector<int>& partons) {
  if (!context().active) return;
  recordHadron(event, index, "one-hadron-collapse");
  const Particle& p = event[index];
  std::ostringstream out;
  out << '\n'; prefix(out);
  out << "HADRON index=" << index << " pdg=" << p.id()
      << " name=\"" << p.name() << '"'
      << " status=" << p.status()
      << " branch=one-hadron-collapse recoil_string=" << recoilString
      << " production_time=unavailable(no_reconstructed_breakup_pair)\n  ";
  momentum(out, p.p()); out << "\n  ";
  vertex(out, "inherited_origin", origin);
  out << "\n  "; vertex(out, "prescription_displacement", p.vProd() - origin);
  out << "\n  "; vertex(out, "saved_vertex", p.vProd());
  out << "\n  "; colourConnection(out, event, partons);
  write(out.str());
}

// Follow decay/copy ancestry for final hadrons created after fragmentation.
// Multiple ancestor strings are retained instead of choosing one arbitrarily.
inline void ancestorStrings(std::ostream& out, const Event& event, int index) {
  std::vector<int> pending(1, index);
  std::set<int> visited;
  std::map<std::pair<int, int>, std::string> sources;
  while (!pending.empty()) {
    int i = pending.back();
    pending.pop_back();
    if (i <= 0 || i >= event.size() || !visited.insert(i).second) continue;
    auto tag = context().tags.find(i);
    if (tag != context().tags.end()) {
      sources[std::make_pair(tag->second.pass, tag->second.string)]
        = tag->second.branch;
    } else {
      std::vector<int> mothers = event[i].motherList();
      pending.insert(pending.end(), mothers.begin(), mothers.end());
    }
  }
  out << " ancestor_strings=[";
  for (const auto& source : sources)
    out << " (pass=" << source.first.first << ",string=" << source.first.second
        << ",branch=" << source.second << ')';
  out << " ]";
}

inline void dumpFinal(const Event& event, bool success) {
  if (!context().active) return;
  context().string = -1;
  std::ostringstream out;
  out << std::defaultfloat << std::setprecision(12);
  out << "\n=== Exported hadrons: time summary [fm/c] ===\n"
      << "Saved t is the formation time in mode 2. Missing production times are marked unavailable.\n"
      << std::left << std::setw(8) << "ID" << std::setw(9) << "PYTHIA"
      << std::setw(12) << "PDG" << std::setw(14) << "Pass/string"
      << std::setw(22) << "prod_time" << std::setw(22) << "saved t" << "Particle\n";
  int summaryIndex = 0;
  for (int i = 0; success && i < event.size(); ++i) {
    const Particle& p = event[i];
    if (!p.isFinal()) continue;
    const int id = summaryIndex++;
    if (!p.isHadron()) continue;
    auto tag = context().tags.find(i);
    std::string source = tag == context().tags.end() ? "see ancestry"
      : std::to_string(tag->second.pass) + "/" + std::to_string(tag->second.string);
    out << std::setw(8) << id << std::setw(9) << i << std::setw(12) << p.id()
        << std::setw(14) << source;
    if (p.hasProductionTime()) out << std::setw(22) << p.tProduction() * MM2FM;
    else out << std::setw(22) << "unavailable";
    out << std::setw(22) << p.tProd() * MM2FM << p.name() << '\n';
  }
  out << "\n=== Final PYTHIA particle record ===\n";
  prefix(out);
  out << "FINAL_EVENT success=" << success << " size=" << event.size() << '\n';
  int outputIndex = 0;
  for (int i = 0; i < event.size(); ++i) {
    const Particle& p = event[i];
    out << "\n  PARTICLE index=" << i << " pdg=" << p.id()
        << " name=\"" << p.name() << '"' << " status=" << p.status()
        << " final=" << p.isFinal() << " hadron=" << p.isHadron()
        << "\n    mothers=" << p.mother1() << ',' << p.mother2()
        << " daughters=" << p.daughter1() << ',' << p.daughter2();
    if (success && p.isFinal()) out << " output_ID=" << outputIndex++;
    auto tag = context().tags.find(i);
    if (tag != context().tags.end())
      out << " source_pass=" << tag->second.pass
          << " source_string=" << tag->second.string
          << " branch=" << tag->second.branch;
    else if (p.isFinal() && p.isHadron()) ancestorStrings(out, event, i);
    out << "\n    colors=" << p.col() << ',' << p.acol() << ' ';
    momentum(out, p.p());
    out << " mass_GeV=" << p.m() << " hasVertex=" << p.hasVertex() << "\n    ";
    vertex(out, "saved_vertex", p.vProd());
    out << "\n    raw_saved_time_mm_per_c=" << p.tProd()
        << " hasProductionTime=" << p.hasProductionTime();
    if (p.hasProductionTime()) {
      out << "\n    prod_time=" << p.tProduction() * MM2FM
          << " raw_prod_time_mm_per_c=" << p.tProduction() << ' ';
      out << "\n    ";
      vertex(out, "stored_breakup_first", p.vBreakupFirst());
      out << "\n    "; vertex(out, "stored_breakup_second", p.vBreakupSecond());
    } else out << "\n    prod_time=unavailable(no_stored_breakup_pair)";
    auto snapshot = context().snapshots.find(i);
    if (snapshot != context().snapshots.end()) {
      const Snapshot& s = snapshot->second;
      out << "\n    vertex_changed_after_fragmentation=" << !sameVertex(s.vertex, p.vProd())
          << " breakup_metadata_changed_after_fragmentation="
          << (s.hasPair != p.hasProductionTime()
              || !sameVertex(s.first, p.vBreakupFirst())
              || !sameVertex(s.second, p.vBreakupSecond()));
    }
    out << '\n';
  }
  write(out.str());
}

} // namespace HadronVertexTrace
} // namespace Pythia8
#endif
