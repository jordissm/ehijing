// Optional hadronization spacetime diagnostics for the eHIJING PYTHIA fork.
// Coordinates are printed as (t, x, y, z), in fm/c and fm respectively.
#ifndef Pythia8_HadronVertexTrace_H
#define Pythia8_HadronVertexTrace_H

#include "Pythia8/Event.h"
#include "Pythia8/FragmentationSystems.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
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
  long long event = -1;
  int attempt = 0, pass = 0, string = -1;
  std::map<int, Tag> tags;
  std::map<int, Snapshot> snapshots;
};

// One context per thread, shared by the wrapper and fragmentation routines.
inline Context& context() {
  static thread_local Context value;
  return value;
}

inline void prefix(std::ostream& out) {
  const Context& c = context();
  out << std::scientific << std::setprecision(16)
      << "[hadronization event=" << c.event << " attempt=" << c.attempt
      << " pass=" << c.pass << " string=" << c.string << "] ";
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

inline void beginEvent(long long event, int attempt) {
  Context& c = context();
  c = Context();
  c.active = true;
  c.event = event;
  c.attempt = attempt;
  std::ostringstream out;
  prefix(out);
  out << "BEGIN coordinates=(t,x,y,z), time=fm/c, position=fm, momentum=GeV\n";
  std::cerr << out.str();
}

class EventScope {
public:
  EventScope(bool enabledIn, long long event, int attempt)
    : enabled(enabledIn) {
    if (enabled) beginEvent(event, attempt);
  }
  ~EventScope() {
    if (!enabled) return;
    std::ostringstream out;
    prefix(out);
    out << "END\n";
    std::cerr << out.str();
    context() = Context();
  }
private:
  bool enabled;
};

// Also support direct PYTHIA callers without an eHIJING event label.
class PassScope {
public:
  PassScope(bool enabledIn, const Event& event)
    : enabled(enabledIn), ownsEvent(false) {
    if (!enabled) return;
    ownsEvent = !context().active;
    if (ownsEvent) beginEvent(-1, 0);
    // A failed attempt may have been rolled back by PYTHIA. Drop tags for
    // entries no longer present before their event indices can be reused.
    context().tags.erase(context().tags.lower_bound(event.size()),
      context().tags.end());
    context().snapshots.erase(context().snapshots.lower_bound(event.size()),
      context().snapshots.end());
    ++context().pass;
    context().string = -1;
  }
  ~PassScope() {
    if (ownsEvent) context() = Context();
  }
private:
  bool enabled, ownsEvent;
};

inline void beginString(int index, const ColSinglet& system,
  const Event& event, const char* branch) {
  context().string = index;
  std::ostringstream out;
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
        << " status=" << p.status() << " colors=" << p.col() << ','
        << p.acol() << " mothers=" << p.mother1() << ',' << p.mother2()
        << " hasVertex=" << p.hasVertex() << ' ';
    momentum(out, p.p());
    out << ' ';
    vertex(out, "absolute_vertex", p.vProd());
    out << '\n';
  }
  std::cerr << out.str();
}

inline void recordHadron(const Event& event, int index, const char* branch) {
  Context& c = context();
  c.tags[index] = Tag{c.pass, c.string, branch};
  const Particle& p = event[index];
  c.snapshots[index] = Snapshot{p.vProd(), p.hasProductionTime(),
    p.vBreakupFirst(), p.vBreakupSecond()};
}

inline void breakups(const char* leg, const std::vector<Vec4>& points) {
  std::ostringstream out;
  for (int i = 0; i < int(points.size()); ++i) {
    prefix(out);
    out << "BREAKUP leg=" << leg << " B" << i << ' ';
    vertex(out, "relative_vertex", points[i], 1.);
    out << '\n';
  }
  std::cerr << out.str();
}

inline void hadron(const Event& event, int index, const char* branch,
  int firstIndex, int secondIndex, const Vec4& origin,
  const Vec4& first, const Vec4& second, const Vec4& displacement,
  double redOsc, double kappa) {
  recordHadron(event, index, branch);
  const Particle& p = event[index];
  Vec4 a = origin * MM2FM + first;
  Vec4 b = origin * MM2FM + second;
  double earliest = std::min(a.e(), b.e());
  std::ostringstream out;
  prefix(out);
  out << "HADRON index=" << index << " pdg=" << p.id()
      << " status=" << p.status() << " branch=" << branch
      << " breaks=B" << firstIndex << ",B" << secondIndex
      << " redOsc=" << redOsc << " kappa_GeV_per_fm=" << kappa << '\n';
  out << "  "; vertex(out, "inherited_origin", origin);
  out << "\n  "; vertex(out, "absolute_breakup_first", a, 1.);
  out << "\n  "; vertex(out, "absolute_breakup_second", b, 1.);
  out << "\n  "; vertex(out, "prescription_displacement", displacement, 1.);
  out << "\n  "; vertex(out, "saved_vertex", p.vProd());
  out << " computed_earliest_breakup=" << earliest
      << " saved_minus_earliest=" << p.tProd() * MM2FM - earliest
      << " hasProductionTime=" << p.hasProductionTime();
  if (p.hasProductionTime()) out << " stored_prod_time=" << p.tProduction() * MM2FM;
  else out << " stored_prod_time=unavailable(no_stored_breakup_pair)";
  out << '\n';
  std::cerr << out.str();
}

inline void collapse(const Event& event, int index, int recoilString,
  const Vec4& origin) {
  recordHadron(event, index, "one-hadron-collapse");
  const Particle& p = event[index];
  std::ostringstream out;
  prefix(out);
  out << "HADRON index=" << index << " pdg=" << p.id()
      << " status=" << p.status()
      << " branch=one-hadron-collapse recoil_string=" << recoilString
      << " production_time=unavailable(no_reconstructed_breakup_pair)\n  ";
  vertex(out, "inherited_origin", origin);
  out << "\n  "; vertex(out, "prescription_displacement", p.vProd() - origin);
  out << "\n  "; vertex(out, "saved_vertex", p.vProd());
  out << '\n';
  std::cerr << out.str();
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
  context().string = -1;
  std::ostringstream out;
  prefix(out);
  out << "FINAL_EVENT success=" << success << " size=" << event.size() << '\n';
  int outputIndex = 0;
  for (int i = 0; i < event.size(); ++i) {
    const Particle& p = event[i];
    out << "  index=" << i << " pdg=" << p.id() << " status=" << p.status()
        << " mothers=" << p.mother1() << ',' << p.mother2()
        << " daughters=" << p.daughter1() << ',' << p.daughter2();
    if (success && p.isFinal()) out << " output_ID=" << outputIndex++;
    auto tag = context().tags.find(i);
    if (tag != context().tags.end())
      out << " source_pass=" << tag->second.pass
          << " source_string=" << tag->second.string
          << " branch=" << tag->second.branch;
    else if (p.isFinal() && p.isHadron()) ancestorStrings(out, event, i);
    out << " colors=" << p.col() << ',' << p.acol() << ' ';
    momentum(out, p.p());
    out << " mass_GeV=" << p.m() << " hasVertex=" << p.hasVertex() << ' ';
    vertex(out, "saved_vertex", p.vProd());
    out << " raw_saved_time_mm_per_c=" << p.tProd()
        << " hasProductionTime=" << p.hasProductionTime();
    if (p.hasProductionTime()) {
      out << " prod_time=" << p.tProduction() * MM2FM
          << " raw_prod_time_mm_per_c=" << p.tProduction() << ' ';
      vertex(out, "stored_breakup_first", p.vBreakupFirst());
      out << ' '; vertex(out, "stored_breakup_second", p.vBreakupSecond());
    } else out << " prod_time=unavailable(no_stored_breakup_pair)";
    auto snapshot = context().snapshots.find(i);
    if (snapshot != context().snapshots.end()) {
      const Snapshot& s = snapshot->second;
      out << " vertex_changed_after_fragmentation=" << !sameVertex(s.vertex, p.vProd())
          << " breakup_metadata_changed_after_fragmentation="
          << (s.hasPair != p.hasProductionTime()
              || !sameVertex(s.first, p.vBreakupFirst())
              || !sameVertex(s.second, p.vBreakupSecond()));
    }
    out << '\n';
  }
  std::cerr << out.str();
}

} // namespace HadronVertexTrace
} // namespace Pythia8
#endif
