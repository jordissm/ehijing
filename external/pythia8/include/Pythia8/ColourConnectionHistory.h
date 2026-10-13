// Diagnostic history of eHIJING medium-stage colour links. No physics is changed.
#ifndef Pythia8_ColourConnectionHistory_H
#define Pythia8_ColourConnectionHistory_H

#include "Pythia8/Event.h"
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace Pythia8 {

class ColourConnectionHistory {
public:
  struct Node {
    int id, pdg, sourceIndex, col, acol;
    bool timed;
    double time;
  };
  struct Change {
    int id, actor;
    std::string cause;
    bool timed, ordered;
    double time;
    std::vector<Node> before, after;
    std::set<int> unresolved;
  };
  struct Link {
    int id, colour, positive, negative, begin, end;
    bool timed, ordered;
    double time;
  };
  struct Connection {
    std::string status = "unavailable";
    bool timed = false, firstTimed = false, hasLatestLink = false;
    double time = 0., firstTime = 0., latestLinkTime = 0.;
    std::string firstStatus = "unavailable";
    std::set<int> sources;
    std::vector<int> links;
  };

  void reset(Event& event) {
    *this = ColourConnectionHistory();
    std::vector<Particle*> initial;
    for (int i = 0; i < event.size(); ++i) {
      event[i].colourHistoryId(0);
      if (event[i].isFinal() && event[i].isParton()) {
        identify(event[i], i, false, 0.);
        initial.push_back(&event[i]);
      }
    }
    apply("pre_medium_state", false, 0., 0, initial);
  }

  int identify(Particle& p, int index, bool timed, double time) {
    const bool newlyCreated = timed;
    timed = timed && std::isfinite(time) && time >= 0.;
    const int id = ++nextId;
    p.colourHistoryId(id);
    nodes[id] = Node{id, p.id(), index, newlyCreated ? 0 : p.col(), newlyCreated ? 0 : p.acol(), timed, time};
    return id;
  }

  void bind(int id, int sourceIndex) {
    if (nodes.count(id)) nodes[id].sourceIndex = sourceIndex;
  }

  // Times are sampled medium-clock times relative to the hard interaction,
  // converted onsite to fm/c. They are not inferred four-vertices.
  void apply(const std::string& cause, bool timed, double time, int actor,
    const std::vector<Particle*>& updated) {
    timed = timed && std::isfinite(time) && time >= 0.;
    Change change{int(changes.size()) + 1, actor, cause, timed, true, time, {}, {}, {}};
    for (const Particle* p : updated) {
      const int id = p->colourHistoryId();
      if (!nodes.count(id)) continue;
      change.before.push_back(nodes[id]);
      if (timed && ((lastTime.count(id) && time < lastTime[id])
          || (nodes[id].timed && time < nodes[id].time))) change.ordered = false;
      nodes[id].col = p->col();
      nodes[id].acol = p->acol();
      change.after.push_back(nodes[id]);
      if (timed) lastTime[id] = lastTime.count(id) ? std::max(time, lastTime[id]) : time;
    }
    std::map<int, std::vector<int>> positive, negative;
    for (const auto& entry : nodes) {
      const Node& n = entry.second;
      if (n.col > 0) positive[n.col].push_back(n.id);
      if (n.acol > 0) negative[n.acol].push_back(n.id);
    }
    std::set<Key> present;
    for (const auto& entry : positive) {
      const int tag = entry.first;
      // Junctions and ambiguous/repeated tags require a different graph model.
      if (entry.second.size() == 1 && negative[tag].size() == 1)
        present.insert(Key(tag, entry.second[0], negative[tag][0]));
    }
    for (const auto& entry : nodes) {
      const Node& n = entry.second;
      if (n.col < 0 || n.acol < 0
          || (n.col && (positive[n.col].size() != 1 || negative[n.col].size() != 1))
          || (n.acol && (positive[n.acol].size() != 1 || negative[n.acol].size() != 1)))
        change.unresolved.insert(n.id);
    }
    changes.push_back(change);
    for (auto it = active.begin(); it != active.end();) {
      if (!present.count(it->first)) {
        links[it->second].end = change.id;
        it = active.erase(it);
      } else ++it;
    }
    for (const Key& key : present) if (!active.count(key)) {
      bool ordered = change.ordered;
      for (int id : {std::get<1>(key), std::get<2>(key)})
        if (timed && ((nodes[id].timed && time < nodes[id].time)
            || (lastTime.count(id) && time < lastTime[id]))) ordered = false;
      links.push_back(Link{int(links.size()) + 1, std::get<0>(key),
        std::get<1>(key), std::get<2>(key), change.id, 0, timed, ordered, time});
      active[key] = int(links.size()) - 1;
    }
  }

  bool sourceIds(const Event& event, int index, std::set<int>& sources) const {
    std::set<int> visited;
    bool complete = true;
    lineage(event, index, visited, sources, complete);
    return complete;
  }

  Connection connection(const Event& event, const std::vector<int>& partons) const {
    Connection result;
    std::set<int> visited;
    bool complete = true;
    for (int index : partons) {
      if (index < 0) { complete = false; continue; }
      lineage(event, index, visited, result.sources, complete);
    }
    if (!complete || result.sources.size() < 2) {
      result.status = "unavailable_missing_or_changed_provenance";
      result.firstStatus = result.status;
      return result;
    }
    std::map<int, std::set<int>> neighbours;
    bool ordered = true, external = false;
    for (const auto& entry : active) {
      const Link& link = links[entry.second];
      const bool pos = result.sources.count(link.positive), neg = result.sources.count(link.negative);
      if (pos != neg) external = true;
      if (!pos || !neg) continue;
      result.links.push_back(link.id);
      neighbours[link.positive].insert(link.negative);
      neighbours[link.negative].insert(link.positive);
      ordered = ordered && link.ordered;
      if (link.timed) {
        result.latestLinkTime = result.hasLatestLink ? std::max(result.latestLinkTime, link.time) : link.time;
        result.hasLatestLink = true;
      }
    }
    std::set<int> reached;
    std::vector<int> pending(1, *result.sources.begin());
    while (!pending.empty()) {
      const int id = pending.back(); pending.pop_back();
      if (!reached.insert(id).second) continue;
      for (int neighbour : neighbours[id]) pending.push_back(neighbour);
    }
    bool unresolved = false;
    if (!changes.empty()) for (int id : result.sources)
      if (changes.back().unresolved.count(id)) unresolved = true;
    if (external || reached != result.sources || unresolved) result.status = "unavailable_incomplete_component";
    else if (!ordered) result.status = "unavailable_nonchronological_model_operations";
    else {
      // The component can first arise through REMOVAL of an external link.
      // The newest surviving internal link alone does not date its creation.
      bool wasConnected = false, seen = false;
      for (const Change& change : changes) {
        bool connected = componentAt(result.sources, change.id);
        if (connected && !wasConnected) {
          bool causal = change.ordered;
          for (int id : result.sources)
            if (change.timed && nodes.at(id).timed && nodes.at(id).time > change.time) causal = false;
          result.timed = change.timed && causal;
          result.time = change.time;
          result.status = !causal ? "unavailable_nonchronological_model_operations"
            : change.cause == "pre_medium_state" ? "unavailable_pre_medium_creation"
            : change.timed ? "recorded_medium_component_creation" : "unavailable_model_clock";
          if (!seen) {
            result.firstTimed = result.timed;
            result.firstTime = result.time;
            result.firstStatus = result.status;
            seen = true;
          }
        }
        wasConnected = connected;
      }
    }
    if (result.firstStatus == "unavailable") result.firstStatus = result.status;
    return result;
  }

  void write(std::ostream& out) const {
    out << std::defaultfloat << std::setprecision(12) << std::boolalpha
        << "\n=== Colour connections before hadronization ===\n"
        << "Times are sampled model-clock times [fm/c] relative to the hard interaction.\n"
        << "Sequence is algorithmic order; it need not equal physical time order.\n"
        << "Pre-medium connections have no recorded creation time; no zero is assigned.\n"
        << "A link is identified by its colour tag AND both stable source IDs.\n";
    for (const auto& entry : nodes) {
      const Node& n = entry.second;
      out << "COLOUR_SOURCE id=" << n.id << " source_index=" << n.sourceIndex
          << " pdg=" << n.pdg << " birth_time=";
      if (n.timed) out << n.time; else out << "unavailable";
      out << '\n';
    }
    for (const Change& change : changes) {
      out << "COLOUR_CHANGE id=" << change.id << " cause=" << change.cause
          << " actor=" << change.actor << " chronological=" << change.ordered << " time=";
      if (change.timed) out << change.time; else out << "unavailable";
      out << " unresolved_sources=(";
      bool first = true;
      for (int id : change.unresolved) { if (!first) out << ','; out << id; first = false; }
      out << ")\n";
      for (size_t i = 0; i < change.after.size(); ++i) {
        const Node& before = change.before[i]; const Node& after = change.after[i];
        out << "  source=" << after.id << " old_colors=" << before.col << ',' << before.acol
            << " new_colors=" << after.col << ',' << after.acol << '\n';
      }
    }
    for (const Link& link : links) {
      out << "COLOUR_LINK id=" << link.id << " colour=" << link.colour
          << " positive_source=" << link.positive << " negative_source=" << link.negative
          << " begin_change=" << link.begin << " end_change=" << link.end
          << " chronological=" << link.ordered << " time=";
      if (link.timed) out << link.time; else out << "unavailable";
      out << '\n';
    }
  }

private:
  typedef std::tuple<int, int, int> Key;
  int nextId = 0;
  std::map<int, Node> nodes;
  std::map<int, double> lastTime;
  std::vector<Change> changes;
  std::vector<Link> links;
  std::map<Key, int> active;

  bool componentAt(const std::set<int>& sources, int sequence) const {
    for (int id : sources) if (changes[sequence - 1].unresolved.count(id)) return false;
    std::map<int, std::set<int>> neighbours;
    for (const Link& link : links) {
      if (link.begin > sequence || (link.end && link.end <= sequence)) continue;
      const bool pos = sources.count(link.positive), neg = sources.count(link.negative);
      if (pos != neg) return false;
      if (pos && neg) {
        if (!link.ordered) return false;
        neighbours[link.positive].insert(link.negative);
        neighbours[link.negative].insert(link.positive);
      }
    }
    std::set<int> reached;
    std::vector<int> pending(1, *sources.begin());
    while (!pending.empty()) {
      int id = pending.back(); pending.pop_back();
      if (!reached.insert(id).second) continue;
      for (int neighbour : neighbours[id]) pending.push_back(neighbour);
    }
    return reached == sources;
  }

  void lineage(const Event& event, int index, std::set<int>& visited,
    std::set<int>& sources, bool& complete) const {
    if (index <= 0 || index >= event.size()) { complete = false; return; }
    if (!visited.insert(index).second) return;
    const Particle& p = event[index];
    const int id = p.colourHistoryId();
    if (id > 0 && nodes.count(id)) {
      sources.insert(id);
      if (nodes.at(id).col != p.col() || nodes.at(id).acol != p.acol()) complete = false;
      return;
    }
    // Native joining entries have two mothers; momentum-only copies preserve
    // the stable ID. Never infer identity by matching flavour or momentum.
    if (p.statusAbs() < 71 || p.statusAbs() > 74 || p.mother1() <= 0
        || p.mother1() >= index || p.mother2() < 0 || p.mother2() >= index) {
      complete = false; return;
    }
    lineage(event, p.mother1(), visited, sources, complete);
    if (p.mother2() > 0) lineage(event, p.mother2(), visited, sources, complete);
  }
};

} // namespace Pythia8
#endif
