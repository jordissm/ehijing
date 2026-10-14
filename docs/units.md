# Unit conventions {#unit_conventions}

This document defines the unit conventions and onsite conversions at the
eHIJING, PYTHIA, and OSCAR/SMASH interfaces. The supported shower workflow has
`PartonVertex:setVertex = off` and imports parton vertices as displacements
relative to the hard interaction. The hadronizer adds the nuclear origin once.

## Conventions by quantity

| Quantity | Unit |
| --- | --- |
| eHIJING nuclear coordinates, path lengths, collision times and formation times | GeV⁻¹ |
| Adapter spacetime quantities after multiplication by `hbarc_gev_fm` | fm for length, fm/c for time |
| Standard PYTHIA particle vertices and stored breakup vertices | mm for length, mm/c for time |
| PYTHIA `PartonVertex` subsystem in this version | fm for length, fm/c for time |
| Local PYTHIA string breakup displacements before conversion to stored vertices | fm for length, fm/c for time |
| OSCAR/SMASH positions and coordinate times | fm for length, fm/c for time |
| Particle momentum, energy and mass | GeV |

eHIJING's core medium calculations use natural units, with hbar = c = 1.
Nuclear density is expressed in GeV³ and nuclear thickness in GeV². Geometry
must retain its natural units while used in those formulas. Converting a
coordinate for a particle vertex does not change the core geometry convention.

PYTHIA's `PartonVertex` exception is documented in the bundled
`external/pythia8/share/Pythia8/xmldoc/VertexInformation.xml`. Standard particle
vertices use mm, but this subsystem assigns values in fm for rope calculations.
The supplied hard-process settings explicitly disable this subsystem. The
separate hadronizer has `ProcessLevel:all = off`, so its ordinary MPI/ISR/FSR
vertex callbacks are not exercised in this workflow. Enabling the subsystem
in the source shower would require an explicit conversion policy during import;
its units cannot be inferred from the `Vec4` type alone.

## Conversion factors

The driver constants live in `ehijing::constants::units`, in
`include/ehijing_constants.hpp`:

| Conversion | Operation |
| --- | --- |
| GeV⁻¹ to fm | Multiply by `hbarc_gev_fm` = 0.19732698 GeV fm |
| fm to mm | Multiply by `fm_to_mm` = 10⁻¹² |
| mm to fm | Multiply by `mm_to_fm` = 10¹² |

With c = 1, the same factors apply to time. In particular, the fm/mm factors
convert fm/c to mm/c and back. A spacetime `Pythia8::Vec4` is constructed as
`(x, y, z, ct)`; its `.e()` accessor returns the fourth coordinate. Scale all
four components together when converting a vertex or displacement. Momentum
vectors also use `Vec4`, but remain in GeV and must not receive these length
conversions.

Conversions are performed onsite at the interface that requires them. A unit
conversion preserves the reference frame and coordinate origin; it is neither
a Lorentz transformation nor a translation.

## Boundary rules

1. Retain eHIJING nuclear coordinates and sampled times in GeV⁻¹ for medium
   calculations. Convert to fm when constructing adapter spacetime quantities.
2. `Modified_FF` converts its newly appended parton displacements from fm to mm.
   Those displacements exclude the nuclear origin. Existing source parton
   vertices must likewise be relative displacements in mm.
3. `Hadronizer` converts the nuclear origin from GeV⁻¹ to mm, copies imported
   displacements in mm, and adds the origin once to each parton. Its reconstructed
   shower displacement, calculated in fm, is also converted onsite to mm.
   Vertices entering standard fragmentation are therefore absolute coordinates
   in mm. PYTHIA converts its local string breakup displacements from fm to mm
   before adding them to that inherited origin.
4. The writer converts the complete formation vertex and the separate production
    time from mm to fm before OSCAR freestreaming. The writer adds no spatial
    smearing to hadron positions. Production
    references can include the inherited string/collapse origin as well as breakup
    vertices; all are stored in mm and transformed as full four-vertices. Spectator
   coordinates generated in fm retain those units in legacy mode. In explicit
   target mode, rest-frame nucleon positions in GeV^-1 are converted to fm,
   boosted as full four-positions, and freestreamed to the output t=0 surface.
   Rest-frame spectator momenta are boosted with the same target four-momentum.
5. Keep velocities dimensionless: `p/E`. With positions in fm and times in fm/c,
   freestream to t = 0 using `x(0) = x(form) - (p/E) * form_time`, with c = 1.

`prod_time` and `form_time` are coordinate times in the particle reference
frame. They are not proper-time intervals. SMASH's configured
`proper_formation_time_` is a separate proper-time quantity. Frame consistency
and the downstream use of imported formation times require their own review.

## Explicit nuclear-target frame handling

The frozen Gaussian density, forward-ray thickness, collision times and exact
elastic recoil are evaluated in the nuclear rest frame. Geometry stays in
GeV^-1, density in GeV^3 and thickness in GeV^2. The low-scale adapter boosts
source momenta and relative four-vertices back to that frame for processing,
then boosts accepted results into the original generator frame. Formation
variables remain ordering variables rather than measured formation events.

Appended parton displacements are converted from GeV^-1 to mm and boosted as
full four-vectors. The hadronizer separately converts and boosts the nuclear
origin `(Rx,Ry,Rz,0)` before adding it once. Its fourth coordinate can become
nonzero after a boost. Diagnostic nuclear snapshots and elastic-collision
observations are already absolute mm four-positions in the generator frame.
Spectator emission uses the same frame conversion and output t=0 surface as
the hadrons. This does not restore missing high-scale shower delays or provide
a dynamical nuclear trajectory model. Density centres are frozen even when
whole-remnant momenta change.

## Diagnostic history interface

`event_XXXXXXXX.history` uses native HepMC3 Asciiv3 ASCII with `U GEV MM`.
Observation vertices retain the current stage's origin/frame convention;
collection performs no Lorentz transformation and adds no nuclear origin.
`stored_vertex_xyzt` preserves the raw `(x,y,z,ct)` value with an explicit unit
label; nuclear geometry coordinates and collision times remain separately
labelled GeV^-1 quantities. String production references and breakup pairs use
absolute mm coordinates. The `oscar_output` snapshot converts the already
computed t = 0 output positions from fm to mm on diagnostic copies only.

A missing particle vertex is flagged by `has_vertex=0`, even though HepMC's
position field contains a zero. Observation vertices (status 900) and auxiliary
serialization anchors (`diagnostic_anchor=1`) are diagnostic bookkeeping, not
new physical interaction assignments. Snapshot sequence numbers are execution
order, not a coordinate-time ordering. See the README for record identity,
source mappings and comment attributes.

## Remaining review items

- Custom source settings or external parton producers must satisfy the import
  convention: relative displacements in mm. `hasVertex()` does not identify
  units or indicate whether a nuclear origin has already been included.
- The hadronizer imports final partons with status 23 and without parent links.
  The subsequent ancestry traversal therefore reconstructs no shower delay:
  `position` and `timeplus` remain zero. Medium displacements assigned by
  `Modified_FF` are still preserved, and PYTHIA still computes hadron breakup
  and formation vertices. Restoring the source shower history is a separate
  change to the spacetime model.
- Custom/legacy coordinate-frame consistency and SMASH's use of the imported
  `form_time` require separate review. Explicit-target boosts are implemented
  at the eHIJING interfaces; downstream transport timing remains separate. Unit conversion does not address either issue.
