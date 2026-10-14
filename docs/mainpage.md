# eHIJING Documentation

eHIJING is the event-generation component that combines PYTHIA DIS event
generation with eHIJING medium modification, hadronization, and OSCAR2013-style
output.

See [Physical assumptions and open questions](assumptions.md) for the current
model scope, known implementation limits and unresolved consistency checks.

The executable workflow is:

1. Parse runtime configuration into `RunConfig`.
2. Configure and initialize the PYTHIA hard process.
3. Reconstruct `DISKinematics` and apply runtime `DISCuts`.
4. Apply low-\f$ Q^2 \f$ medium modification with `Modified_FF`.
5. Convert the modified partonic event into hadrons with `Hadronizer`.
6. Write OSCAR2013 event records and JSON metadata.

## Event generation flow

The diagram follows the executable's control flow, with the default explicit
nuclear target. The blue block runs inside `pythia.next()`; DIS cuts are applied
after that complete generation step. The teal block is the low-scale adapter,
the purple block is the dedicated hadronizer, and green marks accepted output.
Dashed red paths discard the complete attempt and start a fresh hard event
and nucleus. The detailed physics equations and source excerpts follow below.

\dot
digraph ehijing_event_flow {
  graph [rankdir=TB, bgcolor="transparent", pad=0.18,
         nodesep=0.28, ranksep=0.32, compound=true,
         fontname="Helvetica", fontsize=14];
  node [shape=box, style="rounded,filled", fillcolor="#F8FAFC",
        color="#94A3B8", fontcolor="#172B4D", fontname="Helvetica",
        fontsize=13, penwidth=1.1, margin="0.15,0.10"];
  edge [color="#64748B", fontcolor="#475569", fontname="Helvetica",
        fontsize=11, arrowsize=0.7, penwidth=1.1];

  setup [label="Configure and initialize\nRunConfig + hard / hadronization / DIS-cut files\nA, Z, K, n, lambda, HT/GHT mode, geometry\nseeds, Qs/GHT tables, output paths"];
  attempt [label="A. Start a fresh generation attempt\naccepted-event counter is unchanged",
           fillcolor="#E2E8F0"];

  subgraph cluster_generation {
    label="PYTHIA generation: pythia.next()";
    labelloc=t; labeljust=l;
    style="rounded,filled"; color="#93C5FD"; fillcolor="#EFF6FF";
    hard [label="Sample the hard DIS process\nper-nucleon beam + nuclear PDF / isospin\nneutral-current exchange; vacuum ISR",
          fillcolor="#DBEAFE"];
    geometry [label="Assign the hard vertex and medium geometry\non: A identities; Woods-Saxon / Hulthen centers\nGaussian nucleons; primary excluded from rescattering\noff: legacy uniform hard-sphere geometry",
              fillcolor="#DBEAFE"];
    collisions [label="For each shower radiator, sample the medium\nforward-ray thickness TA -> self-consistent Qs\ngluon TMD -> transfers, times, azimuths\nexplicit mode also stores target identities",
                fillcolor="#DBEAFE"];
    shower [label="Develop the high-scale parton shower\npT ordering; choose either HT or GHT + LPM weights\nmedium term above Qsep; vacuum shower continues below\nsave final partons, collision histories, separation scale",
            fillcolor="#DBEAFE"];
    hard -> geometry -> collisions -> shower;
  }

  cuts [shape=diamond, style=filled, fillcolor="#FEF3C7",
        color="#D97706", label="Configured DIS cuts pass?\nxB, y, nu, Q2, W2"];

  subgraph cluster_medium {
    label="Low-scale adapter: Modified_FF::sample_ff_partons()";
    labelloc=t; labeljust=l;
    style="rounded,filled"; color="#5EEAD4"; fillcolor="#F0FDFA";
    radiation [label="Sample additional medium-induced gluons\nfinal light partons with sampled collisions\nincreasing formation variable; HT or GHT vetoes\nmatch below stored Qsep; new gluons sample later hits",
               fillcolor="#CCFBF1"];
    recoil [label="Apply sampled elastic hits in time order\nnuclear rest frame; shared remnant per identity\nsolve exact two-body recoil for the prescribed transfer\nretain the existing low-scale radiation approximation",
            fillcolor="#CCFBF1"];
    remnants [label="Build the remnant color system\nexchange existing color lines on repeated hits\nsplit each touched nonprimary nucleon once into q + diquark\nkeep untouched identities for spectator output",
              fillcolor="#CCFBF1"];
    radiation -> recoil -> remnants;
  }

  hadronize [label="Hadronizer: separate PYTHIA instance\nimport final partons, colors and relative vertices\nadd the hard origin once; Lund string fragmentation\nassign production / formation vertices; configured decays",
             fillcolor="#EDE9FE", color="#A78BFA"];
  balance [shape=diamond, style=filled, fillcolor="#FEF3C7",
           color="#D97706",
           label="Explicit target ledger agrees?\nfinal baryon number / charge\nvs. untouched target identities"];
  output [label="Write the accepted event\nOSCAR: hadrons + spectators on the t = 0 surface\nJSON: DIS metadata + nuclear bookkeeping\noptional HepMC history / fragmentation trace\nadvance the accepted-event counter",
          fillcolor="#DCFCE7", color="#4ADE80"];
  enough [shape=diamond, style=filled, fillcolor="#FEF3C7",
          color="#D97706", label="Requested accepted\nevent count reached?"];
  finish [label="Finish eHIJING\nreport generation and rejection counts",
          fillcolor="#DCFCE7", color="#4ADE80"];
  smash [label="ELECTRA downstream: SMASH\nhadronic transport using the exported event\nseparate stage with its own formation / interaction model",
         style="rounded,dashed,filled", fillcolor="#F1F5F9"];

  reject [label="Reject the attempt\nno output; keep event ID",
          fillcolor="#FEE2E2", color="#F87171", fontcolor="#991B1B"];
  retry_a [shape=circle, label="A", width=0.35, fixedsize=true,
           fillcolor="#FEE2E2", color="#F87171", fontcolor="#991B1B"];
  next_a [shape=circle, label="A", width=0.35, fixedsize=true,
          fillcolor="#DCFCE7", color="#4ADE80"];

  setup -> attempt -> hard;
  shower -> cuts [label="generated"];
  cuts -> radiation [label="yes"];
  remnants -> hadronize;
  hadronize -> balance [label="fragmentation succeeded"];
  balance -> output [label="yes (legacy: skip identity check)"];
  output -> enough;
  enough -> finish [label="yes"];
  enough -> next_a [label="no: next event", color="#16A34A",
                    fontcolor="#15803D"];
  finish -> smash [style=dashed, label="downstream handoff"];

  edge [style=dashed, color="#DC2626", fontcolor="#B91C1C", constraint=false];
  shower -> reject;
  cuts -> reject [label="no"];
  recoil -> reject;
  remnants -> reject;
  hadronize -> reject;
  balance -> reject [label="no"];
  reject -> retry_a [label="retry"];
  {rank=same; balance; reject;}
  {rank=same; output; retry_a;}
  {rank=same; finish; next_a;}
}
\enddot

**How to read the flow.** Within the shower, collisions and radiation share
the sampled medium information; the diagram groups that work into two boxes
for readability. Each circular **A** connector returns to the start of a fresh
generation attempt at the top. HT and GHT are alternative kernels, and the high-scale
shower and low-scale adapter cover different emission regions. The geometry
is frozen even when momenta recoil. Radiation is sampled before the collected
elastic updates are sorted, so the diagram's arrows describe execution order,
rather than a complete ordering of physical formation events.

**Acceptance and optional branches.** Only a successful attempt writes an
event or its optional diagnostics. The history collector observes existing
states without changing the random-number sequence. The driver's retries
are grouped into 15-attempt batches; exhausting a batch starts another one.
Explicit-target runs stop with a diagnostic after 10000 consecutive attempts
without an accepted event. Initialization, a missing explicit target, or output
I/O failures terminate the run instead of following the red retry paths.
The B/Q identity check and shared exact remnant recoil apply to explicit-target
mode; the legacy branch uses the prescriptions described in the physics
section. SMASH is a separate ELECTRA stage after eHIJING output.

## Physics included in eHIJING

This section describes the physics implemented in this checkout, including
the ELECTRA extensions to the original eHIJING model. The theoretical reference
is W. Ke, Y.-Y. Zhang, H. Xing and X.-N. Wang,
[eHIJING: an Event Generator for Jet Tomography in Electron-Ion Collisions](https://arxiv.org/html/2304.10779v2),
[Phys. Rev. D 110, 034001 (2024)](https://doi.org/10.1103/PhysRevD.110.034001).
Equation numbers cited below refer to that paper; equations displayed here use
one consistent notation and identify numerical choices made in the source.
In particular, the explicit nuclear target, shared remnant recoil and spacetime
interfaces described here must be read together with the implementation, rather
than inferred from the published calculation alone.

The simulation connects a short-distance lepton-parton scattering to a QCD
shower, parton interactions with cold nuclear matter, a color-connected remnant
system and string fragmentation. Its medium is a sampled, frozen nuclear
configuration. It has no temperature field, hydrodynamic evolution or dynamically
evolving nuclear density. Hadronic transport belongs to the subsequent SMASH
stage in ELECTRA.

The implementation excerpts below are taken from this checkout. Source
paths are relative to `src/ehijing/`; original variable names are retained.
Surrounding setup and control flow are omitted where indicated by
`// ...`, and indentation is adjusted for the manual. The excerpts are
reading aids, rather than standalone programs; update them alongside their
source when a physics prescription changes.

### Kinematics, notation and units

Let \f$k^\mu\f$ and \f$k'^\mu\f$ be the incoming and scattered lepton
momenta, and \f$P^\mu\f$ the incoming target momentum **per nucleon** used by
the hard generator. This documentation follows `compute_dis_kinematics()` and
defines the momentum transferred to the hadronic system by

\f[
 q^\mu=k^\mu-k'^\mu,\qquad Q^2=-q^2>0,\qquad
 x_B=\frac{Q^2}{2P\cdot q},\qquad
 y=\frac{P\cdot q}{P\cdot k}.
\f]

The target-rest-frame photon energy and the hadronic invariant mass are

\f[
 \begin{aligned}
 \nu&=\frac{P\cdot q}{\sqrt{P^2}},\\
 W^2&=(P+q)^2=m_N^2+2m_N\nu-Q^2\\
    &=m_N^2+Q^2\left(\frac{1}{x_B}-1\right).
 \end{aligned}
\f]

The expressions involving \f$m_N\f$ assume \f$P^2=m_N^2\f$. The code
reconstructs these invariants from the event's leptons and target, then applies
the configured `DISCuts`; cuts select events and do not change the interaction
kernel. For a detected hadron, a useful invariant energy fraction is
\f$z_h=P\cdot p_h/(P\cdot q)=E_h/\nu\f$ in the nuclear rest frame.
Hadron transverse momentum in SIDIS is measured relative to the exchanged
boson direction. It should be distinguished from the transverse variables used
internally for parton collisions and branchings.

Throughout the medium discussion, \f$\boldsymbol{\kappa}_i\f$ denotes the
transverse momentum transferred in collision \f$i\f$ and
\f$\boldsymbol{\ell}\f$ denotes the transverse momentum in a radiative
splitting. Thus the paper's collision momentum \f$\mathbf{k}_\perp\f$ is
written as \f$\boldsymbol{\kappa}\f$ here to avoid confusion with the lepton
momentum. The emitted gluon carries fraction \f$z\f$ of its parent's energy
in the high-energy approximation; \f$z\f$ is distinct from \f$z_h\f$.

The metric is \f$p^2=E^2-\mathbf{p}^2\f$. For comparison with the paper,
\f$p^\pm=(E\pm p_z)/\sqrt{2}\f$ gives
\f$p\cdot r=p^+r^-+p^-r^+-\mathbf{p}_\perp\cdot\mathbf{r}_\perp\f$.
The numerical medium calculations use the nuclear rest frame and natural
units, \f$\hbar=c=1\f$:

\f[
 \begin{gathered}
 [E]=[p]=\mathrm{GeV},\qquad [r]=[t]=\mathrm{GeV}^{-1},\\
 [\rho_N]=\mathrm{GeV}^{3},\qquad [T_A]=\mathrm{GeV}^{2},\qquad
 [\hat q]=\mathrm{GeV}^{3}.
 \end{gathered}
\f]

Multiplication by \f$\hbar c=0.19732698\,\mathrm{GeV\,fm}\f$ converts
natural-unit lengths to fm. A transport coefficient in \f$\mathrm{GeV}^3\f$
is divided by this factor to obtain \f$\mathrm{GeV}^2/\mathrm{fm}\f$.
See @ref unit_conventions for the additional mm conversions at PYTHIA interfaces.

**Implementation: DIS reconstruction.** `src/dis_kinematics.cpp`, `compute_dis_kinematics()`.

```cpp
kin.pGamma  = kin.peIn - kin.peOut; // four-momentum of the virtual photon/Z^0/W^±

const double P_dot_q = kin.pProton * kin.pGamma;
const double P_dot_k = kin.pProton * kin.peIn;

kin.nu = P_dot_q / std::sqrt(kin.pProton * kin.pProton);
kin.Q2 = -kin.pGamma.m2Calc();
kin.Q  = std::sqrt(kin.Q2);
kin.W2 = (kin.pProton + kin.pGamma).m2Calc();
kin.W  = std::sqrt(kin.W2);
kin.bjorken_x = kin.Q2 / (2.0 * P_dot_q);
kin.y  = P_dot_q / P_dot_k;
```

Here `kin.pGamma` is the transferred four-momentum, `kin.pProton` is the
per-nucleon target momentum, and PYTHIA's `Vec4` multiplication evaluates
the Minkowski scalar product used in the definitions above.

### Hard scattering, nuclear PDFs and the vacuum baseline

PYTHIA generates the hard process and the vacuum initial- and final-state
showers. The supplied HERMES configuration enables neutral-current
\f$\gamma^*/Z\f$ exchange and disables charged-current exchange. In the
photon-dominated, unpolarized, leading-order limit, the familiar parton-model
normalization per target nucleon is

\f[
 \begin{gathered}
 \frac{d^2\sigma^{eA}}{dx_B\,dQ^2}
 =\frac{2\pi\alpha_{\rm EM}^2}{x_B Q^4}
   \left[1+(1-y)^2\right]F_2^A(x_B,Q^2),\\
 F_2^A=x_B\sum_f e_f^2\left[f_{f/A}+f_{\bar f/A}\right].
 \end{gathered}
\f]

This expression neglects lepton and target masses, the longitudinal structure
function and weak-exchange terms. It establishes notation, rather than replacing
PYTHIA's process implementation. The incoming hadron is a per-nucleon beam;
the explicit target's \f$A\f$ nucleons are not generated as \f$A\f$
independent hard DIS events.

For interpreting the nuclear PDF input, write

\f[
 f_{i/A}(x,Q^2)=\frac{Z}{A}f_{i/p|A}(x,Q^2)
              +\frac{A-Z}{A}f_{i/n|A}(x,Q^2),\qquad
 f_{i/p|A}=R_i^A f_{i/p}.
\f]

The isospin mixture and the fitted nuclear modification \f$R_i^A\f$ affect
the flavor and momentum fraction selected at the hard vertex. They are separate
from the gluon TMD that controls final-state rescattering. An nPDF does not
itself supply a sampled nucleon position, a Fermi momentum or a recoil history.
The local driver sets `PDF:nPDFSetA = 3` (EPPS16 NLO) when \f$Z>2\f$ and
`0` (isospin only) otherwise, overriding that setting from the hard-process
file. The bundled PYTHIA also supports EPS09 LO/NLO choices, but changing the
driver's prescription is a separate physics choice. The HERMES file enables
`PDF:useHardNPDFA`; using an NLO PDF set does not make the generated hard
matrix element an NLO calculation.

The vacuum baseline includes perturbative radiation and its recoil,
primordial parton transverse momentum from PYTHIA's configured model, and the
transverse momentum generated during string fragmentation. Consequently a
hadron's transverse momentum is not just the sum of medium kicks. In the
current HERMES bundle, `SpaceShower:dipoleRecoil = on` and
`TimeShower:interleave = off`. The bundled `TimeShower` explicitly retains FSR
when eHIJING is enabled. These are the settings of this checkout; the paper's
choice of global recoil must not be assumed to describe them.

**Implementation: nuclear PDF selection.** `src/main.cpp`, `main()`.

```cpp
int nuclear_modification = (atomic_number > 2) ? 3 : 0;

// ...

add_arg<int>(pythia, "PDF:nPDFSetA", nuclear_modification);
add_arg<int>(pythia, "PDF:nPDFBeamA", nucleus_PDG_id);
```

The two excerpts show the species-dependent PDF choice and its application
after reading the hard-process settings. The numerical target identifier
`nucleus_PDG_id` is constructed from the requested atomic and mass numbers.

### Nuclear geometry and the path through matter

The ELECTRA executable defaults to `--nuclear-target on`. It samples exactly
\f$A\f$ identified nucleons, including exactly \f$Z\f$ protons. A primary
species is chosen with probability \f$Z/A\f$, then an identity uniformly
within that species; its center is the hard-interaction position
\f$\mathbf r_0\f$. This geometrical assignment does not resample the hard
PDF or give the primary nucleon a Fermi momentum. The primary is excluded
from subsequent medium collisions by identity.

For \f$A>2\f$, the sampled centers follow a spherical Woods--Saxon law,
including the radial volume measure:

\f[
 \begin{gathered}
 \rho_{\rm WS}(r)\propto\frac{1}{1+\exp[(r-R)/a]},\\
 p(r)\,dr\propto\frac{r^2\,dr}{1+\exp[(r-R)/a]},\qquad
 0\leq r\leq R+10a.
 \end{gathered}
\f]

The radius \f$R\f$ and surface diffuseness \f$a\f$ have tabulated defaults
for the supported heavy targets, listed in the README, or must both be supplied
explicitly. The centers are independent, without a minimum separation,
recentering, deformation or short-range correlation model. For deuterium,
\f$A=2,Z=1\f$, the separation distribution is Hulthén:

\f[
 p(r)\,dr\propto
 \left(e^{-\alpha_H r}-e^{-\beta_H r}\right)^2dr,\qquad
 \alpha_H=0.228\,\mathrm{fm}^{-1},\quad
 \beta_H=1.18\,\mathrm{fm}^{-1}.
\f]

The proton and neutron are placed at opposite half-separations. For
\f$A=1\f$, the only center is at the origin; excluding the primary leaves
no rescattering medium.

Each nonprimary identity contributes a normalized Gaussian of width
\f$w\f$:

\f[
 \rho_j(\mathbf r)=\frac{1}{(2\pi w^2)^{3/2}}
  \exp\!\left[-\frac{|\mathbf r-\mathbf r_j|^2}{2w^2}\right],\qquad
 \rho_N(\mathbf r)=\sum_{j\ne j_{\rm primary}}\rho_j(\mathbf r).
\f]

The default `--nucleon-width-fm 0.5` is a geometry benchmark, without an
associated physics fit in this checkout. A parton starting at \f$\mathbf r_s\f$
samples the forward ray \f$\mathbf r(s)=\mathbf r_s+s\mathbf n\f$,
where \f$\mathbf n=\mathbf v/|\mathbf v|\f$ and
\f$\mathbf v=\mathbf p/E\f$. Define
\f$d_j=(\mathbf r_j-\mathbf r_s)\cdot\mathbf n\f$ and
\f$b_j^2=|\mathbf r_j-\mathbf r_s|^2-d_j^2\f$. The code integrates each
Gaussian analytically:

\f[
 \begin{aligned}
 T_j&=\int_0^\infty ds\,\rho_j(\mathbf r_s+s\mathbf n)\\
 &=\frac{\exp[-b_j^2/(2w^2)]}{4\pi w^2}
   \operatorname{erfc}\!\left(-\frac{d_j}{\sqrt{2}w}\right),\\
 T_A&=\sum_{j\ne j_{\rm primary}}T_j.
 \end{aligned}
\f]

Collision identity and distance are jointly sampled with
\f$\Pr(j)=T_j/T_A\f$ and
\f$p(s|j)=\rho_j(\mathbf r_s+s\mathbf n)/T_j\f$, for \f$s\geq0\f$.
The corresponding nuclear-rest-frame time is \f$t=s/|\mathbf v|\f$.
Gaussian tails eliminate a sharp nuclear boundary. The frozen centers stay
fixed even when remnants acquire momentum, and the originally sampled rays
are not recalculated after kicks. This is an eikonal geometry prescription,
rather than a propagation of a dynamically recoiling nucleus.

**Implementation: Gaussian ray thickness.** `external/pythia8/include/eHIJING/NuclearTarget.h`, `NuclearTarget::path_weights()`.

```cpp
std::vector<double> path_weights(const std::array<double, 3>& start,
                                 const std::array<double, 3>& velocity) const {
    const double speed = norm(velocity);
    std::vector<double> weights(nucleons_.size(), 0.);
    if (!(speed > 0.)) return weights;
    for (size_t j = 0; j < nucleons_.size(); ++j) {
        if (static_cast<int>(j) == primary_) continue;
        double longitudinal = 0., distance2 = 0.;
        for (int k = 0; k < 3; ++k) {
            const double d = nucleons_[j].position[k] - start[k];
            longitudinal += d * velocity[k] / speed;
            distance2 += d * d;
        }
        const double transverse2 = std::max(0., distance2 - longitudinal * longitudinal);
        weights[j] = std::exp(-transverse2 / (2. * width_ * width_))
            * .5 * std::erfc(-longitudinal / (std::sqrt(2.) * width_))
            / (2. * pi * width_ * width_);
    }
    return weights;
}
```

`longitudinal` is the projected distance to a nucleon center, `transverse2`
is its squared distance from the ray, and `weights[j]` is the integrated
contribution to the forward thickness. The explicit `primary_` check
implements exclusion of the hard-struck identity.

With `--nuclear-target off`, `NuclearGeometry` instead samples the hard vertex
uniformly inside a sphere, with \f$R_A=r_0A^{1/3}\f$,
\f$r_0=1.2\,\mathrm{fm}\f$, and
\f$\rho_0=0.17\,\mathrm{fm}^{-3}\f$. For a unit ray direction its remaining
distance and full chord thickness are

\f[
 \begin{aligned}
 L&=-\mathbf r_0\cdot\mathbf n+
   \sqrt{R_A^2-|\mathbf r_0|^2+(\mathbf r_0\cdot\mathbf n)^2},\\
 T_{\rm chord}&=2\rho_0
   \sqrt{R_A^2-|\mathbf r_0|^2+(\mathbf r_0\cdot\mathbf n)^2}.
 \end{aligned}
\f]

Legacy callers use both \f$\rho_0L\f$ and `compute_TA()`'s full chord; the
collision sampler rescales the tabulated \f$Q_s^2(T_{\rm chord})\f$ by
\f$\rho_0L/T_{\rm chord}\f$. Thus the legacy `Thickness` variable is not
the explicit target's forward-ray thickness. Legacy collision times also
exclude an initial interval \f$t_0=r_0\f$ equal to the radius coefficient in
natural units. Explicit-target mode uses primary-identity exclusion instead.

### Gluon TMD, collision rate and saturation scale

The medium kernel is implemented in `EHIJING::MultipleCollision` and
`alphas_PhiG()` in the bundled `external/pythia8/src/eHIJING.cc`. Its color
factors are

\f[
 \begin{gathered}
 N_c=3,\qquad d_A=N_c^2-1=8,\\
 C_F=\frac{N_c^2-1}{2N_c}=\frac43,\qquad C_A=N_c=3.
 \end{gathered}
\f]

For a parton in representation \f$R\f$, the differential eikonal collision
rate is related to an effective gluon TMD per nucleon by

\f[
 \begin{gathered}
 \frac{d\Gamma_R}{d^2\boldsymbol{\kappa}}
 =\frac{C_R}{d_A}\rho_N\,
   \frac{\alpha_s\phi_g(x_g,\kappa^2;Q_s^2)}{\kappa^2},\\
 x_g=\frac{x_B\kappa^2}{Q^2},\qquad
 \alpha_s\phi_g=K\frac{(1-x_g)^n x_g^\lambda}{\kappa^2+Q_s^2}.
 \end{gathered}
\f]

These are the rate and parametrization in paper Eqs. (9), (19). The function
stored in the code is the product \f$\alpha_s\phi_g\f$; multiplying it by
another running coupling would change the model. In this normalization
\f$K\f$ supplies the momentum-squared dimension needed for
\f$\alpha_s\phi_g\f$ to have units \f$\mathrm{GeV}^{-2}\f$; the CLI
accepts its numerical value in the source's GeV convention. The input
`--K` sets `eHIJING:TransverseMomentumDistributionConstant` in both shower
stages. The HERMES settings use \f$n=4\f$ and \f$\lambda=-0.25\f$;
the XML default for \f$\lambda\f$ is \f$-0.2\f$ when no override is read.
Positive \f$K\f$ normalizes the interaction strength, \f$n\f$ suppresses
the large-\f$x_g\f$ endpoint, and negative \f$\lambda\f$ increases the
small-\f$x_g\f$ density. This parametrization supplies no independent
small-\f$x\f$ evolution of the target.

**Implementation: gluon TMD.** `external/pythia8/src/eHIJING.cc`, `EHIJING::alphas_PhiG()`.

```cpp
double alphas_PhiG(double x, double kT2, double Qs2,
           double powerG, double lambdaG, double K){
    return K * std::pow(1. - x, powerG) * std::pow(x, lambdaG) / (kT2 + Qs2);
}
```

`x` corresponds to the exchanged gluon fraction, `kT2` to the collision
transfer squared, and `Qs2` to the screening scale. The return value already
contains the product of the coupling and the TMD.

Using \f$d^2\boldsymbol{\kappa}=\pi\,d\kappa^2\f$ for an azimuthally
symmetric kernel, the saturation scale is determined implicitly:

\f[
 \begin{aligned}
 Q_s^2&=\frac{C_A}{d_A}T_A
       \int^{Q^2/x_B}d^2\boldsymbol{\kappa}\,
         \alpha_s\phi_g(x_g,\kappa^2;Q_s^2)\\
 &=\frac{\pi C_A K T_A}{d_A}
       \int_{\kappa_{\min}^2}^{Q^2/x_B}
       d\kappa^2\,\frac{(1-x_B\kappa^2/Q^2)^n
          (x_B\kappa^2/Q^2)^\lambda}{\kappa^2+Q_s^2}.
 \end{aligned}
\f]

The first integral denotes the self-consistency relation of paper Eq. (20);
the second specifies the numerical lower bound used here:
\f$\kappa_{\min}^2=0.01\mu^2\f$, with \f$\mu=0.25\,\mathrm{GeV}\f$.
The code changes variable to \f$u=\ln(1+\kappa^2/Q_s^2)\f$ and solves the
resulting equation by bracketing and bisection. `Qs.dat` tabulates the result
on \f$(T_A,\ln x_B,\ln Q^2)\f$; explicit-target queries outside the table
solve the same equation directly. The solver can return roots below
\f$\mu^2\f$; the explicit-target high-scale **separation scale** is separately
floored at \f$\mu^2\f$. A failure to bracket a root warns and returns
\f$\mu^2\f$. This distinction matters when comparing to the paper's stated
minimum saturation scale.

The tables depend on \f$K,n,\lambda\f$ but their filenames do not encode
these parameters. A parameter scan must use matching table contents and
separate table directories as appropriate; loading an existing table does not
refit its physics parameters.

**Implementation: self-consistent saturation equation.** `external/pythia8/src/eHIJING.cc`, `MultipleCollision::Qs2_self_consistent_eq()`.

```cpp
double prefactor = piCAoverdA * TA;

// Integrand as a function of the log variable u
auto integrand_in_log_variable = [this, Qs2, xB, Q2](double u) {
    double q2 = Qs2 * (std::exp(u) - 1);
    double jacobian = Qs2 + q2;
    double xg = q2 * xB / Q2; // Gluon's light-cone momentum fraction
    return alphas_PhiG(xg, q2, Qs2, this->powerG_, this->lambdaG_, this->Kfactor_) * jacobian;
};
double error;

// Limits of integration in the original variable q^2
double q2_min = 0.01 * mu2;
double q2_max = Q2 / xB;

// Limits of integration in the log variable u
double u_min = std::log(1. + q2_min / Qs2);
double u_max = std::log(1. + q2_max / Qs2);

// Perform the integration using the 1D quadrature integrator
double res =  prefactor * quad_1d(integrand_in_log_variable, {u_min, u_max}, error);

return res - Qs2;
```

The integration variable `u` absorbs the screened denominator through its
Jacobian. Returning `res - Qs2` supplies the residual whose zero is located
by the bisection solver. The cutoff and upper kinematic limit are explicit
in `q2_min` and `q2_max`.

For a fixed path, the ideal dilute independent-scattering measure is

\f[
 \begin{gathered}
 \langle N_R\rangle=
 \frac{C_R T_A}{d_A}\int_{\kappa_{\min}^2}^{Q^2/x_B}
 \pi\,d\kappa^2\,
 \frac{\alpha_s\phi_g(x_B\kappa^2/Q^2,\kappa^2;Q_s^2)}{\kappa^2},\\
 \Pr(N_R)=\frac{\langle N_R\rangle^{N_R}}{N_R!}
             e^{-\langle N_R\rangle}.
 \end{gathered}
\f]

The implementation samples a decreasing-transfer sequence with proposal and
veto factors, storing a time, azimuth and, in explicit mode, target identity
for each collision. This proposal order is not time order. Its infrared
resolution is \f$\kappa_{\min}=0.1\mu=0.025\,\mathrm{GeV}\f$; the finite
collision count depends on that resolution. Event rejection by subsequent
exact recoil further conditions the accepted distribution, so the Poisson
expression is not a claim about the final accepted-event multiplicity.

**Implementation: sampling explicit-target collisions.** `external/pythia8/src/eHIJING.cc`, `MultipleCollision::sample_all_qt2_target()`.

```cpp
while (q2 > q2min) {
    const double lnr = std::log(std::max(1e-15, flat_gen(gen)));
    const double xg = q2 / q2max;
    q2 *= std::pow(1. + (lambdaG_ - 1.) * lnr * q2 / tildeTA
        * std::pow(xg, -lambdaG_), 1. / (lambdaG_ - 1.));
    // Keep the original eHIJING veto kernel and momentum-transfer proposal.
    if (flat_gen(gen) < std::pow(1. - xg, powerG_) * q2 / (q2 + qs2)) {
        double time;
        int id = target.sample_collision(start, velocity,
            [this]() { return flat_gen(gen); }, time);
        q2_list.push_back(q2); t_list.push_back(time); targets.push_back(id);
        phi_list.push_back(2. * M_PI * flat_gen(gen));
    }
}
return q2_list.size();
```

`q2` evolves downward through transfer proposals; accepted entries receive
a sampled position/time, identity and azimuth. This excerpt preserves the
actual ordering of the proposal and veto operations, rather than replacing
them with a schematic independent draw from the rate.

### Transport coefficient and transverse momentum broadening

The second moment of the rate defines the local broadening coefficient:

\f[
 \hat q_R(s)=\int d^2\boldsymbol{\kappa}\,
       \kappa^2\frac{d\Gamma_R(s)}{d^2\boldsymbol{\kappa}}
 =\frac{C_R}{d_A}\rho_N(s)
       \int d^2\boldsymbol{\kappa}\,\alpha_s\phi_g.
\f]

For the same integration bounds and path-dependent \f$Q_s\f$ used in the
self-consistency equation,

\f[
 \int ds\,\hat q_R(s)=\frac{C_R}{C_A}Q_s^2,\qquad
 \overline{\hat q}_R=\frac{C_R}{C_A}\frac{Q_s^2}{L}.
\f]

The second expression is a path average for a specified length \f$L\f$;
there is no unique geometric edge length for a Gaussian target. With the
present color factors, gluon broadening is \f$C_A/C_F=9/4\f$ times quark
broadening for an otherwise identical path and kernel. The second moment is
less sensitive to unresolved soft transfers than the collision count because
its \f$\kappa^2\f$ weight cancels the rate's explicit
\f$1/\kappa^2\f$. In this TMD model the zero-cutoff second moment converges
for \f$\lambda>-1\f$.

These equations define the transport coefficient, rather than promising that
every similarly named legacy helper returns it with this normalization.
For example, `qhatA()` includes `rho0 * Qs2 / TA`, while `qhatF()` in this
checkout returns `(CF/CA) * Qs2 / TA` without that density factor. Derive
reported \f$\hat q\f$ from a dimensionally consistent rate or path integral.
The observable hadron broadening additionally depends on fragmentation,
selection, radiative recoil and its correlation with energy loss.

**Implementation: legacy transport helpers.** `external/pythia8/include/eHIJING/eHIJING.h`, `MultipleCollision::qhatA() and qhatF()`.

```cpp
double qhatA(double TA, double xB, double Q2) {
    return rho0*Qs2(TA, xB, Q2)/TA;
};
// compute the qhat of a quark
double qhatF(double TA, double xB, double Q2) {
    return CF/CA*Qs2(TA, xB, Q2)/TA;
};
```

The missing `rho0` factor in `qhatF()` is visible directly. These helpers
are included to explain the normalization caveat above; the mathematical
definition of the transport coefficient remains the second moment of the
collision rate.

### Radiation kernels: vacuum, HT and GHT

The medium correction acts on gluon emission from energetic light quarks,
antiquarks and gluons. Define a vacuum density per \f$dz\,d\ell^2\f$ by

\f[
 \mathcal V_{a}(z,\ell^2)
 =\frac{\alpha_s(\ell^2)}{2\pi}\frac{P_{ga}(z)}{\ell^2},\qquad
 P_{gq}(z)=C_F\frac{1+(1-z)^2}{z},
\f]

\f[
 P_{gg}(z)=C_A\left[
   \frac{1+(1-z)^3}{z}+\frac{1+z^3}{1-z}\right].
\f]

Here \f$P\f$ excludes the coupling, unlike some conventions in the paper.
The two gluon terms correspond to PYTHIA's two color-dipole contributions;
the low-scale sampler uses the emitted-gluon branch and its finite-\f$z\f$
correction. These are real-emission kernels. Virtual corrections and
no-emission probabilities enter through the Sudakov construction, rather
than through a separately sampled plus distribution.

In the local medium routines the coupling is

\f[
 \alpha_s(\mathcal Q^2)=
 \frac{2\pi}{b_0\max[1,\ln(\mathcal Q^2/\mu^2)]},\qquad
 b_0=\frac{11N_c-2n_f}{6}=\frac92,\qquad n_f=3.
\f]

The implementation uses `2.71828` for the freeze boundary
\f$\mathcal Q^2/\mu^2\simeq e\f$. This is the medium routine's coupling,
with a potentially large frozen value; it does not override all independent
PYTHIA shower-coupling settings.

**Implementation: medium running coupling.** `external/pythia8/src/eHIJING.cc`, `EHIJING::alphas()`.

```cpp
double alphas(double Q2) {
    double scale_ratio_squared = Q2 / mu2;
    if (scale_ratio_squared < 2.71828) {
        return twoPioverb0;
    } else {
        return twoPioverb0 / std::log(scale_ratio_squared);
    }
}
```

`twoPioverb0` is the prefactor in the coupling equation above. The low-scale
adapter evaluates this function at the stored separation scale and holds
that value fixed during its radiation sampling.

The following stochastic kernels are useful mathematical representations of
the collision-conditioned modification. Let the sampled collision history
of a parent in representation \f$R\f$ be
\f$\{(t_i,\boldsymbol{\kappa}_i)\}_{i=1}^{N}\f$ with parent energy
\f$E\f$ in the nuclear rest frame. The interference function is

\f[
 \Phi(u)=1-\cos u,\qquad
 \Phi(u\ll1)=\frac{u^2}{2}+\mathcal O(u^4).
\f]

For `--medium-modification-mode 0`, the collinear higher-twist (HT) form is

\f[
 \begin{gathered}
 \frac{dN_a^{\rm HT}}{dz\,d\ell^2}
 =\mathcal V_a(z,\ell^2)
 \left[1+\frac{C_A}{C_R}\sum_{i=1}^{N}
    \Theta(\ell^2-\kappa_i^2)\frac{2\kappa_i^2}{\ell^2}
    \Phi\!\left(\frac{t_i}{\tau_{\rm HT}}\right)\right],\\
 \tau_{\rm HT}=\frac{2z(1-z)E}{\ell^2}.
 \end{gathered}
\f]

The step function implements the restriction of the collinear expansion to
transfers below the radiation scale. In a continuum notation its medium term
is proportional to
\f$2\mathcal V_a\int ds\,\hat q_A^{\rm rad}(s,\ell^2)
\Phi(s/\tau_{\rm HT})/\ell^2\f$, where

\f[
 \hat q_A^{\rm rad}(s,\ell^2)
 =\rho_N(s)\frac{C_A}{d_A}
   \int_{\kappa_{\min}^2}^{\ell^2}
      \pi\,d\kappa^2\,\alpha_s\phi_g.
\f]

This radiative coefficient has upper limit \f$\ell^2\f$, whereas the
collisional coefficient above uses \f$Q^2/x_B\f$. They cannot generally
be interchanged or fitted as though they were the same parameter.

**Implementation: HT collision-conditioned weight.** `external/pythia8/src/eHIJING.cc`, `eHIJING::sample_z_stochastic(), mode 0`.

```cpp
z = zmin*std::pow(zmax/zmin, flat_gen(gen));
double tauf = 2.*z*(1.0-z)*E/kt2;
double w = 1.0, wmax = 1.0;
for (int i=0; i<Ncolls; i++) {
    double q2 = qt2s[i], t = ts[i];
    if (q2 > kt2) continue;
    w += 2.*q2/kt2 * CA / CR * (1.-cos(t/tauf));
    wmax += 4.*q2/kt2 * CA / CR;
}
acceptance = w/wmax;
```

Within the HT trial loop, `kt2` represents the radiation scale and
`qt2s[i]` the collision transfer squared. The `q2 > kt2` veto implements
the collinear restriction. The initial `w = 1.0` retains the vacuum term;
each accepted collision contributes its color factor and LPM phase.
`wmax` is the overestimate used to accept or reject the trial.

For `--medium-modification-mode 1`, the generalized higher-twist (GHT)
soft-gluon kernel retains the transverse transfer in the propagator and phase:

\f[
 \frac{dN_a^{\rm GHT}}{dz\,d\ell^2}
 =\mathcal V_a(z,\ell^2)
 \left[1+\frac{C_A}{C_R}\sum_{i=1}^{N}
 \frac{2\boldsymbol{\kappa}_i\cdot\boldsymbol{\ell}}
      {|\boldsymbol{\ell}-\boldsymbol{\kappa}_i|^2}
 \Phi\!\left(
  \frac{|\boldsymbol{\ell}-\boldsymbol{\kappa}_i|^2t_i}
       {2z(1-z)E}\right)\right].
\f]

This collision-history representation corresponds to paper Eqs. (44)--(45);
the high-scale implementation performs angular integration/approximation
using `GHT.dat`, Bessel functions and veto weights. It does not sample every
term of the unrestricted finite-\f$z\f$ generalized twist-four expression.
In particular, retaining \f$1-z\f$ and finite-\f$z\f$ splitting numerators
is not a derivation beyond the soft-gluon approximation. An individual angular
interference term can have either sign; it should not be interpreted as an
independent positive emission probability.

In GHT the relevant formation denominator is
\f$|\boldsymbol{\ell}-\boldsymbol{\kappa}_i|^2\f$; in HT it is
\f$\ell^2\f$. The apparent GHT collinear singularity is regulated by the
phase: as \f$D=|\boldsymbol{\ell}-\boldsymbol{\kappa}_i|^2\to0\f$,
\f$\Phi(Dt_i/[2z(1-z)E])/D\propto D\f$. For formation times much longer
than the sampled collision times, the LPM factor suppresses medium radiation.
For rapidly varying phases, its average tends to one. With \f$N=0\f$, both
forms return to the vacuum kernel. HT and GHT are alternative approximations;
their corrections are not added together.

**Implementation: GHT angular interference.** `external/pythia8/src/eHIJING.cc`, `eHIJING::compute_GHT_Angular_Table()`.

```cpp
double eHIJING::compute_GHT_Angular_Table(double X, double Y) {
    double A = Y / (1.0 - X);
    double B = Y * X /(1.0 - X);
    auto dfdphi = [X, Y, A, B](double phi) {
        double cosphi = std::cos(phi);
        double xcphi = X * cosphi;
        return xcphi / (1.0 - xcphi) * (1.0 - std::cos(A - B * cosphi));
    };
    double error;
    double result = quad_1d(dfdphi, {0., M_PI}, error);
    result /= M_PI;
    return result;
}
```

For a collision in `sample_z_stochastic()`, the table variables are
`X = 2*sqrt(kt2*q2)/(kt2+q2)` and
`Y = (kt2+q2)*t*(1-X)/(2*z*(1-z)*E)`. Thus `A - B*cosphi` is the
transfer-dependent formation phase, while `xcphi/(1-xcphi)` is the
angular propagator factor. The division by pi normalizes the angular
average; the calling sampler also has approximation branches described
above.

### Matching the high- and low-scale showers

The executable combines two emission samplers. In the modified PYTHIA
`TimeShower`, a schematic conditional no-emission probability between squared
transverse scales \f$u_1>u_2\f$ is

\f[
 \begin{gathered}
 \mathcal R_a(z,u)=\mathcal V_a(z,u)+\mathcal M_a(z,u)
                    \Theta(u-Q_{\rm sep}^2),\\
 \Delta_a(u_2,u_1)=\exp\!\left[
 -\int_{u_2}^{u_1}du\int_{z_{\min}(u)}^{z_{\max}(u)}dz\,
 \mathcal R_a(z,u)\right].
 \end{gathered}
\f]

where \f$\mathcal M_a\f$ is the selected medium correction. Trial scales
are drawn using this Sudakov structure and accepted with the relevant vetoes
and PYTHIA kinematic checks. Vacuum showering continues to PYTHIA's cutoff;
the high-scale medium term is disabled below \f$Q_{\rm sep}^2\f$.
The explicit-target implementation evaluates a ray-dependent saturation scale,
sets \f$Q_{\rm sep}^2=\max(Q_s^2,\mu^2)\f$, and stores it in the event's
`SeparationScale()`. This value can be updated by shower evolution; the
low-scale adapter consumes the stored event value, rather than an additional
fixed user-chosen matching parameter.

**Implementation: storing the separation scale.** `external/pythia8/src/TimeShower.cc`, `TimeShower::pT2nextQCD()`.

```cpp
if (target) {
  Thickness = target->thickness(start, velocity);
  Qs2 = Thickness > 1e-12 ? eHIJING_Gen->Qs2_target(Thickness, xB, Q20) : EHIJING::mu2;
  Qs2 = std::max(Qs2, EHIJING::mu2);
}
event.setSeparationScale(Qs2);
```

`Thickness` is evaluated along the current explicit-target ray. The
`max` operation supplies the nonperturbative floor, and
`event.setSeparationScale()` is the handoff to the low-scale adapter.

`Modified_FF::sample_ff_partons()` then acts on final light partons with sampled
collisions. It skips heavy quarks and diquarks and samples additional medium
gluons with increasing formation variable. In the source, `kt2` is the
denominator scale inferred from

\f[
 \tau=\frac{2z(1-z)E}{k_T^2},\qquad
 \tau_{\min}=\frac{1}{\sqrt{Q_{\rm sep}^2}},\qquad
 \tau_{\max}=\frac{E}{\mu^2}.
\f]

This starting value follows the local adapter, and should not be replaced by
\f$1/Q\f$ simply because a high-scale starting time is discussed in the
paper. The adapter holds \f$\alpha_s(Q_{\rm sep}^2)\f$ fixed during this
stage. Its energy resolution is `min_parton_energy_gev = 0.2`, with
\f$z_{\min}=\min(0.2\,\mathrm{GeV}/E,0.4)\f$,
\f$z_{\max}=1-z_{\min}\f$, and a radiation loop requiring
\f$E>0.4\,\mathrm{GeV}\f$. Finite-\f$z\f$ veto factors are
\f$[1+(1-z)^2]/2\f$ for a quark and
\f$[1+(1-z)^3]/2\f$ for a gluon.

**Implementation: low-scale formation ordering.** `src/modified_ff.cpp`, `Modified_FF::sample_ff_partons()`.

```cpp
double formation_time_max = particle.e() / EHIJING::mu2;
// Start from the minimum formation time ~ 1/Qs, which is the separation scale;
double formation_time_min = 1.0 / std::sqrt(pythia.event.SeparationScale());
double formation_time = formation_time_min;

// ...

if (mode_ == 0) {
    double coeff = alpha_bar * maxdiffz * 4. * sumq2 / 2. / particle.e();
    formation_time = formation_time + std::log(1. / r) / coeff;
// If mode = 1, use GHT
} else if (mode_ == 1) {
    double invrpower = alpha_bar * log_z_max * 4. * n_collisions;
    double step_factor = std::pow(1. / r, 1. / invrpower);
    formation_time = formation_time * step_factor;
}
```

The first excerpt initializes the adapter's formation-variable interval.
The second shows how an HT trial advances additively and a GHT trial
multiplicatively before the subsequent LPM and kinematic vetoes. Neither
step by itself assigns a measured formation vertex.

In HT, the adapter keeps transfers \f$\kappa_i^2<k_T^2\f$, samples the
emission azimuth uniformly, and uses `lt2 = kt2`. In GHT it uses supporting
collisions with \f$k_T^2<\kappa_i^2\f$, samples a correlated azimuth, and
constructs

\f[
 \ell_T^2=k_T^2+\kappa_i^2
                 +2\sqrt{k_T^2\kappa_i^2}\cos\Delta\phi.
\f]

The matching veto is imposed on `lt2`:
\f$\ell_T^2\leq Q_{\rm sep}^2\f$. The generated gluon's momentum is
constructed with transverse magnitude \f$k_T\f$ and energy \f$zE\f$,
requiring \f$k_T\leq zE\f$. Thus `kt2`, `lt2` and the PYTHIA evolution
scale should not be identified indiscriminately. This is a matching
prescription with vetoes, not a proof that all higher-order overlap terms
have been subtracted.

**Implementation: vetoing overlap and unphysical emission kinematics.** `src/modified_ff.cpp`, `Modified_FF::sample_ff_partons()`.

```cpp
if (lt2 > kt2max_now) {
    continue;
}

// Now, there is a radiation,
// Hard parton splits into p -> p-k & k
double kt = std::sqrt(kt2), k0 = z * particle.e();
if (kt > k0) {
    continue;
}
```

`kt2max_now` was read from `event.SeparationScale()` before the parton
loop. The overlap veto uses `lt2`, while the on-shell emitted-gluon
construction uses `kt = sqrt(kt2)` and requires sufficient gluon energy.

In explicit-target mode a new gluon samples its own subsequent collisions
from the remaining density, starting at its formation-variable displacement
along the parent's ray. Its collision times receive the corresponding offset.
The newly appended gluons undergo elastic recoil but are not recursively
passed through another full low-scale radiation cascade in the same call.
Radiative sampling and exact recoil are separate operations; ordering the
elastic updates in time does not turn the entire algorithm into a
time-resolved coupled transport equation.

### Elastic recoil, color exchange and nuclear remnants

In explicit-target mode each struck nonprimary nucleon owns one shared,
on-shell remnant momentum, initially \f$K_j^\mu=(m_N,\mathbf0)\f$ in the
nuclear rest frame. Sampled elastic hits from all processed partons are sorted
by their rest-frame times. Each hit updates that same momentum, including
when several partons hit the same identity.

The solver `ehijing::recoil_two_body()` retains the sampled transverse
transfer and determines its longitudinal and energy components from

\f[
 \begin{gathered}
 p'=p+\delta,\qquad K_j'=K_j-\delta,\\
 \boldsymbol{\delta}_\perp=
    \sqrt{\kappa^2}(\cos\phi,\sin\phi),\\
 (p')^2=p^2,\qquad (K_j')^2=K_j^2.
 \end{gathered}
\f]

The transverse plane is relative to the incoming jet direction at that update.
In this axis frame define \f$a=p^0+K_j^0\f$,
\f$b=p_z+K_{jz}\f$ and
\f$c=\mathbf K_{j\perp}\cdot\boldsymbol{\delta}_\perp\f$.
The on-shell constraints imply

\f[
 \begin{gathered}
 \delta^0=\frac{b\delta_z+c}{a},\\
 \left(\frac{b^2}{a^2}-1\right)\delta_z^2
 +2\left(\frac{bc}{a^2}+\frac{p^0b}{a}-p_z\right)\delta_z
 +\frac{c^2}{a^2}+\frac{2p^0c}{a}-\kappa^2=0.
 \end{gathered}
\f]

**Implementation: on-shell elastic recoil equation.** `include/nuclear_recoil.hpp`, `ehijing::recoil_two_body()`.

```cpp
const double qx = std::sqrt(qt2) * std::cos(phi), qy = std::sqrt(qt2) * std::sin(phi);
const double a = p.e() + k.e(), b = p.pz() + k.pz();
const double c = k.px() * qx + k.py() * qy;
const double A = (b/a - 1.) * (b/a + 1.);
const double B = 2. * (b*c/(a*a) + p.e()*b/a - p.pz());
const double C = c*c/(a*a) + 2.*p.e()*c/a - qt2;
const double D = B*B - 4.*A*C;
if (!(D >= 0.) || !std::isfinite(D) || A == 0.) return false;
const double stable = -.5 * (B + std::copysign(std::sqrt(D), B));
const double roots[2] = {stable/A, stable == 0. ? 0. : C/stable};
```

The jet has already been rotated onto the longitudinal axis in this
excerpt; `k` is the whole target-remnant four-momentum. `A`, `B` and `C`
are the coefficients of the longitudinal-transfer quadratic above. A
negative discriminant rejects the collision, and the remaining function
checks positive energies, the forward branch and both invariant masses.

Among physical roots the code selects the forward solution with smallest
\f$|\delta_z|\f$, connected to the zero-transfer branch, and verifies the
masses. Four-momentum is conserved at each accepted elastic update; energy
transferred to the target therefore gives collisional energy loss without an
additional independent drag term. A transfer without a physical forward
solution rejects the whole generation attempt. It is not clipped or redrawn
within that event. This conditions the accepted sample and makes the rejection
rate part of assessing the implementation's phenomenology.

After the hits, each touched nonprimary nucleon is split **once** into a
constituent quark and diquark in its own final rest frame. Their opposite
momentum magnitude is

\f[
 \begin{gathered}
 p_* =\frac{\sqrt{[m_N^2-(m_q+m_{qq})^2]
                         [m_N^2-(m_q-m_{qq})^2]}}{2m_N},\\
 m_q=0.33\,\mathrm{GeV},\quad m_{qq}=0.57933\,\mathrm{GeV}.
 \end{gathered}
\f]

An isotropic direction is sampled and both daughter four-momenta are boosted
with the final remnant momentum. The fixed identity determines whether the
valence content is proton-like or neutron-like; the source uses the
corresponding quark/diquark choices with a 3:1 scalar/vector ratio for an
\f$ud\f$ pair. The assigned displacement is the last sampled collision
point. Repeated hits exchange existing color lines, and medium gluons are
inserted into a supporting remnant's string. Neither an additional hit nor a
medium emission creates an additional target baryon in this mode.

**Implementation: shared remnant updates.** `src/modified_ff.cpp`, `Modified_FF::sample_ff_partons()`.

```cpp
std::stable_sort(kicks.begin(), kicks.end(), [](const Kick& a,const Kick& b) { return a.time < b.time; });

// ...

auto& jet = hit.appended ? new_particles.at(hit.source) : pythia.event[hit.source];
auto& r = remnants[hit.target];
auto momentum = jet.p();
if (!ehijing::recoil_two_body(momentum,r.p,hit.qt2,hit.phi)) {
    std::cerr << "WARNING: rejecting event: no physical recoil; target_id=" << hit.target
              << " time_GeV_inv=" << hit.time << " qt2_GeV2=" << hit.qt2 << '\n';
    return false;
}
ensure_remnant(hit.target,hit.appended ? new_sources.at(hit.source) : hit.source,hit.vertex);
jet.p(momentum); ++target->nucleons()[hit.target].hits;
r.vertex = hit.vertex; // Final split starts at the last sampled interaction.
// Exchange an existing color line; repeated hits never create a new pair.
if (jet.col() && (!jet.acol() || uniform_dist_(rng_) < .5)) {
    const int old = jet.col(); jet.col(r.col); r.col = old;
} else {
    const int old = jet.acol(); jet.acol(r.acol); r.acol = old;
}
```

The first line orders the complete elastic-hit list. Inside the hit loop,
`remnants[hit.target]` refers to one persistent nucleon identity. Successful
recoil updates its momentum in place; the color exchange shown here reuses its
existing string endpoints, including after repeated hits.

The primary remnant is handled with the hard-process system. Untouched
identities become spectator nucleons. For supported neutral-current events,
the final non-spectator baryon number \f$B_{\rm final}\f$ and charge
\f$Q_{\rm final}\f$ imply

\f[
 N_p^{\rm spec}=Z-Q_{\rm final},\qquad
 N_n^{\rm spec}=A-B_{\rm final}-N_p^{\rm spec}.
\f]

The scattered DIS lepton is excluded from this hadronic budget; signed final
quantum numbers include antiparticles. In explicit mode the required counts
must match the actual untouched identities or the attempt is rejected.
Spectators carry their sampled positions and initially stationary rest-frame
momenta, boosted to the output frame. There is no binding-energy, evaporation,
coalescence or excitation/decay model for the residual nucleus.

The legacy `--nuclear-target off` branch uses independent constituent recoil
constructions and newly sampled spectator counts. It does not share a target
identity across hits and does not implement the exact whole-remnant recoil
prescription above. The paper's optional simplified collisional-loss formula
should not be interpreted as a separate selectable exact-drag model in the
current ELECTRA CLI.

Exact elastic conservation also does not imply exact global conservation of
the full algorithm. In the existing low-scale radiative step the gluon's
three-momentum is subtracted from the parent and the parent's energy is then
restored to its mass shell independently:

\f[
 \mathbf p'=\mathbf p-\mathbf p_g,\qquad
 E'=\sqrt{|\mathbf p'|^2+m^2},\qquad
 \Delta E_{\rm rad}=E'+E_g-E.
\f]

**Implementation: parent energy after low-scale radiation.** `src/modified_ff.cpp`, `Modified_FF::sample_ff_partons()`.

```cpp
particle.p(particle.p() - kmu);
particle.e(std::sqrt(particle.pAbs2() + particle.m2()));
kmu.e(std::sqrt(kmu.pAbs2()));
```

`kmu` is the emitted gluon. The parent spatial momentum is reduced, then
its energy is independently recomputed on shell. This is the exact code
operation behind the radiative energy-balance limitation in the equation
above.

In general \f$\Delta E_{\rm rad}\ne0\f$. This remaining radiative
approximation is not repaired by the subsequent exact two-body elastic solver.
The explicit-target changes have synthetic tests; the README records that
they have not yet been validated with a physics production sample.

### Lund hadronization and the meaning of formation coordinates

`Hadronizer` imports final partons and their color connections into a separate
PYTHIA instance with `ProcessLevel:all = off`, then runs its string
fragmentation. The vacuum Lund model supplies the nonperturbative conversion
of the color system to hadrons. For light-flavor fragmentation the symmetric
longitudinal law has the form

\f[
 f(z_L)\propto\frac{(1-z_L)^{a_L}}{z_L}
       \exp\!\left(-\frac{b_Lm_{\perp h}^2}{z_L}\right),\qquad
 m_{\perp h}^2=m_h^2+p_{\perp h}^2.
\f]

Here \f$z_L\f$ is a string light-cone fraction, distinct from the parton
emission fraction \f$z\f$ and the SIDIS fraction \f$z_h\f$.
The parameters and flavor production probabilities come from the
hadronization settings and PYTHIA tune; this equation is not a complete
description of massive endpoints, junctions or low-mass string collapse.
The local HERMES settings use `Tune:pp = 19`,
`StringFragmentation:stopMass = 0.0` and `HadronLevel:mStringMin = 0.5`.
The dedicated hadronizer also synchronizes selected particle masses with
ELECTRA/SMASH. Decays follow its settings; the supplied bundle disables
hadron decays for downstream handling.

With `Fragmentation:setVertices = on` and `HadronVertex:mode = 2`, the saved
hadron vertex follows the late formation prescription. For an ordinary light
string segment without endpoint corrections or smearing, its local form is

\f[
 x_{\rm form}^{\mu}=x_{\rm origin}^{\mu}
       +\frac{x_i^\mu+x_{i+1}^\mu}{2}
       +\frac{p_h^\mu}{2\kappa_{\rm str}},
\f]

where \f$x_i^\mu,x_{i+1}^\mu\f$ are relative breakup vertices in fm and
\f$\kappa_{\rm str}\f$ is the string tension in GeV/fm. It is unrelated
to a collision transfer \f$\boldsymbol{\kappa}\f$. The actual PYTHIA
implementation contains massive-endpoint factors, special string branches
and collapse prescriptions; its recorded production references and breakup
pairs should be used for diagnostics instead of extrapolating this simple
formula to all hadrons.

**Implementation: late Lund formation prescription.** `external/pythia8/src/StringFragmentation.cc`, `StringFragmentation::setHadronVertices()`.

```cpp
if (hadronVertex == 0) prodPoints = middlePoint;
else if (hadronVertex == 1 || hadronVertex == 2)
  prodPoints = middlePoint + 0.5 * redOsc * pHad / kappaVtx;
```

`middlePoint` is the average of the local neighboring breakup vertices,
`pHad` the hadron momentum and `kappaVtx` the string tension. `redOsc` is
one for the light-endpoint limit displayed above; massive endpoints can
modify it. The saved absolute vertex later adds `prodPoints * FM2MM` to
the inherited origin.

The hard origin \f$(\mathbf r_0,0)\f$ and relative medium displacements
are converted from natural units to PYTHIA's mm units. In explicit-target
mode they are boosted as four-vectors to the generator frame. The hadronizer
adds the nuclear origin exactly once. Import retains final partons with
status 23 but does not retain their parent links, so its ancestry traversal
does **not** reconstruct a high-scale shower delay. Assigned medium
displacements remain present, and Lund breakup vertices are still generated.
A sampled formation variable is therefore not a complete measured shower
formation event.

**Implementation: importing final partons.** `src/hadronizer.cpp`, `Hadronizer::hadronize()`.

```cpp
pythia.event.append(particle.id(),
                    constants::pythia::hard_process_outgoing_status,
                    particle.col(),
                    particle.acol(),
                    particle.px(),
                    particle.py(),
                    particle.pz(),
                    particle.e(),
                    particle.m());

auto& appended = pythia.event.back();
if (history) source_indices.push_back(i);
// Input vertices are displacements from the hard interaction in mm.
appended.vProd(particle.vProd());
```

This `append` overload supplies momentum, mass, status and color, with no
parent indices. The following `vProd` assignment preserves the imported
relative displacement in mm. Together these explain why the medium
displacement survives while the source shower ancestry does not.

For OSCAR output, the writer converts the saved full formation vertex to fm
and extrapolates freely to its \f$t=0\f$ input surface:

\f[
 \mathbf x(0)=\mathbf x_{\rm form}
          -\frac{\mathbf p_h}{E_h}\,t_{\rm form}.
\f]

**Implementation: freestreaming to the output surface.** `src/event_output.cpp`, `write_final_hadrons()` (called by `write_event_output()`).

```cpp
if (std::isfinite(form_time) && std::isfinite(e) && e > 0.0) {
    const double dt_to_zero = -form_time;
    x_hadron += p.px() / e * dt_to_zero;
    y_hadron += p.py() / e * dt_to_zero;
    z_hadron += p.pz() / e * dt_to_zero;
}
```

`form_time` and the three position variables have already been converted
to fm-based output units. The factor `p/E` is dimensionless; no additional
boost factor or formation-time conversion is introduced here.

It separately writes `prod_time` and `form_time` as coordinate times in that
same frame. These are not proper-time intervals to be multiplied by a new
Lorentz factor. The appended dimensionless `prod_time_status` is `0` for
unavailable production information (a zero placeholder), `1` for an available
production reference including zero, and `2` for not applicable spectator
nucleons. SMASH uses zero for status `0` and warns on stderr; this is an
explicit fallback, not a reconstructed production event. The extrapolated
\f$t=0\f$ position is not a new physical
interaction or formation vertex. The writer adds no spatial smearing.
eHIJING stops at this output interface; any subsequent hadron rescattering or
formation-dependent interaction profile is part of SMASH's configured model.
The optional HepMC history likewise records observations and available
attributes, rather than supplying missing shower worldlines.

### Regime of applicability and interpreting observables

The paper motivates a localized hard vertex and fragmentation outside the
target by the hierarchies

\f[
 \frac{1}{Q}\ll L,\qquad
 \frac{\nu}{Q^2}=\frac{1}{2m_Nx_B}\ll L,\qquad
 \tau_h\sim\frac{z_h\nu}{\Lambda_{\rm QCD}^2}\gg L.
\f]

They are physical scope estimates, not automatic cuts or replacements for
the generated Lund vertices. The code does not dynamically switch to coherent
photon-dipole scattering when \f$x_B\f$ becomes small, nor to a heavy-quark
medium-radiation calculation when a heavy flavor is produced. Fitted nPDF
shadowing does not add a coherent multi-nucleon hard-scattering geometry.
The medium kernels are dilute, eikonal, soft-gluon approximations with repeated
sampling; they do not implement an unrestricted dense-medium resummation.

Useful SIDIS observables include the multiplicity ratio and hadron broadening,
evaluated with identical binning and selection in both samples:

\f[
 \begin{gathered}
 R_A^h(z_h,p_{T,h}^2)=
 \frac{\left[(1/N_{\rm DIS})\,dN_h/(dz_h\,dp_{T,h}^2)\right]_A}
      {\left[(1/N_{\rm DIS})\,dN_h/(dz_h\,dp_{T,h}^2)\right]_D},\\
 \Delta\langle p_{T,h}^2\rangle=
   \langle p_{T,h}^2\rangle_A-\langle p_{T,h}^2\rangle_D.
 \end{gathered}
\f]

Here \f$D\f$ denotes the chosen deuterium baseline. These definitions do
not assert an existing histogram or analysis in the generator. Both radiation
and broadening depend on the same sampled collision history. Selecting a
large-\f$z_h\f$ hadron also selects its energy-loss history, so the observed
broadening is a conditional quantity and cannot generally be equated to
\f$z_h^2(C_F/C_A)Q_s^2\f$ without additional approximations. Hadronization,
the nPDF flavor mixture, recoil prescription, acceptance and any subsequent
SMASH interactions also enter an ELECTRA-level comparison.

**Implementation: configured kinematic selection.** `src/dis_kinematics.cpp`, `is_valid_dis_event()`.

```cpp
bool is_valid_dis_event(const DISKinematics& kin, const DISCuts& cuts) {
    return passes_range(kin.y, cuts.y_min, cuts.y_max) &&
           passes_range(kin.bjorken_x, cuts.bjorken_x_min, cuts.bjorken_x_max) &&
           passes_range(kin.nu, cuts.nu_min_gev, cuts.nu_max_gev) &&
           passes_range(kin.Q2, cuts.q2_min_gev2, cuts.q2_max_gev2) &&
           passes_range(kin.W2, cuts.w2_min_gev2, cuts.w2_max_gev2);
}
```

The selection checks only the configured ranges in reconstructed DIS
variables. The localization and hadron-formation hierarchies discussed
above do not appear as additional implicit acceptance conditions.

For reproducible interpretation, retain the HT/GHT mode, \f$K,n,\lambda\f$,
the table directory, geometry parameters, hard/shower and hadronization
settings, DIS cuts and component revision. The main implementation locations
are `src/dis_kinematics.cpp`, `src/modified_ff.cpp`,
`include/nuclear_recoil.hpp`, `src/hadronizer.cpp` and
`src/event_output.cpp`, together with the bundled PYTHIA files
`src/eHIJING.cc`, `src/TimeShower.cc`, `src/PartonLevel.cc` and
`include/eHIJING/NuclearTarget.h`. The physics equations should be checked
against those sources when changing the model; this section records their
present prescriptions and introduces no new physics settings.

## Runtime Inputs

The generator is configured through command-line arguments and input files
rather than recompilation-sensitive constants.

- Hard-process PYTHIA settings are loaded from `RunConfig::hard_process_config_path`.
- Hadronization PYTHIA settings are loaded from `RunConfig::hadronization_config_path`.
- DIS event-selection cuts are loaded from `RunConfig::dis_cuts_config_path`.
- Tabulated eHIJING data are read from `RunConfig::tabulation_path`.

Keep the three settings files together under `input/ehijing/experiments/<experiment>/`
using the names `hard_process.setting`, `hadronization.setting`, and
`dis_cuts.setting`. The Slurm submit script selects the directory with
the required `EXPERIMENT` variable; each config path can still be overridden
independently when needed. Use a separate experiment directory when detector
acceptance or reconstruction assumptions require different settings, including
hadronization settings. Add bundles only for experiments with configured files;
the repository currently supplies the HERMES bundle.

The DIS cuts file is parsed by `load_dis_cuts()` and supports the optional keys
`yMin`, `yMax`, `xBMin`, `xBMax`, `nuMin`, `nuMax`, `Q2Min`, `Q2Max`, `W2Min`,
and `W2Max`.

## Public API Map

- `RunConfig` and `parse_args()` define the executable runtime interface.
- `DISKinematics`, `DISCuts`, `compute_dis_kinematics()`, and `trigger()` define
  DIS reconstruction and event selection.
- `Modified_FF` owns the low-\f$ Q^2 \f$ medium-modification sampler.
- `Hadronizer` owns the dedicated PYTHIA instance used for string fragmentation.
- `write_event_headers()` and `write_event_output()` write OSCAR2013 particle
  records and per-event JSON metadata.
- `make_event_paths()` and related helpers define deterministic sharded output
  paths for large productions.
- `ehijing::constants` centralizes compile-time conventions that are not meant
  to change run by run.

## Conventions

The code relies on a few explicit PYTHIA/eHIJING conventions:

- PYTHIA DIS event-record indices are defined in `ehijing::constants::pythia`.
- PDG Monte Carlo particle IDs live in `ehijing::constants::pdg` and use the
  `_id` suffix.
- Masses and numerical cutoffs that are still compile-time constants are grouped
  in `ehijing::constants::mass`, `ehijing::constants::numeric`, and related
  namespaces.

See @ref unit_conventions for spacetime units, conversion factors, interface
requirements, and remaining review items.

## Generating The HTML

From `src/ehijing`, either run Doxygen directly:

```bash
doxygen Doxyfile
```

or use the CMake target:

```bash
cmake -S . -B build
cmake --build build --target docs
```

If you only want the documentation target and do not need to configure the
generator executable, use:

```bash
cmake -S . -B build-docs -DEHIJING_BUILD_DOCS_ONLY=ON
cmake --build build-docs --target docs
```

The generated entry point is:

```text
docs/doxygen/html/index.html
```
