# Learning the physics of eHIJING

This tutorial is for a new graduate student who knows special relativity and
has encountered quarks and gluons, but has not studied deep-inelastic
scattering or parton propagation through nuclei. The goal is to understand
what an eHIJING event represents, why its equations have their particular
form, and how to read the implementation critically.

Read Sections 1–5 first for the vacuum DIS picture, Sections 6–12 for the
medium, and Sections 13–17 for the connection to hadrons and measurements.
The final sections give a source-reading route and exercises with answers.
You can run every Python block independently with Python 3.10 or newer and
only the standard library. They calculate illustrative quantities without
generating events or reading simulation data. C++ blocks are excerpts from
this checkout, with surrounding code omitted; they are not complete programs.
All source paths below are relative to `src/ehijing/` unless explicitly stated.

This page explains existing prescriptions. Numerical examples are teaching
examples, not fits or recommended production settings. For the detailed
implementation reference, see the [physics manual](mainpage.md),
[assumptions review](assumptions.md), and [unit conventions](units.md).
The [eHIJING paper](https://doi.org/10.1103/PhysRevD.110.034001)
is the primary reference for the published model; this checkout also contains
ELECTRA-specific extensions described below.

## 1. What question is the generator answering?

Imagine sending an energetic electron into a nucleus. The electron exchanges
a virtual photon with a quark. The struck quark radiates and interacts with
the rest of the nucleus, then becomes part of a collection of hadrons. We
observe the scattered electron and some of those hadrons. Can we infer the
properties of the nuclear matter traversed by the quark?

An event generator makes a statistical prediction for that measurement. One
generated event is one sampled final state and its internal bookkeeping.
Histograms over many events approximate distributions; a single event does
not display an experimentally measurable microscopic history in full.

**Deep-inelastic scattering (DIS)** means that the lepton resolves the
target's partonic structure at a sufficiently short distance and produces a
hadronic system above the resonance region. **Semi-inclusive DIS (SIDIS)**
also identifies a hadron or jet in that final state. A jet is a reconstructed
collimated collection of particles, associated with an energetic quark or
gluon and its radiation; an observed jet is not itself a single quark.

The chain in this application is

\f[
 \begin{aligned}
 e+A&\longrightarrow e'+\text{hard parton + remnant}\\
    &\longrightarrow e'+\text{shower + nuclear interactions}\\
    &\longrightarrow e'+\text{color-connected partons}\\
    &\longrightarrow e'+\text{hadrons}.
 \end{aligned}
\f]

PYTHIA supplies the hard process, vacuum shower machinery, and Lund string
fragmentation. eHIJING adds a model of nuclear parton scattering and modified
radiation. In ELECTRA, the resulting hadrons may subsequently enter SMASH
for hadronic transport. This last stage is a separate model with its own
interaction and formation prescriptions.

The target medium here is cold nuclear matter: a sampled, frozen nuclear
configuration. It has no hydrodynamic temperature or evolving quark–gluon
plasma field. A final-state interaction can be partonic before fragmentation
or hadronic in downstream transport; these describe different stages.

## 2. Start with the measured electron: DIS kinematics

### 2.1 Four-momenta and the virtual photon

Let \f$k^\mu\f$ be the incoming electron momentum, \f$k'^\mu\f$ the
scattered electron momentum, and \f$P^\mu\f$ the target momentum **per
nucleon** used by the hard generator. Use the metric
\f$p^2=(p^0)^2-\mathbf p^2\f$. The momentum delivered to the hadronic
system is

\f[
 q^\mu=k^\mu-k'^\mu,\qquad Q^2=-q^2>0.
\f]

The photon is spacelike: its four-momentum squared is negative. It is an
internal exchange, so it need not obey the on-shell condition of a real
photon. \f$Q=\sqrt{Q^2}\f$ is called the hard scale or resolution scale.
The approximate transverse resolution is \f$1/Q\f$ in natural units.

Three useful invariants are

\f[
 x_B=\frac{Q^2}{2P\cdot q},\qquad
 y=\frac{P\cdot q}{P\cdot k},\qquad
 W^2=(P+q)^2.
\f]

At leading order in the collinear parton picture, \f$x_B\f$ is the momentum
fraction of the struck parton. Beyond that approximation it is a measured
kinematic variable, not an exact event-by-event parton label.
\f$y\f$ is the inelasticity: in the target rest frame it measures the
fraction of electron energy transferred. \f$W\f$ is the invariant mass
of the hadronic system associated with the hard per-nucleon reaction. It
is not the invariant mass of all \f$A\f$ explicit target nucleons.

### 2.2 Derive the target-rest-frame formulas

For \f$P^\mu=(m_N,\mathbf0)\f$, define the energy transfer
\f$\nu=q^0=E_e-E'_e\f$. Then \f$P\cdot q=m_N\nu\f$, so

\f[
 \nu=\frac{Q^2}{2m_Nx_B},\qquad y=\frac{\nu}{E_e}.
\f]

Expanding \f$(P+q)^2\f$ gives

\f[
 W^2=P^2+2P\cdot q+q^2
     =m_N^2+2m_N\nu-Q^2
     =m_N^2+Q^2\left(\frac1{x_B}-1\right).
\f]

Neglecting the electron mass, the electron scattering angle \f$\theta\f$
also gives \f$Q^2=4E_eE'_e\sin^2(\theta/2)\f$. For a fixed incoming
energy, \f$x_B,Q^2,y\f$ are therefore related. They are not three
independently adjustable variables.

As an illustration, choose \f$E_e=27.6\,\mathrm{GeV}\f$,
\f$Q^2=4\,\mathrm{GeV}^2\f$, \f$x_B=0.2\f$ and
\f$m_N=0.938\,\mathrm{GeV}\f$:

```python
from math import asin, degrees, sqrt

Ee, Q2, xB, mN = 27.6, 4.0, 0.2, 0.938
nu = Q2 / (2.0 * mN * xB)
y = nu / Ee
Ee_out = Ee - nu
W2 = mN**2 + Q2 * (1.0 / xB - 1.0)
sin2_half_theta = Q2 / (4.0 * Ee * Ee_out)
assert 0.0 < y < 1.0 and 0.0 <= sin2_half_theta <= 1.0
theta = 2.0 * asin(sqrt(sin2_half_theta))
print(f"nu = {nu:.3f} GeV, y = {y:.3f}")
print(f"W = {sqrt(W2):.3f} GeV, theta = {degrees(theta):.3f} degrees")
```

This gives \f$\nu\simeq10.661\,\mathrm{GeV}\f$,
\f$y\simeq0.386\f$, \f$W\simeq4.109\,\mathrm{GeV}\f$, and
\f$\theta\simeq5.30^\circ\f$. A large \f$Q^2\f$ alone does not
guarantee a large \f$W\f$; inspect both when selecting DIS.

The implementation uses Lorentz scalar products, which also work when the
nucleus moves. In `src/dis_kinematics.cpp`, `compute_dis_kinematics()` contains:

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

The `Vec4` multiplication is the Minkowski scalar product. The function
assumes the application's DIS event-record indexing. Reusing it for a
different PYTHIA process requires checking those particle identities first.
Configured `DISCuts` filter these reconstructed variables after
`pythia.next()` has completed its generation step.

### 2.3 Keep three different transverse momenta separate

We will write \f$\boldsymbol\kappa\f$ for a momentum transfer from the
medium, \f$\boldsymbol\ell\f$ for a radiative splitting momentum, and
\f$\mathbf p_{T,h}\f$ for a measured hadron momentum perpendicular to
the virtual-photon axis. Their transverse planes need not be identical:
internal kicks are defined relative to the local parton direction.

Similarly, the following fractions are different:

| Symbol | Meaning |
| --- | --- |
| \f$x_B\f$ | Lepton-reconstructed Bjorken variable. |
| \f$x_g\f$ | Momentum fraction assigned to an exchanged medium gluon. |
| \f$z\f$ | Emitted-gluon energy fraction in the high-energy splitting approximation. |
| \f$z_h=P\cdot p_h/(P\cdot q)\f$ | Measured SIDIS hadron fraction; \f$E_h/\nu\f$ in the target rest frame. |
| \f$z_L\f$ | Light-cone fraction used by string fragmentation. |

## 3. Units: do the dimensional analysis before the physics

Setting \f$\hbar=c=1\f$ identifies inverse energy with time and length.
It does not make a number quoted in fm equal to the same number in
\f$\mathrm{GeV}^{-1}\f$. The conversion is

\f[
 \hbar c=0.19732698\,\mathrm{GeV\,fm},\qquad
 1\,\mathrm{fm}\simeq5.068\,\mathrm{GeV}^{-1}.
\f]

The medium code uses positions and times in \f$\mathrm{GeV}^{-1}\f$,
density in \f$\mathrm{GeV}^3\f$, and line-integrated density, called
thickness, in \f$\mathrm{GeV}^2\f$. Standard stored PYTHIA particle
vertices use mm and mm/c. OSCAR/SMASH uses fm and fm/c. Thus

\f[
 t_{\rm fm/c}=(\hbar c)t_{\rm GeV^{-1}},\qquad
 x_{\rm mm}=10^{-12}x_{\rm fm}.
\f]

For example, \f$0.17\,\mathrm{fm}^{-3}\f$ becomes
\f$0.17(0.19732698)^3\simeq0.001306\,\mathrm{GeV}^3\f$.
A broadening coefficient in \f$\mathrm{GeV}^3\f$ must instead be
**divided** by \f$\hbar c\f$ to obtain
\f$\mathrm{GeV}^2/\mathrm{fm}\f$.

For a spacetime `Pythia8::Vec4`, construction is `(x, y, z, ct)` and `.e()`
accesses the time coordinate. For a momentum `Vec4`, `.e()` accesses energy.
The type alone cannot tell you which units or physical meaning it contains.
Convert all four components of a vertex together. A conversion of units
does not change a reference frame or add a nuclear origin.
The [unit conventions](units.md) also document the separate `PartonVertex`
subsystem's fm exception; the supplied shower settings disable it.

## 4. How the hard collision is sampled

### 4.1 PDFs are distributions of incoming partons

A parton distribution function (PDF) \f$f_{i/N}(x,\mu_F^2)\f$ encodes the
longitudinal parton content resolved at factorization scale \f$\mu_F\f$.
In the leading-order sampling picture, \f$f_i(x)dx\f$ supplies a parton
number density in a small momentum-fraction interval. It is not normalized
like the probability distribution of a single constituent quark. PDFs evolve
with resolution; the momentum sum rule involves \f$x f_i\f$ summed over
all partons, including gluons.

In the photon-exchange, massless, leading-order limit, the per-nucleon DIS
cross section can be written

\f[
 \frac{d^2\sigma}{dx_B\,dQ^2}
 =\frac{2\pi\alpha_{\rm EM}^2}{x_BQ^4}
  [1+(1-y)^2]F_2(x_B,Q^2),\qquad
 F_2=x_B\sum_f e_f^2(f_f+f_{\bar f}).
\f]

Here \f$\alpha_{\rm EM}\f$ is the electromagnetic coupling and
\f$e_f\f$ is the quark charge in units of the positron charge. The squared
charge weights which flavors the photon sees. The \f$Q^{-4}\f$ dependence
reflects the squared photon propagator. This expression omits target/lepton
mass effects, longitudinal structure-function contributions and weak terms;
it introduces the factorized picture rather than reproducing PYTHIA's full
hard-process implementation. In factorization, unresolved target structure
is in PDFs and short-distance scattering is in a calculable hard coefficient.

The supplied HERMES settings enable neutral-current photon/Z exchange and
disable charged-current W exchange. A nuclear event has one sampled hard
per-nucleon reaction, followed by the nuclear final-state model. It does not
have \f$A\f$ separately generated hard DIS collisions.

### 4.2 Nuclear PDFs and the outgoing medium describe different physics

The per-nucleon nuclear mixture is conventionally written

\f[
 f_{i/A}=\frac ZA f_{i/p|A}+\frac{A-Z}{A}f_{i/n|A},\qquad
 f_{i/p|A}=R_i^A f_{i/p}.
\f]

\f$A\f$ is the nucleon number, \f$Z\f$ the proton number, and
\f$R_i^A\f$ a fitted nuclear modification of the bound-proton PDF.
Neutron PDFs introduce the corresponding isospin mixture. Nuclear PDFs can
encode shadowing, antishadowing, the EMC effect and high-\f$x\f$ Fermi-motion
effects at the distribution level. They do not provide sampled nucleon Fermi
momenta or the positions and recoils of individual nucleons.

The driver in `src/main.cpp` selects EPPS16 NLO (`PDF:nPDFSetA = 3`) for
\f$Z>2\f$, and isospin only (`0`) otherwise. This overrides that setting
from the hard-process file. An NLO PDF input does not upgrade a leading-order
hard process to an NLO event generator.

The outgoing parton's scatterings use a separate effective gluon distribution
introduced in Section 7. Changing the hard nPDF and changing the rescattering
strength answer different questions; neither is a substitute for the other.

## 5. Why the outgoing quark produces a shower

### 5.1 Real branchings and their enhanced regions

A sufficiently energetic colored parton can radiate. In a collinear
approximation, quark-to-quark-plus-gluon radiation has real-emission density

\f[
 d\mathcal P_{q\to qg}
 =\frac{\alpha_s(\ell^2)}{2\pi}P_{gq}(z)
   dz\,\frac{d\ell^2}{\ell^2},\qquad
 P_{gq}(z)=C_F\frac{1+(1-z)^2}{z},\quad C_F=\frac43.
\f]

The emitted gluon carries fraction \f$z\f$. The \f$1/z\f$ factor enhances
soft radiation; \f$d\ell^2/\ell^2\f$ enhances collinear radiation.
\f$\alpha_s\f$ is the strong coupling. Integrating over a broad range
of scales creates large logarithms, so repeated branchings can be important
even when a single hard emission has a small probability.

PYTHIA approximates this cascade with a shower ordered in a transverse
evolution scale, including kinematic constraints, recoil, and competing
branchings. The displayed density is a useful limit, not a complete shower
algorithm or a normalized probability over all \f$z,\ell^2\f$.

### 5.2 Why one needs a probability of no branching

Suppose a resolved splitting rate is \f$\mathcal R(z,u)\f$, where
\f$u\f$ is a squared evolution scale. Let \f$u_1>u_2\f$. The probability
of no resolved splitting in that interval is the Sudakov factor

\f[
 \Delta(u_2,u_1)=\exp\left[
  -\int_{u_2}^{u_1}du\int_{z_{\min}(u)}^{z_{\max}(u)}dz\,
       \mathcal R(z,u)\right].
\f]

To understand the exponential, divide the interval into small pieces with
branching probabilities \f$d\mathcal P_i\f$. No branching requires
\f$\prod_i(1-d\mathcal P_i)\f$. Taking its logarithm and using
\f$\ln(1-d\mathcal P_i)\simeq-d\mathcal P_i\f$ produces the integral.
This makes the alternatives “branch” and “do not branch” consistent.

A veto algorithm samples from an easier overestimate and accepts candidates
with a rate ratio. Rejected trials are part of sampling a probability law;
they are not additional emitted particles. In eHIJING the selected medium
correction modifies the splitting kernel, and therefore its corresponding
no-emission construction, rather than merely adding a fixed number of gluons.

Initial-state radiation belongs to the incoming shower; final-state radiation
belongs to the outgoing shower. Recoil prescriptions decide how other momenta
change to conserve momentum during a branching. The current HERMES bundle
uses `SpaceShower:dipoleRecoil = on` and `TimeShower:interleave = off`;
the modified code retains final-state showering when eHIJING is enabled.
These settings matter when comparing the vacuum baseline with the paper.

## 6. Where is the quark created, and what does it cross?

### 6.1 Centers, density and the primary nucleon

The ELECTRA executable defaults to `--nuclear-target on`. It samples exactly
\f$A\f$ nucleon identities and exactly \f$Z\f$ protons. It chooses the
primary species with probability \f$Z/A\f$, then an identity uniformly
within that species. Its center becomes the hard vertex \f$\mathbf r_0\f$.
The primary identity is excluded from subsequent medium scattering.

For heavy targets, centers follow a spherical Woods–Saxon distribution:

\f[
 \begin{gathered}
 \rho_{\rm WS}(r)\propto\frac{1}{1+\exp[(r-R)/a]},\\
 p(r)dr\propto\frac{r^2dr}{1+\exp[(r-R)/a]},\qquad 0\leq r\leq R+10a.
 \end{gathered}
\f]

\f$R\f$ is a nuclear radius and \f$a\f$ the surface diffuseness.
The \f$r^2\f$ factor comes from the spherical volume element
\f$d^3r=r^2dr\,d\Omega\f$. Sampling just the Woods–Saxon profile
as a one-dimensional radius distribution would overweight the center.
This implementation samples independent centers without recentering,
deformation, or a minimum-distance correlation.

For deuterium, the proton–neutron separation follows a Hulthén radial law,
\f$p(r)\propto(e^{-\alpha_Hr}-e^{-\beta_Hr})^2\f$ with
\f$\alpha_H=0.228\,\mathrm{fm}^{-1}\f$ and
\f$\beta_H=1.18\,\mathrm{fm}^{-1}\f$. The centers lie at opposite
half-separations. For \f$A=1\f$ the only identity is the primary, so
there is no nonprimary rescattering medium.

Each remaining center \f$\mathbf r_j\f$ contributes a normalized Gaussian:

\f[
 \rho_j(\mathbf r)=\frac{e^{-|\mathbf r-\mathbf r_j|^2/(2w^2)}}
 {(2\pi w^2)^{3/2}},\qquad
 \rho_N(\mathbf r)=\sum_{j\ne j_{\rm primary}}\rho_j(\mathbf r).
\f]

The integral of \f$\rho_N\f$ over space is \f$A-1\f$. The default
\f$w=0.5\,\mathrm{fm}\f$ is a benchmark, not a fit established for
this implementation. Woods–Saxon specifies the distribution of centers;
the Gaussian specifies the density around each center. They play different
roles and introduce different lengths.

The primary species assignment is not explicitly conditioned on the generated
hard flavor's proton/neutron PDF contributions. Agreement of the average
proton fraction with \f$Z/A\f$ does not establish correct conditional
flavor–species correlations; the assumptions review tracks this open question.

### 6.2 Derive the forward-ray thickness

In the eikonal approximation the fast parton samples a straight ray,
\f$\mathbf r(s)=\mathbf r_s+s\mathbf n\f$ with \f$s\geq0\f$ and
\f$|\mathbf n|=1\f$. For one nucleon define the distance along the ray
\f$d_j=(\mathbf r_j-\mathbf r_s)\cdot\mathbf n\f$ and the squared
distance perpendicular to it
\f$b_j^2=|\mathbf r_j-\mathbf r_s|^2-d_j^2\f$. Then

\f[
 |\mathbf r(s)-\mathbf r_j|^2=b_j^2+(s-d_j)^2.
\f]

Substitution into the Gaussian separates a transverse factor and a
one-dimensional integral. Using
\f$u=(s-d_j)/(\sqrt2w)\f$ gives

\f[
 T_j=\int_0^\infty ds\,\rho_j(\mathbf r_s+s\mathbf n)
 =\frac{e^{-b_j^2/(2w^2)}}{4\pi w^2}
    \operatorname{erfc}\left(-\frac{d_j}{\sqrt2w}\right),\qquad
 T_A=\sum_{j\ne j_{\rm primary}}T_j.
\f]

The complementary error function is
\f$\operatorname{erfc}(v)=2\int_v^\infty e^{-u^2}du/\sqrt\pi\f$.
It enforces the forward starting boundary. A nucleon far ahead contributes
almost its full line integral; one far behind contributes only its Gaussian
tail. A larger transverse separation suppresses either contribution.

The joint prescription samples an identity with probability \f$T_j/T_A\f$
and then a distance with density
\f$p(s|j)=\rho_j(\mathbf r_s+s\mathbf n)/T_j\f$.
For speed \f$|\mathbf v|\f$, the elapsed rest-frame time is
\f$t=s/|\mathbf v|\f$. This is evaluated by
`NuclearTarget::path_weights()` and `sample_collision()` in
`external/pythia8/include/eHIJING/NuclearTarget.h`.

Here is an independent calculation for three illustrative centers relative
to a ray starting at the origin along positive \f$z\f$. Using fm throughout
this block gives thickness in \f$\mathrm{fm}^{-2}\f$:

```python
from math import erfc, exp, pi, sqrt

w_fm = 0.5
centers_fm = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 1.0)]
weights = []  # All three are nonprimary centers in this example.
for x, y, z in centers_fm:
    b2, d = x*x + y*y, z
    Tj = exp(-b2 / (2.0 * w_fm**2)) * erfc(-d / (sqrt(2.0) * w_fm))
    weights.append(Tj / (4.0 * pi * w_fm**2))
TA_fm2 = sum(weights)
TA_GeV2 = TA_fm2 * 0.19732698**2
assert abs(sum(Tj / TA_fm2 for Tj in weights) - 1.0) < 1e-12
print("identity probabilities:", [round(Tj / TA_fm2, 3) for Tj in weights])
print(f"TA = {TA_fm2:.4f} fm^-2 = {TA_GeV2:.5f} GeV^2")
```

The probabilities are approximately `0.863, 0.020, 0.117`. The center behind
the origin contributes weakly, and the off-axis center is suppressed. These
are path weights, not a sample of a physical nucleus.

### 6.3 What the eikonal approximation leaves fixed

The sampled density centers remain frozen even when remnant momenta change.
The original rays are not resampled after kicks. This is appropriate to the
intended fast-parton approximation; it does not implement slow, curved
trajectories through an evolving nucleus. Gaussian density has tails, so
there is no unique sharp exit point or path length \f$L\f$.

The legacy `--nuclear-target off` mode instead uses a uniform hard sphere,
\f$R_A=1.2A^{1/3}\,\mathrm{fm}\f$, and continuous density.
Its thickness and initial excluded interval differ from explicit identity
exclusion. The [physics manual](mainpage.md#nuclear-geometry-and-the-path-through-matter)
explains those differences; switching modes changes a physics prescription.

## 7. From the gluon distribution to a collision rate

### 7.1 What a rate differential in transverse transfer means

For a fast parton in color representation \f$R\f$, the local model is

\f[
 \frac{d\Gamma_R}{d^2\boldsymbol\kappa}
 =\frac{C_R}{d_A}\rho_N\frac{\alpha_s\phi_g(x_g,\kappa^2;Q_s^2)}{\kappa^2},
 \qquad x_g=\frac{x_B\kappa^2}{Q^2}.
\f]

Multiplying by \f$ds\,d^2\boldsymbol\kappa\f$ gives the infinitesimal
scattering probability in the relativistic path-length description.
For QCD's three colors, \f$d_A=N_c^2-1=8\f$,
\f$C_R=C_F=4/3\f$ for a quark or antiquark, and \f$C_R=C_A=3\f$ for
a gluon. These Casimirs quantify color-charge strength. For an otherwise
identical path and kernel the gluon rate is \f$9/4\f$ of the quark rate.

The effective transverse-momentum-dependent (TMD) gluon distribution is
modeled through the product

\f[
 \alpha_s\phi_g(x_g,\kappa^2;Q_s^2)
 =K\frac{(1-x_g)^n x_g^\lambda}{\kappa^2+Q_s^2}.
\f]

\f$K\f$ normalizes the interaction, \f$n\f$ suppresses transfers near
\f$x_g=1\f$, and negative \f$\lambda\f$ increases the small-\f$x_g\f$
weight. \f$Q_s\f$ sets the transverse turnover scale. The HERMES file
sets \f$n=4,\lambda=-0.25\f$; the XML default for \f$\lambda\f$ is
\f$-0.2\f$. The ansatz does not supply independent dynamical target
small-\f$x\f$ evolution.

Check units: \f$[\rho_N]=\mathrm{GeV}^3\f$,
\f$[\alpha_s\phi_g]=\mathrm{GeV}^{-2}\f$, and
\f$[d\Gamma/d^2\kappa]=\mathrm{GeV}^{-1}\f$.
Multiplication by \f$ds\f$ and \f$d^2\kappa\f$ is dimensionless.
In this normalization \f$K\f$ is dimensionless; the CLI supplies its
numerical value in the source convention.

`EHIJING::alphas_PhiG()` in `external/pythia8/src/eHIJING.cc` is:

```cpp
double alphas_PhiG(double x, double kT2, double Qs2,
           double powerG, double lambdaG, double K){
    return K * std::pow(1. - x, powerG) * std::pow(x, lambdaG) / (kT2 + Qs2);
}
```

The function already returns \f$\alpha_s\phi_g\f$. Adding another
factor of \f$\alpha_s\f$ would change the scattering model. Its `kT2`
argument is a collision transfer squared, not the lepton virtuality or an
observed hadron's transverse momentum.

### 7.2 A finite collision count needs an infrared resolution

For an azimuthally symmetric integrand,
\f$d^2\boldsymbol\kappa=\kappa\,d\kappa\,d\varphi
\to\pi\,d\kappa^2\f$. Define the resolved mean number for a fixed path:

\f[
 \overline N_R=\frac{C_RT_A}{d_A}
 \int_{\kappa_{\min}^2}^{Q^2/x_B}\pi\,d\kappa^2\,
       \frac{\alpha_s\phi_g}{\kappa^2}.
\f]

The upper limit enforces \f$x_g\leq1\f$. The numerical lower limit is
\f$\kappa_{\min}^2=0.01\mu^2\f$, with
\f$\mu=0.25\,\mathrm{GeV}\f$, hence
\f$\kappa_{\min}=0.025\,\mathrm{GeV}\f$.
The screened denominator \f$\kappa^2+Q_s^2\f$ alone does not render the
collision count finite: the rate also has \f$1/\kappa^2\f$ and its
small-\f$x_g\f$ factor. “How many collisions?” depends on resolution.

An ideal independent-scattering process gives
\f$\Pr(N)=e^{-\overline N}\overline N^N/N!\f$. The implementation
uses decreasing-transfer proposals with vetoes and then attaches times and
identities. Transfer order is not chronological order. Moreover, later
rejection of a whole event by exact recoil conditions the accepted sample.
The ideal Poisson expression is therefore an explanation of the dilute
measure, not a guarantee for the accepted collision-count histogram.

## 8. Why the saturation scale is an implicit equation

The distribution depends on \f$Q_s^2\f$, and \f$Q_s^2\f$ is determined
from the same distribution and traversed thickness:

\f[
 Q_s^2=\frac{C_A}{d_A}T_A
 \int_{\kappa_{\min}^2}^{Q^2/x_B}\pi\,d\kappa^2\,
       \alpha_s\phi_g(x_B\kappa^2/Q^2,\kappa^2;Q_s^2).
\f]

This is self-consistency, not a formula obtained by inserting a universal
fixed saturation momentum. An increase in \f$T_A\f$ increases the right
side for fixed \f$Q_s\f$; the resulting scale also changes the denominator
inside the integral. You must find a root for the requested
\f$T_A,x_B,Q^2,K,n,\lambda\f$.

The dimensions are instructive: the momentum integral of
\f$\alpha_s\phi_g\f$ is dimensionless, so multiplication by \f$T_A\f$
gives \f$\mathrm{GeV}^2\f$, as required. A larger nucleus often offers
a larger thickness, but neither \f$Q_s^2\f$ nor its path fluctuations are
determined by \f$A\f$ alone.

The code solves a residual \f$F(Q_s^2)=\text{right side}-Q_s^2\f$ by
bracketing and bisection. To integrate efficiently it substitutes

\f[
 u=\ln(1+\kappa^2/Q_s^2),\qquad
 d\kappa^2=(Q_s^2+\kappa^2)du.
\f]

The Jacobian cancels the screened denominator. `MultipleCollision` caches
solutions in `Qs.dat` on thickness and logarithmic kinematic grids;
explicit-target queries outside the grid solve directly. The root can lie
below \f$\mu^2\f$. A bracketing failure warns and returns \f$\mu^2\f$;
the shower separation scale has its own floor, discussed in Section 12.
Do not identify that floor with every saturation root.

Table filenames do not encode \f$K,n,\lambda\f$. A parameter change
requires matching table contents and separate directories as appropriate;
reusing an old `Qs.dat` does not recompute its physics. This tutorial does
not create or modify any such tables.

## 9. Broadening and the transport coefficient

### 9.1 A random-walk derivation

After \f$N\f$ transverse kicks, an idealized parton has
\f$\Delta\mathbf p_\perp=\sum_i\boldsymbol\kappa_i\f$. Squaring gives

\f[
 \left\langle\left|\sum_i\boldsymbol\kappa_i\right|^2\right\rangle
 =\sum_i\langle\kappa_i^2\rangle
  +2\sum_{i<j}\langle\boldsymbol\kappa_i\cdot\boldsymbol\kappa_j\rangle.
\f]

For independent isotropic kicks the cross terms vanish. The accumulated
mean squared momentum is the sum of the individual second moments. This
motivates the local transport coefficient

\f[
 \hat q_R(s)=\int d^2\boldsymbol\kappa\,\kappa^2
               \frac{d\Gamma_R(s)}{d^2\boldsymbol\kappa}
 =\frac{C_R}{d_A}\rho_N(s)\int d^2\boldsymbol\kappa\,\alpha_s\phi_g.
\f]

It measures transverse momentum squared acquired per unit length in this
approximation. It has natural units \f$\mathrm{GeV}^3\f$, commonly
reported as \f$\mathrm{GeV}^2/\mathrm{fm}\f$. Unlike the collision
count, its explicit \f$\kappa^2\f$ weight cancels the rate's
\f$1/\kappa^2\f$ factor. For finite positive \f$Q_s^2\f$, the
zero-cutoff second moment converges if \f$\lambda>-1\f$.

For identical integration bounds and the same path-dependent \f$Q_s\f$
in the saturation equation, integrating along the ray yields

\f[
 \int ds\,\hat q_R(s)=\frac{C_R}{C_A}Q_s^2.
\f]

For a quark this coefficient is \f$4/9\f$. If one specifies a finite
effective length \f$L\f$, one can define a path average
\f$\overline{\hat q}_R=(C_R/C_A)Q_s^2/L\f$. A Gaussian target has no
unique edge, so the chosen definition of \f$L\f$ must accompany it.
Also inspect helper functions before reporting their output: the current
legacy `qhatA()` and `qhatF()` helpers have inconsistent density factors,
as documented in the [physics manual](mainpage.md#transport-coefficient-and-transverse-momentum-broadening).

### 9.2 Why a hadron's broadening is a different quantity

The outgoing shower radiates and fragments. An approximate collinear picture
may assign a hadron a fraction of its parent's transverse momentum, but
fragmentation adds its own transverse momentum and energy loss changes which
hadrons enter a \f$z_h\f$ bin. Vacuum recoil and primordial transverse
momentum also contribute. A measured
\f$\Delta\langle p_{T,h}^2\rangle\f$ is therefore a conditional hadron
observable, not automatically \f$(4/9)Q_s^2\f$ or
\f$z_h^2(4/9)Q_s^2\f$.

## 10. Formation time and quantum interference

### 10.1 Derive the splitting energy denominator

Consider a high-energy, nearly massless parent of energy \f$E\f$. Its
quark and gluon daughters carry longitudinal fractions \f$1-z\f$ and
\f$z\f$, and opposite small transverse momenta of magnitude \f$\ell\f$.
Using \f$\sqrt{p_z^2+\ell^2}\simeq p_z+\ell^2/(2p_z)\f$, the
two-particle energy excess over the collinear parent is

\f[
 \Delta E\simeq\frac{\ell^2}{2zE}
                   +\frac{\ell^2}{2(1-z)E}
             =\frac{\ell^2}{2z(1-z)E}.
\f]

Relative quantum phases grow as \f$\Delta E\,t\f$. A characteristic
coherence or formation time is therefore

\f[
 \tau_f\sim\frac1{\Delta E}
          =\frac{2z(1-z)E}{\ell^2}.
\f]

This estimate assumes relativistic collinear light partons. Larger energy
lengthens coherence; larger transverse momentum shortens it. It is a
coordinate-time scale in the chosen nuclear rest frame, not a proper
lifetime of an on-shell unstable particle.

### 10.2 The LPM factor is an interference weight

Radiation associated with different scattering amplitudes interferes. The
Landau–Pomeranchuk–Migdal (LPM) factor appearing here is

\f[
 \Phi(t/\tau_f)=1-\cos(t/\tau_f).
\f]

When \f$t\ll\tau_f\f$, the amplitudes have not developed a substantial
relative phase, and \f$\Phi\simeq t^2/(2\tau_f^2)\f$. Treating each
collision as an independent fully formed radiation source would overcount
this regime. When phases vary rapidly over the sampled ensemble, the cosine
averages toward zero and \f$\langle\Phi\rangle\f$ approaches one.
At a specific time the factor oscillates between zero and two; its
pointwise value does not monotonically approach one.

```python
from math import cos

hbarc = 0.19732698
E, z, ell2 = 10.0, 0.2, 0.5  # GeV, dimensionless, GeV^2
tau_GeV_inv = 2.0 * z * (1.0 - z) * E / ell2
tau_fm = tau_GeV_inv * hbarc
print(f"tau = {tau_GeV_inv:.3f} GeV^-1 = {tau_fm:.3f} fm/c")
for t_fm in (0.1, 1.0, 3.0):
    phase = (t_fm / hbarc) / tau_GeV_inv
    print(f"t = {t_fm:.1f} fm/c: 1-cos(t/tau) = {1.0-cos(phase):.5f}")
```

Here \f$\tau_f=6.4\,\mathrm{GeV}^{-1}\simeq1.263\,\mathrm{fm}/c\f$.
The factors at 0.1, 1 and 3 fm/c are about 0.0031, 0.2975 and 1.72.
A sampled formation variable helps order and weight emissions. It does
not by itself establish a complete physical shower worldline.

## 11. What HT and GHT keep in the radiation kernel

### 11.1 Higher twist: expand in a small transfer

“Twist” classifies power-suppressed operator contributions in a hard-scale
expansion. In the nuclear problem, multiple scattering can enhance such
contributions. The practical HT kernel here uses a collinear expansion in
the medium transfer relative to the radiation momentum.

Define the vacuum density
\f$\mathcal V_a(z,\ell^2)=\alpha_s(\ell^2)P_{ga}(z)/(2\pi\ell^2)\f$.
For a sampled parent collision history
\f$\{(t_i,\boldsymbol\kappa_i)\}_{i=1}^N\f$, the HT representation is

\f[
 \begin{gathered}
 \frac{dN_a^{\rm HT}}{dz\,d\ell^2}
 =\mathcal V_a(z,\ell^2)\left[1+\mathcal W_{\rm HT}\right],\\
 \mathcal W_{\rm HT}=\frac{C_A}{C_R}\sum_i
 \Theta(\ell^2-\kappa_i^2)\frac{2\kappa_i^2}{\ell^2}
 \left(1-\cos\frac{t_i}{\tau_{\rm HT}}\right),\\
 \tau_{\rm HT}=\frac{2z(1-z)E}{\ell^2}.
 \end{gathered}
\f]

\f$\Theta\f$ is a step function selecting transfers below the radiation
scale. The `1` is the vacuum contribution. Every medium term contains the
relative transfer strength \f$\kappa_i^2/\ell^2\f$, a color factor and
the formation interference weight. HT is selected with
`--medium-modification-mode 0`.

In `external/pythia8/src/eHIJING.cc`, `sample_z_stochastic()` evaluates the
HT trial weight with:

```cpp
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

In this function `kt2` is the radiative scale, `q2` a collision transfer,
and `wmax` an overestimate using \f$1-\cos\leq2\f$. `w` itself can
exceed one: it is a kernel weight, not a probability. The rate ratio used
by the veto sampler is the candidate acceptance probability.

Replacing the collision sum with a density integral introduces a radiative
transport moment whose transfer upper bound is \f$\ell^2\f$. Section 9's
collisional moment used \f$Q^2/x_B\f$. Thus a reported \f$\hat q\f$
needs its scale and integration convention; these two moments are not
interchangeable constants.

### 11.2 Generalized higher twist: retain transfer dependence

The GHT soft-gluon representation retains transverse transfer in both
the propagator and phase:

\f[
 \frac{dN_a^{\rm GHT}}{dz\,d\ell^2}
 =\mathcal V_a(z,\ell^2)\left[1+\frac{C_A}{C_R}\sum_i
 \frac{2\boldsymbol\kappa_i\cdot\boldsymbol\ell}
 {|\boldsymbol\ell-\boldsymbol\kappa_i|^2}
 \left(1-\cos\frac{|\boldsymbol\ell-\boldsymbol\kappa_i|^2t_i}
 {2z(1-z)E}\right)\right].
\f]

The phase now depends on the relative vector
\f$\boldsymbol\ell-\boldsymbol\kappa_i\f$. This retains information
lost by the HT transfer expansion, especially when the two momenta are
comparable. An individual angular interference contribution can be negative;
it is not an independently sampled positive emission probability.
The high-scale implementation uses angular integrations, `GHT.dat`, Bessel
functions and approximation branches. The equation is a collision-conditioned
representation, not a literal list of independently generated emissions.

The apparent denominator singularity is softened by interference. Set
\f$D=|\boldsymbol\ell-\boldsymbol\kappa_i|^2\f$ and
\f$a=t_i/[2z(1-z)E]\f$. As \f$D\to0\f$,

\f[
 \frac{1-\cos(aD)}{D}\simeq\frac{a^2D}{2}\longrightarrow0.
\f]

GHT is selected with `--medium-modification-mode 1`. The implemented
soft-gluon approximation is not the unrestricted finite-\f$z\f$
generalized twist-four calculation. Keeping some \f$1-z\f$ factors
or finite-\f$z\f$ numerators does not remove that limitation. HT and
GHT are alternative approximations; their corrections are not added together.
With no sampled collisions, their medium terms vanish.

## 12. How two radiation stages share the work

### 12.1 High-scale shower

The modified PYTHIA shower contains vacuum radiation plus the chosen medium
correction above a separation scale. Schematically,

\f[
 \mathcal R_a(z,u)=\mathcal V_a(z,u)
      +\mathcal M_a(z,u)\Theta(u-Q_{\rm sep}^2).
\f]

Vacuum radiation continues to PYTHIA's shower cutoff. In explicit-target
mode the stored separation scale is floored:

\f[
 Q_{\rm sep}^2=\max(Q_s^2,\mu^2),\qquad
 \mu^2=0.0625\,\mathrm{GeV}^2.
\f]

The event stores it as `SeparationScale()`, a **squared momentum**, despite
the accessor's name. Shower evolution can update it. The low-scale stage
consumes the value stored in the event, not an independent fixed scale
selected in this tutorial.

### 12.2 Low-scale medium radiation

`Modified_FF::sample_ff_partons()` processes final light partons with sampled
collisions after the high-scale shower. It samples additional medium gluons
in increasing formation variable, starting with

\f[
 \tau_{\min}=\frac1{\sqrt{Q_{\rm sep}^2}},\qquad
 \tau_{\max}=\frac{E}{\mu^2},\qquad
 k_T^2=\frac{2z(1-z)E}{\tau}.
\f]

These are the local adapter's bounds, not simply \f$1/Q\f$ and the
nuclear exit time. Radiation is still suppressed by its collision-conditioned
weights. The adapter skips heavy quarks and diquarks; it does not provide
a heavy-quark energy-loss prescription when PYTHIA produces heavy flavor.

It uses a 0.2 GeV parton-energy resolution, requires parent energy above
0.4 GeV for the loop, and chooses
\f$z_{\min}=\min(0.2\,\mathrm{GeV}/E,0.4)\f$,
\f$z_{\max}=1-z_{\min}\f$. Its strong coupling is held at
\f$\alpha_s(Q_{\rm sep}^2)\f$. The local medium coupling in
`EHIJING::alphas()` is a frozen one-loop form,

\f[
 \alpha_s(\mathcal Q^2)=\frac{2\pi}{b_0\max[1,\ln(\mathcal Q^2/\mu^2)]},
 \qquad b_0=\frac92\quad(n_f=3).
\f]

The frozen value is about 1.40, so low-scale perturbative accuracy should
not be inferred merely from using a running-coupling formula. PYTHIA has
its own shower-coupling settings; this function does not control them all.

For HT the matching variable `lt2` equals `kt2`. In the low-scale GHT
construction it is instead

\f[
 \ell_T^2=k_T^2+\kappa_i^2
                   +2\sqrt{k_T^2\kappa_i^2}\cos\Delta\varphi.
\f]

The adapter keeps only \f$\ell_T^2\leq Q_{\rm sep}^2\f$, while the
new gluon is constructed with transverse magnitude \f$k_T\f$ and energy
\f$zE\f$, requiring \f$k_T\leq zE\f$. The distinction is visible in
`src/modified_ff.cpp`:

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

The first veto avoids the chosen scale overlap; the second makes the
massless emitted-gluon longitudinal momentum real. Matching is an algorithmic
prescription with vetoes, not proof of subtraction of every higher-order
overlap contribution.

New gluons sample later collisions from their assigned formation-variable
displacement and undergo elastic recoil. They do not recursively undergo a
second complete low-scale radiation cascade in the same call. Radiative
ordering and chronological elastic updates remain separate operations.

## 13. Elastic recoil, color and what remains of the nucleus

### 13.1 A kick transfers energy as well as transverse momentum

Explicit-target mode gives each nonprimary nucleon one shared on-shell
remnant momentum, initially \f$K_j^\mu=(m_N,\mathbf0)\f$ in the
nuclear rest frame. It sorts sampled elastic hits by rest-frame time,
including hits from different partons on the same identity.

For a parton \f$p\f$ and remnant \f$K_j\f$, an elastic update solves

\f[
 p'=p+\delta,\qquad K'_j=K_j-\delta,\qquad
 p'^2=p^2,\qquad K_j'^2=K_j^2,
 \qquad\boldsymbol\delta_\perp=\boldsymbol\kappa.
\f]

The sampled transverse transfer alone does not determine a valid collision.
The longitudinal transfer and energy must also satisfy both mass shells.
`ehijing::recoil_two_body()` in `include/nuclear_recoil.hpp` solves these
constraints, selecting the physical forward solution connected to zero
transfer. Positive energy and invariant masses are checked.

For intuition, consider a massless incoming parton of energy \f$E\f$
along positive \f$z\f$, hitting a stationary nucleon of mass \f$M\f$.
Let the parton lose energy \f$\omega\f$ and longitudinal momentum
\f$u>0\f$, while acquiring transverse magnitude \f$\kappa\f$. Its
final momentum is \f$p'=(E-\omega,\boldsymbol\kappa,E-u)\f$ in
energy-first notation. Adding the two on-shell equations gives
\f$u=(1+M/E)\omega\f$. Substituting into the target equation gives

\f[
 \left[\left(1+\frac ME\right)^2-1\right]\omega^2
 -2M\omega+\kappa^2=0.
\f]

The small root behaves as \f$\omega\simeq\kappa^2/(2M)\f$ for a
sufficiently small transfer. Even a purely transverse *sampled* transfer
causes energy loss once exact recoil is imposed. A separate drag term would
need its own justification to avoid counting that transfer twice.
For repeated hits the remnant moves, so this simple stationary-target
quadratic must be replaced by the full solver.

If a transfer has no physical forward solution, the application rejects the
entire generation attempt and regenerates its hard event and nucleus.
Transfers are not clipped inside accepted events. Therefore the accepted
distribution is conditioned by recoil feasibility. Rejection rates matter
when judging physical predictions, not merely execution speed.

### 13.2 Why color connections and shared identities matter

A parton also exchanges color with the target. Color labels tell the string
hadronizer which endpoints and gluon kinks form a connected singlet system.
They are not spatial coordinates or measured particle identities. A medium
kick can alter which partons are connected even for similar outgoing momenta.

After all hits, each touched nonprimary identity is split once into a quark
and diquark with its fixed proton/neutron valence content. A diquark here
is an effective two-quark string endpoint. Repeated hits update the same
remnant, exchange existing color lines, and do not create one extra baryon
per collision. Medium gluons are inserted into a supporting remnant's string.

The target bookkeeping distinguishes the primary, touched nonprimary
identities and untouched spectators. In neutral-current DIS, signed final
hadronic baryon number \f$B_{\rm final}\f$ and charge \f$Q_{\rm final}\f$
imply required spectator counts

\f[
 N_p=Z-Q_{\rm final},\qquad
 N_n=A-B_{\rm final}-N_p.
\f]

The scattered DIS lepton is excluded from that budget; final decay leptons
are counted with the hadronic system as appropriate. The code counts signed
quantum numbers, including antiparticles. In explicit mode the required
counts must match actual untouched identities or the attempt is rejected.
This conserves the target's baryon/charge ledger; it is not a model of
nuclear binding, de-excitation or evaporation.

### 13.3 Local elastic conservation does not establish global closure

Exact elastic recoil preserves four-momentum for that two-body update, and
the later remnant split preserves its four-momentum. The low-scale radiation
step still subtracts the gluon's spatial momentum and restores the parent
on shell independently. In `src/modified_ff.cpp`:

```cpp
particle.p(particle.p() - kmu);
particle.e(std::sqrt(particle.pAbs2() + particle.m2()));
```

After the second line, the new parent energy is determined by its mass shell,
so the parent-plus-gluon energy need not equal the incoming parent energy.
For a massless collinear parent, any nonzero opening angle already makes the
sum of the two daughter energies exceed the magnitude of their total spatial
momentum. Thus an on-shell two-body split needs a recoil or virtuality
prescription to close exactly. This implementation does not claim exact
global event energy conservation. Checking baryon number and charge alone
cannot detect that difference.

## 14. How partons become hadrons

At low scales the perturbative parton picture must be replaced by a
nonperturbative model. In the Lund picture, a color field between endpoints
forms a string with approximately constant energy per length. Gluons create
kinks. Quark–antiquark production breaks the string into color-singlet hadrons.

For a light-flavor string, the symmetric longitudinal fragmentation law has
the form

\f[
 f(z_L)\propto\frac{(1-z_L)^{a_L}}{z_L}
   \exp\left[-\frac{b_Lm_{\perp h}^2}{z_L}\right],\qquad
 m_{\perp h}^2=m_h^2+p_{\perp h}^2.
\f]

\f$a_L\f$ shapes the longitudinal distribution and \f$b_L\f$ has
inverse momentum-squared units. Flavor production, transverse momentum,
baryon formation, massive endpoints, junctions and small-mass string
treatment require additional rules. The fragmentation settings and tune
therefore affect measured hadrons even if the partonic medium kernel is
held fixed. An in-medium-modified parton system passed to vacuum string
fragmentation is not a perturbative calculation of confinement in a nucleus.

`Hadronizer` imports final partons and their color connections into a
separate PYTHIA instance with `ProcessLevel:all = off`. It preserves imported
medium displacements and adds the hard nuclear origin once. The import uses
status 23 without the source parent links. Its subsequent ancestry traversal
therefore does not reconstruct the missing high-scale shower delay.
Lund breakup and formation vertices are still generated. The supplied
hadronization bundle disables hadron decays for downstream handling;
custom settings must be inspected separately.

## 15. Four different meanings of “time”

Confusion about times can change a transport interpretation by orders of
magnitude. Identify the object before using a number:

| Quantity | Meaning in this workflow |
| --- | --- |
| Sampled collision time \f$t_i\f$ | Rest-frame time assigned along a frozen eikonal ray. |
| Radiation formation variable \f$\tau\f$ | Coherence/order variable in the emission sampler, with assigned medium displacements. |
| Hadron production reference | Available string breakup/origin information, with an explicit availability status. |
| Stored hadron formation vertex | PYTHIA's selected fragmentation vertex prescription, passed to output. |

The supplied vertex settings use the late hadron prescription. For an
ordinary light string segment, neglecting endpoint corrections and smearing,
its local form can be written

\f[
 x_{\rm form}^\mu=x_{\rm origin}^\mu
       +\frac{x_i^\mu+x_{i+1}^\mu}{2}
       +\frac{p_h^\mu}{2\kappa_{\rm str}}.
\f]

The breakup displacements \f$x_i,x_{i+1}\f$ are in fm and the string
tension \f$\kappa_{\rm str}\f$ in GeV/fm. It is unrelated to a momentum
transfer \f$\boldsymbol\kappa\f$. Special fragmentation branches use
additional prescriptions, so recorded references should be used instead
of extending this example formula to every particle.

The output writer converts full vertices to fm, then extrapolates a free
trajectory to a common \f$t=0\f$ surface:

\f[
 \mathbf x(0)=\mathbf x_{\rm form}
             -\frac{\mathbf p_h}{E_h}t_{\rm form}.
\f]

This is an input representation, not a claim that the hadron physically
formed at \f$t=0\f$. For example, a hadron formed at \f$x=2\,\mathrm{fm}\f$
and \f$t=3\,\mathrm{fm}/c\f$ with \f$v_x=0.6\f$ is represented at
\f$x(0)=0.2\,\mathrm{fm}\f$. Propagating that representation forward
freely for 3 fm/c returns 2 fm.

`prod_time` and `form_time` are output-frame coordinate times, not proper
time intervals requiring another Lorentz factor. `prod_time_status = 0`
means a missing production reference with a zero placeholder; `1` means
an available reference, including genuine zero; `2` means not applicable
for pre-existing spectators. SMASH uses zero and warns for status `0`.
That fallback does not establish when a missing physical production event
occurred. Its configured formation profile is an additional downstream model.

The optional HepMC history contains diagnostic observations, not inferred
interaction vertices for every shower branching. Its snapshot order follows
execution, which mixes shower scale ordering and separate time-ordered
elastic processing. Do not interpret snapshot sequence as a universal clock.

## 16. How nuclear effects appear in measurements

### 16.1 Multiplicity ratios

In an identical bin and acceptance for target \f$A\f$ and a reference
\f$D\f$ (often deuterium), define the per-DIS hadron yield

\f[
 M_A^h(z_h)=\frac1{N_{{\rm DIS},A}}\frac{dN_{h,A}}{dz_h},\qquad
 R_A^h(z_h)=\frac{M_A^h(z_h)}{M_D^h(z_h)}.
\f]

\f$N_{\rm DIS}\f$ counts selected DIS events, while \f$N_h\f$ counts
hadrons and may include several from one event. Dividing by DIS count
removes the trivial difference in sample sizes. In equal-width bins the
bin width cancels in the ratio. Radiation can deplete leading high-\f$z_h\f$
hadrons and populate lower fractions, but flavor mixtures, acceptance,
hadronization, remnant contributions and hadronic transport can also alter
the ratio. Suppression is not by itself a unique determination of energy loss.

For example, 240 hadrons in 1000 nuclear DIS events and 400 in 2000 reference
events give a ratio of \f$0.24/0.20=1.2\f$. Comparing the raw counts would
lead to the wrong conclusion.

### 16.2 Transverse broadening and a simple analysis block

Define

\f[
 \Delta\langle p_{T,h}^2\rangle
 =\langle p_{T,h}^2\rangle_A-\langle p_{T,h}^2\rangle_D.
\f]

Both samples must use the same hadron species, DIS cuts, \f$z_h\f$ bin,
photon-axis definition and detector acceptance. If events have statistical
weights, use weighted event and hadron sums consistently. The following
block uses deliberately tiny **unweighted synthetic** samples, already
selected into one common \f$z_h\f$ bin. Each tuple is
`(number_of_DIS_events, list_of_hadron_pT_squared_in_GeV2)`:

```python
from statistics import mean

samples = {
    "A": (5, [0.10, 0.30, 0.50, 0.70]),
    "D": (10, [0.10, 0.20, 0.30, 0.40, 0.50]),
}
def multiplicity(sample):
    n_dis, pt2 = sample
    assert n_dis > 0
    return len(pt2) / n_dis

ratio = multiplicity(samples["A"]) / multiplicity(samples["D"])
broadening = mean(samples["A"][1]) - mean(samples["D"][1])
print(f"multiplicity ratio = {ratio:.3f}")
print(f"broadening = {broadening:.3f} GeV^2")
```

The results are 1.600 and \f$0.100\,\mathrm{GeV}^2\f$. These values
demonstrate definitions; they contain no prediction or uncertainty estimate.
For real samples, account for event-level correlations between hadrons when
estimating errors; independently resampling each hadron can misstate them.
Use matched revisions/settings and distinguish eHIJING-only predictions from
the changes introduced by SMASH.

## 17. When is this physical picture appropriate?

The intended localized hard scattering and long-lived fragmenting system
are motivated by scale hierarchies such as

\f[
 \frac1Q\ll L,\qquad
 \ell_c\sim\frac{\nu}{Q^2}=\frac1{2m_Nx_B}\ll L,\qquad
 \tau_h\sim\frac{z_h\nu}{\Lambda_{\rm QCD}^2}\gg L.
\f]

The first compares transverse resolution with the nucleus. The second is a
longitudinal coherence/localization estimate: decreasing \f$x_B\f$ lets
the hard fluctuation extend over several nucleons. The third estimates
whether hadronization occurs well outside the target. Here
\f$\Lambda_{\rm QCD}\f$ denotes a characteristic nonperturbative QCD
scale, and \f$L\f$ a specified characteristic nuclear length. The last
estimate is not the generated Lund formation time for an individual hadron.

```python
hbarc, mN = 0.19732698, 0.938
Q2, xB, zh, Lambda = 4.0, 0.2, 0.5, 0.25
nu = Q2 / (2.0 * mN * xB)
resolution_fm = hbarc / Q2**0.5
coherence_fm = hbarc * nu / Q2
hadron_time_fm = hbarc * zh * nu / Lambda**2
print(f"transverse resolution = {resolution_fm:.3f} fm")
print(f"coherence estimate = {coherence_fm:.3f} fm")
print(f"hadron formation estimate = {hadron_time_fm:.3f} fm/c")
```

This gives roughly 0.099 fm, 0.526 fm and 16.83 fm/c. Against an illustrative
5 fm path the scales suggest short-distance resolution and substantial
hadronization delay, but “much less” or “much greater” must be assessed
over the actual event and hadron distribution. One example does not validate
an entire acceptance or prove the treatment of low-energy hadrons.

These conditions are scope estimates, not implicit generator cuts. The
code does not switch to a coherent photon–dipole hard-scattering model at
small \f$x_B\f$, or to a dense-medium radiation resummation, or to a
heavy-quark medium kernel. Fitted nPDF shadowing does not supply that
missing hard-event geometry. The published approximations and the added
explicit-target/recoil implementation require separate validation.

The [assumptions review](assumptions.md) records open questions about sampling,
matching, flavor conditioning, energy closure and spacetime. Synthetic
component tests check useful identities and interfaces; they do not
establish agreement with measured nuclear observables. The explicit-target
extensions have not yet been validated with a physics production sample.

## 18. A guided route through the implementation

Read one event in this order. The paths identify the operative code rather
than promising each intermediate record is a complete physical time slice.

| Step | Source to read | Question to answer |
| --- | --- | --- |
| Configuration and attempt loop | `src/main.cpp`, `src/cli.cpp` | Which settings are overridden, and which failures discard an attempt? |
| Lepton reconstruction | `src/dis_kinematics.cpp` | Which four-vectors define \f$x_B,Q^2,y,W^2\f$? |
| Target geometry | `external/pythia8/include/eHIJING/NuclearTarget.h` | Which identity is primary, and which density contributes along a ray? |
| Transfer and saturation sampling | `external/pythia8/src/eHIJING.cc` | What are the bounds, table inputs, and proposal/veto variables? |
| High-scale medium shower | `external/pythia8/src/TimeShower.cc` | Where does the medium term enter, and how is `SeparationScale()` stored? |
| Additional radiation | `src/modified_ff.cpp` | What do `kt2`, `lt2`, and `formation_time` mean at this point? |
| Elastic collisions | `include/nuclear_recoil.hpp`, `src/modified_ff.cpp` | Does the update conserve both momenta and masses, and which remnant is reused? |
| Strings and vertices | `src/hadronizer.cpp`, `external/pythia8/src/StringFragmentation.cc` | Which color/vertex information survives import, and which histories do not? |
| Output and diagnostics | `src/event_output.cpp`, `src/spacetime_history.cpp` | Which frame, units, time surface, availability flags and spectator ledger are written? |

For each variable, write down its physical quantity, units, reference frame
and coordinate origin. For each random draw, identify its proposal density
and acceptance condition. For each conservation claim, identify the exact
subset of particles and stage to which it applies. These habits make it
possible to distinguish an approximation from a coding error.

To inspect settings without generating events, run this from the ELECTRA
repository root:

```bash
rg -n 'WeakBosonExchange|SpaceShower|TimeShower|PDF:|eHIJING:' \
  input/ehijing/experiments/hermes_2007/hard_process.setting
rg -n 'alphas_PhiG|Qs2_self_consistent_eq|sample_z_stochastic' \
  src/ehijing/external/pythia8/src/eHIJING.cc
rg -n 'SeparationScale|formation_time|lt2|recoil_two_body' \
  src/ehijing/src/modified_ff.cpp
```

For production interpretation record the component/library revisions,
hard and hadronization settings, DIS cuts, \f$A,Z,K,n,\lambda\f$, target
mode and geometry parameters, table provenance, and random-stream policy.
Use the existing run guides for invocation. Choosing a different geometry,
kernel or formation prescription is a physics comparison to specify before
changing the model, not an exercise silently performed by this page.

## 19. Exercises with worked answers

### Exercise 1: Can these kinematics occur?

At \f$E_e=10\,\mathrm{GeV}\f$, request
\f$Q^2=4\,\mathrm{GeV}^2\f$, \f$x_B=0.1\f$ and
\f$m_N=0.938\,\mathrm{GeV}\f$. Is the point allowed in the
massless-electron target-rest-frame picture?

**Answer.** \f$\nu=4/(2\times0.938\times0.1)=21.32\,\mathrm{GeV}\f$,
so \f$y=2.13>1\f$ and \f$E'_e<0\f$. This is impossible. A plausible
\f$Q^2\f$ and \f$x_B\f$ considered separately need not form a
physical pair at a given beam energy.

### Exercise 2: Why does a Woods–Saxon radius need a Jacobian?

Near the center, take the spatial density as approximately constant. Compare
the probabilities in two equally narrow shells at radii \f$r\f$ and
\f$2r\f$.

**Answer.** Their shell volumes are proportional to \f$r^2dr\f$ and
\f$4r^2dr\f$, so the second contains approximately four times the
probability. Sampling a radius uniformly would miss that volume factor.

### Exercise 3: What broadening follows from a specified saturation scale?

For a teaching example, suppose the consistently defined path scale is
\f$Q_s^2=0.36\,\mathrm{GeV}^2\f$. What is the integrated quark
transport moment? What average follows if an effective length of 4 fm is
explicitly chosen?

**Answer.** The moment is
\f$(C_F/C_A)Q_s^2=(4/9)\times0.36=0.16\,\mathrm{GeV}^2\f$.
The specified average is \f$0.04\,\mathrm{GeV}^2/\mathrm{fm}\f$.
Neither number alone is the broadening of selected final hadrons. The
result requires the same transfer bounds as the saturation equation.

### Exercise 4: How does energy change formation interference?

Keep \f$z,\ell^2\f$ and a short collision time \f$t\ll\tau_f\f$
fixed, and double the parent energy. What happens to the LPM factor?

**Answer.** \f$\tau_f\propto E\f$ doubles, the phase halves, and
\f$1-\cos(t/\tau_f)\simeq t^2/(2\tau_f^2)\f$ falls by a factor
of four. This statement applies in the small-phase limit and does not
predict the full energy dependence after all phase-space integrations.

### Exercise 5: Do several hits imply several target baryons?

Three partons hit the same nonprimary nucleon. How many remnant quark–diquark
pairs does explicit-target mode create for that identity?

**Answer.** One. All hits update the same remnant momentum in chronological
order, and it is split once. Counting one pair per hit would add spurious
target baryon number and change the string system.

### Exercise 6: Diagnose two tempting conclusions

An output event has the expected final baryon number and charge. Its hadrons
are listed at \f$t=0\f$. Does this establish exact energy conservation
and instantaneous physical hadron formation?

**Answer.** Neither follows. Baryon/charge constraints do not test the
low-scale radiation energy residual. The \f$t=0\f$ positions are free
extrapolations from stored formation vertices; the separate production and
formation fields retain coordinate-time information and availability status.

## 20. Further reading

Use the references for the published theory and the source for this checkout's
actual algorithm. The tutorial's derivations and illustrative examples are
not new model settings or new validation results.

- [eHIJING: an Event Generator for Jet Tomography in Electron-Ion Collisions](https://doi.org/10.1103/PhysRevD.110.034001),
  [arXiv version](https://arxiv.org/abs/2304.10779): published model and phenomenology.
- [The local physics manual](mainpage.md#physics-included-in-ehijing): detailed
  kernels, source excerpts and implementation conventions.
- [Physical assumptions and open questions](assumptions.md): the current
  model/implementation distinction, unresolved issues and validation scope.
- [Unit conventions](units.md): the authoritative coordinate, unit and
  interface rules for this application.
- The bundled PYTHIA documentation under
  `external/pythia8/share/Pythia8/xmldoc/`, particularly
  `SpacelikeShowers.xml`, `TimelikeShowers.xml`, `Fragmentation.xml` and
  `VertexInformation.xml`: shower, string and vertex settings for the
  bundled version.
