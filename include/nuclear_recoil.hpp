#pragma once
#include "Pythia8/Basics.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ehijing {
// Exact elastic recoil for a prescribed transverse transfer in the incoming
// jet's axis frame. Choose the forward solution continuously connected to q=0.
// No clipping or resampling of the prescribed transverse momentum is allowed.
inline bool recoil_two_body(Pythia8::Vec4& jet, Pythia8::Vec4& target,
                            double qt2, double phi) {
    if (!(qt2 >= 0.) || !std::isfinite(qt2) || !std::isfinite(phi)
        || !(jet.e() > 0.) || !(target.e() > 0.) || !(jet.pAbs() > 0.)) return false;
    const double theta = jet.theta(), azimuth = jet.phi();
    auto p = jet, k = target;
    p.rot(0., -azimuth); p.rot(-theta, 0.);
    k.rot(0., -azimuth); k.rot(-theta, 0.);
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
    double best = std::numeric_limits<double>::infinity();
    Pythia8::Vec4 chosen;
    for (double z : roots) {
        const double energy = (b*z + c)/a;
        if (p.e()+energy <= 0. || k.e()-energy <= 0. || p.pz()+z <= 0.) continue;
        Pythia8::Vec4 q(qx, qy, z, energy);
        const auto pp = p+q, kk = k-q;
        const double jet_tolerance = 1e-10 * std::max({1., p.e()*p.e(), pp.e()*pp.e()});
        const double target_tolerance = 1e-10 * std::max({1., k.e()*k.e(), kk.e()*kk.e()});
        if (std::abs(pp.m2Calc()-p.m2Calc()) > jet_tolerance
            || std::abs(kk.m2Calc()-k.m2Calc()) > target_tolerance) continue;
        if (std::abs(z) < best) { best = std::abs(z); chosen = q; }
    }
    if (!std::isfinite(best)) return false;
    chosen.rot(theta, 0.); chosen.rot(0., azimuth);
    jet += chosen; target -= chosen;
    return true;
}
} // namespace ehijing
