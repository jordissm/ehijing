#ifndef EHIJING_NUCLEARTARGET_H
#define EHIJING_NUCLEARTARGET_H

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace EHIJING {

// All positions, widths and times use GeV^-1; path thickness uses GeV^2.
// Geometry follows TRENTo's Woods-Saxon/Hulthen prescriptions, not its
// heavy-ion entropy deposition or nucleon-nucleon interaction model.
struct TargetNucleon {
    bool proton = false;
    std::array<double, 3> position{};
    int hits = 0;
};

class NuclearTarget {
public:
    static constexpr double hbarc = 0.19732698;
    NuclearTarget(int A, int Z, double width_fm,
                  double radius_fm = 0., double surface_fm = 0.)
        : nucleons_(A > 0 ? A : 0), Z_(Z), width_(width_fm / hbarc) {
        if (A < 1 || Z < 0 || Z > A || !(width_fm > 0.) || !std::isfinite(width_fm))
            throw std::invalid_argument("Invalid explicit nuclear target");
        if (A > 2) {
            if (radius_fm == 0. && surface_fm == 0.) {
                switch (A) {
                case 63: radius_fm = 4.20; surface_fm = .596; break;
                case 129: radius_fm = 5.36; surface_fm = .590; break;
                case 197: radius_fm = 6.38; surface_fm = .535; break;
                case 208: radius_fm = 6.62; surface_fm = .546; break;
                case 238: radius_fm = 6.81; surface_fm = .600; break;
                default: throw std::invalid_argument(
                    "No TRENTo radius/surface for this A: supply --nuclear-radius-fm and --nuclear-surface-fm");
                }
            }
            if (!(radius_fm > 0.) || !(surface_fm > 0.)
                || !std::isfinite(radius_fm) || !std::isfinite(surface_fm))
                throw std::invalid_argument("Woods-Saxon radius and surface must be positive");
        }
        radius_ = radius_fm / hbarc;
        surface_ = surface_fm / hbarc;
        if (A == 2 && Z != 1)
            throw std::invalid_argument("Hulthen target requires A=2, Z=1");
    }

    void sample(const std::function<double()>& random) {
        for (size_t j = 0; j < nucleons_.size(); ++j) {
            nucleons_[j] = TargetNucleon{};
            nucleons_[j].proton = static_cast<int>(j) < Z_;
        }
        if (nucleons_.size() == 2) {
            // r is the proton-neutron separation, not either COM radius.
            const double a = .228 * hbarc, b = 1.18 * hbarc;
            double r;
            do { r = -std::log(open(random())) / (2. * a); }
            while (random() > std::pow(1. - std::exp(-(b - a) * r), 2));
            nucleons_[0].position = direction(.5 * r, random);
            for (int k = 0; k < 3; ++k) nucleons_[1].position[k] = -nucleons_[0].position[k];
        } else if (nucleons_.size() > 2) {
            // Volume proposal implements the r^2 Jacobian. No minimum-distance
            // correlations or recentering distortion are introduced.
            const double maximum = radius_ + 10. * surface_;
            for (auto& nucleon : nucleons_) {
                double r;
                do { r = maximum * std::cbrt(random()); }
                while (random() > 1. / (1. + std::exp((r - radius_) / surface_)));
                nucleon.position = direction(r, random);
            }
        }
        // Select species with Z/A, then uniformly select an identity of that species.
        const bool proton = random() < static_cast<double>(Z_) / nucleons_.size();
        const int first = proton ? 0 : Z_;
        const int count = proton ? Z_ : static_cast<int>(nucleons_.size()) - Z_;
        primary_ = first + std::min(count - 1, static_cast<int>(random() * count));
        nucleons_[primary_].hits = 1;
    }

    int primary() const { return primary_; }
    int Z() const { return Z_; }
    double width() const { return width_; }
    double radius() const { return radius_; }
    double surface() const { return surface_; }
    const std::vector<TargetNucleon>& nucleons() const { return nucleons_; }
    std::vector<TargetNucleon>& nucleons() { return nucleons_; }

    double density(const std::array<double, 3>& x) const {
        double sum = 0.;
        for (size_t j = 0; j < nucleons_.size(); ++j) {
            if (static_cast<int>(j) == primary_) continue;
            double distance2 = 0.;
            for (int k = 0; k < 3; ++k) distance2 += std::pow(x[k] - nucleons_[j].position[k], 2);
            sum += std::exp(-distance2 / (2. * width_ * width_));
        }
        return sum / std::pow(2. * pi * width_ * width_, 1.5);
    }

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

    double thickness(const std::array<double, 3>& start,
                     const std::array<double, 3>& velocity) const {
        double result = 0.;
        for (double w : path_weights(start, velocity)) result += w;
        return result;
    }

    // Joint sample of identity and distance along the positive ray. This
    // samples exactly rho_j(start + s*n) ds / integrated_total_thickness.
    int sample_collision(const std::array<double, 3>& start,
                         const std::array<double, 3>& velocity,
                         const std::function<double()>& random, double& time) const {
        const auto weights = path_weights(start, velocity);
        double total = 0.;
        for (double w : weights) total += w;
        if (!(total > 0.)) throw std::runtime_error("Cannot sample collision on zero thickness");
        double choice = random() * total;
        int j = -1;
        for (size_t i = 0; i < weights.size(); ++i) {
            if (weights[i] == 0.) continue;
            j = static_cast<int>(i);
            choice -= weights[i];
            if (choice <= 0.) break;
        }
        const double speed = norm(velocity);
        double centre = 0.;
        for (int k = 0; k < 3; ++k)
            centre += (nucleons_[j].position[k] - start[k]) * velocity[k] / speed;
        const double tail = open(random()) * std::erfc(-centre / (std::sqrt(2.) * width_));
        double low = 0., high = std::max(0., centre) + 16. * width_;
        for (int iteration = 0; iteration < 80; ++iteration) {
            const double mid = .5 * (low + high);
            if (std::erfc((mid - centre) / (std::sqrt(2.) * width_)) > tail) low = mid;
            else high = mid;
        }
        time = .5 * (low + high) / speed;
        return j;
    }

private:
    static constexpr double pi = 3.14159265358979323846;
    static double open(double u) { return std::max(1e-15, std::min(1. - 1e-15, u)); }
    static double norm(const std::array<double, 3>& v) {
        return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    static std::array<double, 3> direction(double r, const std::function<double()>& random) {
        const double cosine = 2. * random() - 1., phi = 2. * pi * random();
        const double sine = std::sqrt(std::max(0., 1. - cosine * cosine));
        return {r * sine * std::cos(phi), r * sine * std::sin(phi), r * cosine};
    }
    std::vector<TargetNucleon> nucleons_;
    int Z_, primary_ = -1;
    double width_, radius_, surface_;
};

// Link-time guard: rebuilding the modified PYTHIA library is required.
int nuclear_target_abi_version();
} // namespace EHIJING
#endif
