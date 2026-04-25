#include "modified_ff.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace Pythia8;

Modified_FF::Modified_FF(int mode,
                         int atomic_number,
                         int mass_number,
                         double k_factor,
                         double xg_n,
                         double xg_lambda,
                         std::string table_path)
    : mode_(mode),
      atomic_number_(atomic_number),
      mass_number_(mass_number),
      z_over_a_(static_cast<double>(atomic_number) / mass_number),
      collision_sampler_(k_factor, xg_n, xg_lambda),
      ehijing_geometry_(mass_number, atomic_number),
      random_device_(),
      rng_(random_device_()),
      uniform_dist_(0.0, 1.0)
{
    collision_sampler_.Tabulate(table_path);
}

// The realization of the modified FF
void Modified_FF::sample_ff_partons(Pythia& pythia, const DISKinematics& kinematics, double& Rx, double& Ry, double& Rz)
{
    // Use fixed coupling at Qs of this event for anything medium-induced below Qs
    double kt2max_now = pythia.event.SeparationScale();
    double alpha_fix = EHIJING::alphas(kt2max_now);
    double alpha_bar = alpha_fix * EHIJING::CA / M_PI;
    double min_energy = .2;

    // This will hold the gluons emitted from the medium-contrbuted FF in this stage
    // as well as recoil partons that keeps the entire system color neutral.
    std::vector<Particle> new_particles;
    new_particles.clear();

    // For each parton from the high-Q shower, sample gluons and recoils from low-Q medium
    // corrections according to either HT or Generlized HT in the soft gluon limit
    std::vector<double> x_arr;
    std::vector<double> y_arr;
    std::vector<double> z_arr;
    std::vector<double> t_arr;

    for (int i = 0; i < pythia.event.size(); ++i)
    {
        
        auto& particle = pythia.event[i];

        // Only consider final state partons
        if (!particle.isFinal()) {
            continue;
        }

        // Only consider partons with energy above 2*min_energy in the nuclear rest frame
        // to avoid non-perturbative region and ensure the validity of the soft gluon approximation
        if (particle.e() < 2 * min_energy) {
            continue;
        }

        // Skip modification of heavy quarks
        /*
        Parton PIDs:
            d: 1    c: 4
            u: 2    b: 5
            s: 3    t: 6
            g: 21
        */
        auto abs_particle_id = particle.idAbs();
        bool is_light_parton = (abs_particle_id == 1) || (abs_particle_id == 2) || (abs_particle_id == 3) || (abs_particle_id == 21);
        if (!is_light_parton) {
            continue;
        }

        // Get its collision history
        std::vector<double> qt2s = particle.coll_qt2s(), ts = particle.coll_ts(), phis = particle.coll_phis();
        int n_collisions = ts.size();
        // If no collision, nothing to do
        if (n_collisions == 0) {
            continue;
        }

        // Compute particle velocity
        double vx = particle.px() / particle.e();
        double vy = particle.py() / particle.e();
        double vz = particle.pz() / particle.e();

        // Path length and nuclear thickness along the particle path
        double path_length = ehijing_geometry_.compute_L(pythia.event.Rx(), pythia.event.Ry(), pythia.event.Rz(), vx, vy, vz);
        double nuclear_thickness = ehijing_geometry_.compute_TA(pythia.event.Rx(), pythia.event.Ry(), pythia.event.Rz(), vx, vy, vz);

        // Location of the hard vertex
        Rx = pythia.event.Rx();
        Ry = pythia.event.Ry();
        Rz = pythia.event.Rz();

        // Useful quantity for the HT approach
        double sumq2 = 0.;
        for (int i = 0; i < n_collisions; ++i) {
            sumq2 += qt2s[i];
        }
        // Neglect too soft momentum kicks
        if (sumq2 < 1e-9) {
            continue;
        }

        // Formation time ordered fragmentation gluon
        // A very large cut off, since the LPM effect will effective regulate the formation time divergence
        double formation_time_max = particle.e() / EHIJING::mu2;
        // Start from the minimum formation time ~ 1/Qs, which is the separation scale;
        double formation_time_min = 1.0 / std::sqrt(pythia.event.SeparationScale());
        double formation_time = formation_time_min;

        double acceptance = 0;
        // z, kt2, lt2, phi_k, and dphiqk = phi_q - phik
        double z, kt2, lt2, phik, dphiqk;

        // Hold fragmented gluons and recoiled beam remnants
        std::vector<Particle> fragmented_gluons, recoiled_beam_remnants;
        fragmented_gluons.clear();
        recoiled_beam_remnants.clear();

        // Formation time loop
        while (formation_time < formation_time_max && particle.e() > 2 * min_energy) {
            double z_min = std::min(min_energy / particle.e(), .4);
            double z_max = 1. - z_min;
            if (z_max < z_min) {
                break;
            }
            double log_z_max = std::log(z_max / z_min);
            double maxdiffz = 1. / z_min - 1. / z_max + 2. * log_z_max;

            // Step 1: Find the next formation_time
            double r = uniform_dist_(rng_);

            // If mode = 0, use HT
            if (mode_ == 0) {
                double coeff = alpha_bar * maxdiffz * 4. * sumq2 / 2. / particle.e();
                formation_time = formation_time + std::log(1. / r) / coeff;
            // If mode = 1, use GHT
            } else if (mode_ == 1) {
                double invrpower = alpha_bar * log_z_max * 4. * n_collisions;
                double step_factor = std::pow(1. / r, 1. / invrpower);
                formation_time = formation_time * step_factor;
            }

            if (formation_time > formation_time_max || formation_time < formation_time_min) {
                break;
            }

            acceptance = 0.;
            // If mode = 0, use HT
            if (mode_ == 0) {
                for (int j = 0; j < n_collisions; ++j) {
                    double lpm_interference_phase_factor = (1. - std::cos(ts[j] / formation_time));
                    acceptance += qt2s[j] * lpm_interference_phase_factor;
                }
                acceptance /= (2. * sumq2);
            // If mode = 1, use GHT
            } else if (mode_ == 1) {
                for (int j = 0; j < n_collisions; ++j) {
                    double lpm_interference_phase_factor = (1. - std::cos(ts[j] / formation_time));
                    double z1mz = formation_time * qt2s[j] / 2. / particle.e();
                    if (z1mz > .25) {
                        acceptance += lpm_interference_phase_factor * log_z_max;
                    }
                    else {
                        double dz = std::sqrt(.25 - z1mz);
                        double z1 = .5 - dz;
                        double z2 = .5 + dz;
                        if (z1 > z_min) {
                            acceptance += lpm_interference_phase_factor * std::log(z1 / z_min);
                        }
                        if (z2 < z_max) {
                            acceptance += lpm_interference_phase_factor * std::log(z_max / z2);
                        }
                    }
                }
                acceptance /= (log_z_max * 2. * n_collisions);
            }


            if (acceptance < uniform_dist_(rng_)) {
                continue;
            }

            // Step 2: sample z, which also determines kt2
            acceptance = 0.;
            // If mode = 0, use HT
            if (mode_ == 0) {
                double N1 = 2 * (1. / z_min - 2.);
                double N2 = -4 * std::log(2. * (1. - z_max));
                double Ntot = N1 + N2;
                double r0 = N1 / Ntot;

                double acceptance = 0.;
                while (acceptance < uniform_dist_(rng_)) {
                    double r = uniform_dist_(rng_);
                    if (r < r0) {
                        z = z_min / (1. - z_min * r * Ntot / 2.);
                        acceptance = .5 / (1. - z);
                    } else {
                        z = 1. - std::exp(-(r * Ntot - N1) / 4.) / 2.;
                        acceptance = .25 / z / z;
                    }
                }

                // Inversion of the formation time formula to get k_T^2
                // Eq. (35) in https://doi.org/10.1103/PhysRevD.110.034001
                kt2 = 2 * (1. - z) * z * particle.e() / formation_time;

                // reject cases where qt2 > kt2 for mode = 0
                double Num = 0., Den = 0.;
                for (int j = 0; j < n_collisions; j++) {
                    if (ts[j] < 0) {
                        continue;
                    }

                    double q2 = qt2s[j], t = ts[j];
                    if (kt2 > q2) {
                        Num += q2 * (1. - std::cos(t / formation_time));
                    }

                    Den += q2 * (1. - std::cos(t / formation_time));
                }
                if (Num / Den < uniform_dist_(rng_)) {
                    continue;
                }
            // If mode = 1, use GHT
            } else if (mode_ == 1) {
                bool ok = false;
                double minimum_q2 = 2 * min_energy / formation_time;
                for (int j = 0; j < n_collisions; j++) {
                    if (qt2s[j] > minimum_q2 && ts[j] > 0) {
                        ok = true;
                    }
                }
                if (!ok) {
                    continue;
                }
                while (acceptance < uniform_dist_(rng_)) {
                    z = z_min * std::pow(z_max / z_min, uniform_dist_(rng_));

                    // Inversion of the formation time formula to get k_T^2
                    // Eq. (35) in https://doi.org/10.1103/PhysRevD.110.034001
                    kt2 = 2 * (1. - z) * z * particle.e() / formation_time;

                    double num = 0.;
                    for (int j = 0; j < n_collisions; j++) {
                        if (kt2 < qt2s[j] && ts[j] > 0) {
                            num += 1.;
                        }
                    }
                    acceptance = num / n_collisions;
                }
            }

            // Step 3: Correct for the splitting function
            // If particle is a gluon
            if (particle.id() == 21 && (1 + std::pow(1. - z, 3)) / 2. < uniform_dist_(rng_)) {
                continue;
            }
            // If particle is a quark
            if (particle.id() != 21 && (1 + std::pow(1. - z, 2)) / 2. < uniform_dist_(rng_)) {
                continue;
            }

            // Finally, sample phikT2 and compute the deflection of the hard parton lt2
            // If mode = 0, use HT
            if (mode_ == 0) {

                // No particular angular structure in the HT expansion
                phik = 2 * M_PI * uniform_dist_(rng_);
                lt2 = kt2;

            // If mode = 1, use GHT
            } else if (mode_ == 1) {
                
                double Psum = 0.;
                std::vector<double> dP;
                dP.resize(n_collisions);
                for (int j = 0; j < n_collisions; j++) {
                    if (kt2 < qt2s[j]) {
                        Psum += (1. - std::cos(ts[j] / formation_time));
                    }
                    dP[j] = Psum;
                }
                for (int j = 0; j < n_collisions; j++) {
                    dP[j] /= Psum;
                }
                double rc = uniform_dist_(rng_);
                int choice = -1;
                for (int j = 0; j < n_collisions; j++) {
                    if (rc < dP[j]) {
                        choice = j;
                        break;
                    }
                }
                // sample phik ~ (1+delta cos) / (1+delta^2 + 2 delta cos)
                double delta = std::sqrt(kt2 / qt2s[choice]);
                acceptance = 0.;
                while (acceptance < uniform_dist_(rng_)) {
                    r = uniform_dist_(rng_);
                    dphiqk = 2. * std::atan(std::tan(M_PI / 2. * r) * (delta + 1) / (delta - 1));
                    acceptance = (1 + delta * std::cos(dphiqk)) / 2.;
                }
                phik = phis[choice] + ((uniform_dist_(rng_) > .5) ? dphiqk : (-dphiqk));
                // lt2 = |(K-q) + q| = |k-q|^2
                lt2 = kt2 + qt2s[choice] + 2. * std::sqrt(kt2 * qt2s[choice]) * std::cos(dphiqk);
            }

            // The virtuality was not allowed to be overlap with the high-Q shower
            if (lt2 > kt2max_now) {
                continue;
            }

            // Now, there is a radiation,
            // Hard parton splits into p -> p-k & k
            double kt = std::sqrt(kt2), k0 = z * particle.e();
            if (kt > k0) {
                continue;
            }

            double kz = std::sqrt(k0 * k0 - kt2);
            Vec4 kmu{kt * std::cos(phik), kt * std::sin(phik), kz, k0};
            kmu.rot(particle.theta(), 0.);
            kmu.rot(0., particle.phi());
            particle.p(particle.p() - kmu);
            particle.e(std::sqrt(particle.pAbs2() + particle.m2()));
            kmu.e(std::sqrt(kmu.pAbs2()));

            // the gluon can continue to collide

            std::vector<double> g_qt2s, g_ts, g_phis;
            collision_sampler_.sample_all_qt2(21,
                                              kmu.e(),
                                              path_length,
                                              nuclear_thickness,
                                              kinematics.bjorken_x,
                                              kinematics.Q2,
                                              g_qt2s,
                                              g_ts,
                                              g_phis);
            Vec4 Qtot{0., 0., 0., 0.};
            double e0 = kmu.e();
            for (int ig = 0; ig < g_ts.size(); ig++) {
                double qt = std::sqrt(g_qt2s[ig]);
                double phi = g_phis[ig];
                Vec4 qmu{qt * std::cos(phi), qt * std::sin(phi), -qt * qt / 4. / e0, 0.0};
                Qtot = Qtot + qmu;
            }
            Qtot.rot(kmu.theta(), 0.);
            Qtot.rot(0., kmu.phi());
            // Elastic broadening
            kmu = kmu + Qtot;
            kmu.e(std::sqrt(kmu.pAbs2()));
            kmu = kmu * e0 / kmu.e();

            // update the color if it is a hard gluon
            // first, the spliting process
            int k_col, k_acol;
            // if the gluon forms inside the nuclei,
            // we consider it will lose color correlation with the original parton,
            // and form a new string with beam remnant
            Particle gluon = Particle(21, 201, i, 0, 0, 0, pythia.event.nextColTag(),
                                        pythia.event.nextColTag(), kmu, 0.0, 0);
            int qid, diqid;
            double mq, mdiq, mn;
            if (uniform_dist_(rng_) < z_over_a_) { // diquark from a proton
                diqid = 2101;
                mdiq = 0.57933;
                qid = 2;
                mq = 0.33;
                mn = 0.93847;
                if (uniform_dist_(rng_) < 2. / 3.) { // take away a u
                    qid = 2;
                    if (uniform_dist_(rng_) < .75) {
                        diqid = 2101;
                    } else {
                        diqid = 2103;
                    }
                } else { // take away the d
                    qid = 1;
                    diqid = 2203;
                }
            } else { // diquark from a neutron
                diqid = 2101;
                mdiq = 0.57933;
                qid = 1;
                mq = 0.33;
                mn = 0.93957;
                if (uniform_dist_(rng_) < 2. / 3.) { // take away a d
                    qid = 1;
                    if (uniform_dist_(rng_) < .75) {
                        diqid = 2101;
                    } else {
                        diqid = 2103;
                    }
                } else {
                    // take away the u
                    qid = 2;
                    diqid = 1103;
                }
            }
            double pabs = std::sqrt((mn * mn - std::pow(mq + mdiq, 2)) *
                                    (mn * mn - std::pow(mq - mdiq, 2))) /
                            (2. * mn);
            double costheta = uniform_dist_(rng_) * 2. - 1.;
            double sintheta = std::sqrt(std::max(1. - costheta * costheta, 1e-9));
            double rphi = 2 * M_PI * uniform_dist_(rng_);
            double Nqz = pabs * costheta, Nqx = pabs * sintheta * std::cos(rphi),
                    Nqy = pabs * sintheta * std::sin(rphi);
            Vec4 pq{Nqx, Nqy, Nqz, 0}, pdiq{-Nqx, -Nqy, -Nqz, 0};
            pq.e(std::sqrt(pq.pAbs2() + mq * mq));
            pdiq.e(std::sqrt(pdiq.pAbs2() + mdiq * mdiq));
            Particle recolQ = Particle(qid, 201, i, 0, 0, 0, gluon.acol(), 0, pq, mq, 0);
            Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0, 0, gluon.col(), pdiq, mdiq, 0);
            // time formation_time
            recoiled_beam_remnants.push_back(gluon);
            x_arr.push_back(formation_time * HBARC * kmu.px() / kmu.e());
            y_arr.push_back(formation_time * HBARC * kmu.py() / kmu.e());
            z_arr.push_back(formation_time * HBARC * kmu.pz() / kmu.e());
            t_arr.push_back(formation_time * HBARC);

            recoiled_beam_remnants.push_back(recolQ);
            x_arr.push_back(formation_time * HBARC * pq.px() / pq.e());
            y_arr.push_back(formation_time * HBARC * pq.py() / pq.e());
            z_arr.push_back(formation_time * HBARC * pq.pz() / pq.e());
            t_arr.push_back(formation_time * HBARC);

            recoiled_beam_remnants.push_back(recoldiQ);
            x_arr.push_back(formation_time * HBARC * pdiq.px() / pdiq.e());
            y_arr.push_back(formation_time * HBARC * pdiq.py() / pdiq.e());
            z_arr.push_back(formation_time * HBARC * pdiq.pz() / pdiq.e());
            t_arr.push_back(formation_time * HBARC);
        } // Comment out this loop to turn off radiation

        // Handle recoil and remnants
        // If there is radiation, recoil goes to radiation.
        // If not, goes to the hard quark
        int Nrad = fragmented_gluons.size();
        for (int j = 0; j < n_collisions; j++) {
            double qT = std::sqrt(qt2s[j]), phiq = phis[j];
            double qx = qT * std::cos(phiq);
            double qy = qT * std::sin(phiq);
            double qz = -qT * qT / 4. / particle.e();
            Vec4 qmu{qx, qy, qz, 0};
            qmu.rot(particle.theta(), 0.);
            qmu.rot(0., particle.phi());
            // turn off the elastic broading Wenbin
            particle.p(particle.p() + qmu);
            particle.e(std::sqrt(particle.pAbs2() + particle.m2()));
            int q_col, q_acol;
            // update color
            if (particle.id() == 21) {
                if (std::rand() % 2 == 0) {
                    q_acol = particle.acol();
                    particle.acol(pythia.event.nextColTag());
                    q_col = particle.acol();
                } else {
                    q_col = particle.col();
                    particle.col(pythia.event.nextColTag());
                    q_acol = particle.col();
                }
            } else if (particle.id() > 0) {
                q_col = particle.col();
                particle.col(pythia.event.nextColTag());
                q_acol = particle.col();
            } else {
                q_acol = particle.acol();
                particle.acol(pythia.event.nextColTag());
                q_col = particle.acol();
            }
            int qid, diqid;
            double mq, mdiq, mn;
            if (uniform_dist_(rng_) < z_over_a_) { // diquark from a proton
                diqid = 2101;
                mdiq = 0.57933;
                qid = 2;
                mq = 0.33;
                mn = 0.93847;
                if (uniform_dist_(rng_) < 2. / 3.) { // take away a u
                    qid = 2;
                    if (uniform_dist_(rng_) < .75) {
                        diqid = 2101;
                    } else {
                        diqid = 2103;
                    }
                } else { // take away the d
                    qid = 1;
                    diqid = 2203;
                }
            } else { // diquark from a neutron
                diqid = 2101;
                mdiq = 0.57933;
                qid = 1;
                mq = 0.33;
                mn = 0.93957;
                if (uniform_dist_(rng_) < 2. / 3.) { // take away a d
                    qid = 1;
                    if (uniform_dist_(rng_) < .75) {
                        diqid = 2101;
                    } else {
                        diqid = 2103;
                    }
                } else { // take away the u
                    qid = 2;
                    diqid = 1103;
                }
            }
            double pabs =
                std::sqrt((mn * mn - std::pow(mq + mdiq, 2)) * (mn * mn - std::pow(mq - mdiq, 2))) /
                (2. * mn);
            double costheta = uniform_dist_(rng_) * 2. - 1.;
            double sintheta = std::sqrt(std::max(1. - costheta * costheta, 1e-9));
            double rphi = 2 * M_PI * uniform_dist_(rng_);
            double Nqz = pabs * costheta, Nqx = pabs * sintheta * std::cos(rphi),
                   Nqy = pabs * sintheta * std::sin(rphi);
            Vec4 pq{Nqx, Nqy, Nqz, 0}, pdiq{-Nqx, -Nqy, -Nqz, 0};
            // decide which object takes the recoil
            if (std::rand() % 2 == 0) {
                pq = pq - qmu;
            } else {
                pdiq = pdiq - qmu;
            }
            pq.e(std::sqrt(pq.pAbs2() + mq * mq));
            pdiq.e(std::sqrt(pdiq.pAbs2() + mdiq * mdiq));
            Particle recolQ = Particle(qid, 201, i, 0, 0, 0, q_col, 0, pq, mq, 0);
            Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0, 0, q_acol, pdiq, mdiq, 0);

            recoiled_beam_remnants.push_back(recolQ);
            x_arr.push_back(ts[j] * HBARC * pq.px() / pq.e());
            y_arr.push_back(ts[j] * HBARC * pq.py() / pq.e());
            z_arr.push_back(ts[j] * HBARC * pq.pz() / pq.e());
            t_arr.push_back(ts[j] * HBARC);

            recoiled_beam_remnants.push_back(recoldiQ);
            x_arr.push_back(ts[j] * HBARC * pdiq.px() / pdiq.e());
            y_arr.push_back(ts[j] * HBARC * pdiq.py() / pdiq.e());
            z_arr.push_back(ts[j] * HBARC * pdiq.pz() / pdiq.e());
            t_arr.push_back(ts[j] * HBARC);
        }

        for (auto& fragmented_gluon : fragmented_gluons) {
            new_particles.push_back(fragmented_gluon);
        }

        for (auto& recoiled_beam_remnant : recoiled_beam_remnants) {
            new_particles.push_back(recoiled_beam_remnant);
        }
    }

    int ixyz = 0;
    for (auto& new_particle : new_particles) {
        
        pythia.event.append(new_particle.id(),
                            201, 
                            new_particle.col(), 
                            new_particle.acol(), 
                            new_particle.px(), 
                            new_particle.py(), 
                            new_particle.pz(), 
                            new_particle.e(), 
                            new_particle.m());
        
        auto& appended = pythia.event.back();
        appended.vProd(Vec4(
            x_arr[ixyz] + Rx * HBARC,
            y_arr[ixyz] + Ry * HBARC,
            z_arr[ixyz] + Rz * HBARC,
            t_arr[ixyz]
        ));

        ixyz++;
    }
}
