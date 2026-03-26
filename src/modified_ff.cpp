#include "modified_ff.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace Pythia8;

Modified_FF::Modified_FF(int mode_,
                         int Z_,
                         int A_,
                         double K,
                         double n,
                         double lambda,
                         std::string TablePath)
    : mode(mode_),
      Z(Z_),
      A(A_),
      ZoverA(Z * 1.0 / A),
      Coll(K, n, lambda),
      eHIJING_Geometry(A, Z),
      rd(),
      gen(rd()),
      dist(0.0, 1.0)
{
    Coll.Tabulate(TablePath);
}

// The realization of the modified FF
void Modified_FF::sample_FF_partons(Event& event, double& Rx, double& Ry, double& Rz)
{
    Vec4 pProton = event[1].p();   // four-momentum of proton
    Vec4 peIn = event[4].p();      // incoming electron
    Vec4 peOut = event[6].p();     // outgoing electron
    Vec4 pGamma = peIn - peOut;    // virtual boson photon/Z^0/W^+-
    double Q20 = -pGamma.m2Calc(); // hard scale square
    double Q0 = std::sqrt(Q20);
    double xB = Q20 / (2. * pProton * pGamma); // Bjorken x
    double nu = pGamma.e();
    double W2 = (pProton + pGamma).m2Calc();
    auto& hardP = event[5];

    // Use fixed coupling at Qs of this event for anything medium-induced below Qs
    double kt2max_now = event.SeparationScale();
    double alpha_fix = EHIJING::alphas(kt2max_now);
    double alphabar = alpha_fix * EHIJING::CA / M_PI;
    double emin = .2;

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
    for (int i = 0; i < event.size(); ++i)
    {

        // only find final-state quarks and gluons with energy above
        // 2*emin in the nuclear rest frame
        auto& p = event[i];
        auto abspid = p.idAbs();
        if (!p.isFinal())
            continue;
        if (p.e() < 2 * emin)
            continue;

        // we are skipping modifications of heavy quarks right now!!!
        bool isLightParton = (abspid == 1) || (abspid == 2) || (abspid == 3) || (abspid == 21);
        if (!isLightParton)
            continue;

        // Get its collision history
        std::vector<double> qt2s = p.coll_qt2s(), ts = p.coll_ts(), phis = p.coll_phis();
        // if no collision, nothing to do
        int Ncolls = ts.size();
        if (Ncolls == 0)
            continue;

        double vx = p.px() / p.e(), vy = p.py() / p.e(), vz = p.pz() / p.e();
        double L = eHIJING_Geometry.compute_L(event.Rx(), event.Ry(), event.Rz(), vx, vy, vz);
        double TA = eHIJING_Geometry.compute_TA(event.Rx(), event.Ry(), event.Rz(), vx, vy, vz);
        Rx = event.Rx();
        Ry = event.Ry();
        Rz = event.Rz();

        double sumq2 = 0.; // useful quantity for H-T approach
        for (int i=0; i < Ncolls; ++i) sumq2 += qt2s[i];
        if (sumq2 < 1e-9)
            continue; // negelect too soft momentum kicks

        // tauf ordered fragmentation gluon
        // A very large cut off, since the LPM effect will effective regulate the tauf divergence
        double taufmax = p.e() / EHIJING::mu2;
        // Start from the minimum tauf ~ 1/Qs, which is the separation scale;
        double taufmin = 1.0 / std::sqrt(event.SeparationScale());
        double tauf = taufmin;

        double acceptance = 0;
        // z, kt2, lt2, phi_k, and dphiqk = phi_q - phik
        double z, kt2, lt2, phik, dphiqk;

        // holds fragmented gluons and recoiled beam remnants
        std::vector<Particle> frag_gluons, recoil_remnants;
        frag_gluons.clear();
        recoil_remnants.clear();

        // Formation time loop
        // Commnent out this loop to turn off the radiation Wenbin // begin of the loop
        while (tauf < taufmax && p.e() > 2 * emin)
        {
            double zmin = std::min(emin / p.e(), .4);
            double zmax = 1. - zmin;
            if (zmax < zmin)
                break;
            double maxlogz = std::log(zmax / zmin);
            double maxdiffz = 1. / zmin - 1. / zmax + 2. * maxlogz;
            // step1: find the next tauf
            double r = dist(gen);
            if (mode == 1)
            {
                double invrpower = alphabar * maxlogz * 4. * Ncolls;
                double step_factor = std::pow(1. / r, 1. / invrpower);
                tauf = tauf * step_factor;
            } else
            {
                double coeff = alphabar * maxdiffz * 4. * sumq2 / 2. / p.e();
                tauf = tauf + std::log(1. / r) / coeff;
            }
            if (tauf > taufmax || tauf < taufmin)
                break;
            acceptance = 0.;
            if (mode == 1)
            {
                for (int j = 0; j < Ncolls; j++)
                {
                    double phase = (1. - std::cos(ts[j] / tauf));
                    double z1mz = tauf * qt2s[j] / 2. / p.e();
                    if (z1mz > .25)
                        acceptance += phase * maxlogz;
                    else
                    {
                        double dz = std::sqrt(.25 - z1mz);
                        double z1 = .5 - dz;
                        double z2 = .5 + dz;
                        if (z1 > zmin)
                            acceptance += phase * std::log(z1 / zmin);
                        if (z2 < zmax)
                            acceptance += phase * std::log(zmax / z2);
                    }
                }
                acceptance /= (maxlogz * 2. * Ncolls);
            } else
            {
                for (int j = 0; j < Ncolls; j++)
                    acceptance += qt2s[j] * (1. - std::cos(ts[j] / tauf));
                acceptance /= (2. * sumq2);
            }
            if (acceptance < dist(gen))
                continue;

            // step 2: sample z, which also determines kt2
            acceptance = 0.;
            if (mode == 0)
            {
                double N1 = 2 * (1. / zmin - 2.);
                double N2 = -4 * std::log(2. * (1. - zmax));
                double Ntot = N1 + N2;
                double r0 = N1 / Ntot;

                double acceptance = 0.;
                while (acceptance < dist(gen))
                {
                    double r = dist(gen);
                    if (r < r0)
                    {
                        z = zmin / (1. - zmin * r * Ntot / 2.);
                        acceptance = .5 / (1. - z);
                    } else
                    {
                        z = 1. - std::exp(-(r * Ntot - N1) / 4.) / 2.;
                        acceptance = .25 / z / z;
                    }
                }
                kt2 = 2 * (1. - z) * z * p.e() / tauf;
                // reject cases where qt2>kt2 for mode=0
                double Num = 0., Den = 0.;
                for (int j = 0; j < Ncolls; j++)
                {
                    if (ts[j] < 0)
                        continue;

                    double q2 = qt2s[j], t = ts[j];
                    if (kt2 > q2)
                        Num += q2 * (1. - std::cos(t / tauf));
                    Den += q2 * (1. - std::cos(t / tauf));
                }
                if (Num / Den < dist(gen))
                    continue;
            } else
            {
                bool ok = false;
                double minimum_q2 = 2 * emin / tauf;
                for (int j = 0; j < Ncolls; j++)
                {
                    if (qt2s[j] > minimum_q2 && ts[j] > 0)
                        ok = true;
                }
                if (!ok)
                    continue;
                while (acceptance < dist(gen))
                {
                    z = zmin * std::pow(zmax / zmin, dist(gen));
                    kt2 = 2 * (1. - z) * z * p.e() / tauf;
                    double num = 0.;
                    for (int j = 0; j < Ncolls; j++)
                        if (kt2 < qt2s[j] && ts[j] > 0)
                            num += 1.;
                    acceptance = num / Ncolls;
                }
            }

            // step 3 correct for the splitting function
            // correct for splitting function
            if (p.id() == 21 && (1 + std::pow(1. - z, 3)) / 2. < dist(gen))
                continue;
            if (p.id() != 21 && (1 + std::pow(1. - z, 2)) / 2. < dist(gen))
                continue;

            // finally, sample phikT2 and compute the deflection of the hard parton lt2
            if (mode == 0)
            {
                // no particular angular structure in the H-T expansion
                phik = 2 * M_PI * dist(gen);
                lt2 = kt2;
            } else
            {
                double Psum = 0.;
                std::vector<double> dP;
                dP.resize(Ncolls);
                for (int j = 0; j < Ncolls; j++)
                {
                    if (kt2 < qt2s[j])
                        Psum += (1. - std::cos(ts[j] / tauf));
                    dP[j] = Psum;
                }
                for (int j = 0; j < Ncolls; j++)
                    dP[j] /= Psum;
                double rc = dist(gen);
                int choice = -1;
                for (int j = 0; j < Ncolls; j++)
                {
                    if (rc < dP[j])
                    {
                        choice = j;
                        break;
                    }
                }
                // sample phik ~ (1+delta cos) / (1+delta^2 + 2 delta cos)
                double delta = std::sqrt(kt2 / qt2s[choice]);
                acceptance = 0.;
                while (acceptance < dist(gen))
                {
                    r = dist(gen);
                    dphiqk = 2. * std::atan(std::tan(M_PI / 2. * r) * (delta + 1) / (delta - 1));
                    acceptance = (1 + delta * std::cos(dphiqk)) / 2.;
                }
                phik = phis[choice] + ((dist(gen) > .5) ? dphiqk : (-dphiqk));
                // lt2 = |(K-q) + q| = |k-q|^2
                lt2 = kt2 + qt2s[choice] + 2. * std::sqrt(kt2 * qt2s[choice]) * std::cos(dphiqk);
            }
            // The virtuality was not allowed to be overlap with the high-Q shower
            if (lt2 > kt2max_now)
                continue;

            // Now, there is a radiation,
            // Hard parton splits into p -> p-k & k
            double kt = std::sqrt(kt2), k0 = z * p.e();
            if (kt > k0)
                continue;
            double kz = std::sqrt(k0 * k0 - kt2);
            Vec4 kmu{kt * std::cos(phik), kt * std::sin(phik), kz, k0};
            kmu.rot(p.theta(), 0.);
            kmu.rot(0., p.phi());
            p.p(p.p() - kmu);
            p.e(std::sqrt(p.pAbs2() + p.m2()));
            kmu.e(std::sqrt(kmu.pAbs2()));

            // the gluon can continue to collide

            std::vector<double> g_qt2s, g_ts, g_phis;
            Coll.sample_all_qt2(21, kmu.e(), L, TA, xB, Q20, g_qt2s, g_ts, g_phis);
            Vec4 Qtot{0., 0., 0., 0.};
            double e0 = kmu.e();
            for (int ig = 0; ig < g_ts.size(); ig++)
            {
                double qt = std::sqrt(g_qt2s[ig]);
                double phi = g_phis[ig];
                Vec4 qmu{qt * std::cos(phi), qt * std::sin(phi), -qt * qt / 4. / e0, 0.0};
                Qtot = Qtot + qmu;
            }
            Qtot.rot(kmu.theta(), 0.);
            Qtot.rot(0., kmu.phi());
            // turn off elastic broading wenbin
            kmu = kmu + Qtot;
            kmu.e(std::sqrt(kmu.pAbs2()));
            kmu = kmu * e0 / kmu.e();

            // update the color if it is a hard gluon
            // first, the spliting process
            int k_col, k_acol;
            // if the gluon forms inside the nuclei,
            // we consider it will lose color correlation with the original parton,
            // and form a new string with beam remnant
            {
                Particle gluon = Particle(21, 201, i, 0, 0, 0, event.nextColTag(),
                                          event.nextColTag(), kmu, 0.0, 0);
                int qid, diqid;
                double mq, mdiq, mn;
                if (dist(gen) < ZoverA)
                { // diquark from a proton
                    diqid = 2101;
                    mdiq = 0.57933;
                    qid = 2;
                    mq = 0.33;
                    mn = 0.93847;
                    if (dist(gen) < 2. / 3.)
                    { // take away a u
                        qid = 2;
                        if (dist(gen) < .75)
                            diqid = 2101;
                        else
                            diqid = 2103;
                    } else
                    { // take away the d
                        qid = 1;
                        diqid = 2203;
                    }
                } else
                { // diquark from a neutron
                    diqid = 2101;
                    mdiq = 0.57933;
                    qid = 1;
                    mq = 0.33;
                    mn = 0.93957;
                    if (dist(gen) < 2. / 3.)
                    { // take away a d
                        qid = 1;
                        if (dist(gen) < .75)
                            diqid = 2101;
                        else
                            diqid = 2103;
                    } else
                    { // take away the u
                        qid = 2;
                        diqid = 1103;
                    }
                }
                double pabs = std::sqrt((mn * mn - std::pow(mq + mdiq, 2)) *
                                        (mn * mn - std::pow(mq - mdiq, 2))) /
                              (2. * mn);
                double costheta = dist(gen) * 2. - 1.;
                double sintheta = std::sqrt(std::max(1. - costheta * costheta, 1e-9));
                double rphi = 2 * M_PI * dist(gen);
                double Nqz = pabs * costheta, Nqx = pabs * sintheta * std::cos(rphi),
                       Nqy = pabs * sintheta * std::sin(rphi);
                Vec4 pq{Nqx, Nqy, Nqz, 0}, pdiq{-Nqx, -Nqy, -Nqz, 0};
                pq.e(std::sqrt(pq.pAbs2() + mq * mq));
                pdiq.e(std::sqrt(pdiq.pAbs2() + mdiq * mdiq));
                Particle recolQ = Particle(qid, 201, i, 0, 0, 0, gluon.acol(), 0, pq, mq, 0);
                Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0, 0, gluon.col(), pdiq, mdiq, 0);
                // time tauf
                recoil_remnants.push_back(gluon);
                x_arr.push_back(tauf * HBARC * kmu.px() / kmu.e());
                y_arr.push_back(tauf * HBARC * kmu.py() / kmu.e());
                z_arr.push_back(tauf * HBARC * kmu.pz() / kmu.e());
                t_arr.push_back(tauf * HBARC);

                recoil_remnants.push_back(recolQ);
                x_arr.push_back(tauf * HBARC * pq.px() / pq.e());
                y_arr.push_back(tauf * HBARC * pq.py() / pq.e());
                z_arr.push_back(tauf * HBARC * pq.pz() / pq.e());
                t_arr.push_back(tauf * HBARC);

                recoil_remnants.push_back(recoldiQ);
                x_arr.push_back(tauf * HBARC * pdiq.px() / pdiq.e());
                y_arr.push_back(tauf * HBARC * pdiq.py() / pdiq.e());
                z_arr.push_back(tauf * HBARC * pdiq.pz() / pdiq.e());
                t_arr.push_back(tauf * HBARC);
            }
        } // Commnent out this loop to turn off the radiation Wenbin // end of the loop

        // Now handles recoil and remannts
        // if there are radiations, recoil goes to radiations
        // else: goes to the hard quark
        int Nrad = frag_gluons.size();
        for (int j = 0; j < Ncolls; j++)
        {
            double qT = std::sqrt(qt2s[j]), phiq = phis[j];
            double qx = qT * std::cos(phiq);
            double qy = qT * std::sin(phiq);
            double qz = -qT * qT / 4. / p.e();
            Vec4 qmu{qx, qy, qz, 0};
            qmu.rot(p.theta(), 0.);
            qmu.rot(0., p.phi());
            // turn off the elastic broading Wenbin
            p.p(p.p() + qmu);
            p.e(std::sqrt(p.pAbs2() + p.m2()));
            int q_col, q_acol;
            // update color
            if (p.id() == 21)
            {
                if (std::rand() % 2 == 0)
                {
                    q_acol = p.acol();
                    p.acol(event.nextColTag());
                    q_col = p.acol();
                } else
                {
                    q_col = p.col();
                    p.col(event.nextColTag());
                    q_acol = p.col();
                }
            } else if (p.id() > 0)
            {
                q_col = p.col();
                p.col(event.nextColTag());
                q_acol = p.col();
            } else
            {
                q_acol = p.acol();
                p.acol(event.nextColTag());
                q_col = p.acol();
            }
            int qid, diqid;
            double mq, mdiq, mn;
            if (dist(gen) < ZoverA)
            { // diquark from a proton
                diqid = 2101;
                mdiq = 0.57933;
                qid = 2;
                mq = 0.33;
                mn = 0.93847;
                if (dist(gen) < 2. / 3.)
                { // take away a u
                    qid = 2;
                    if (dist(gen) < .75)
                        diqid = 2101;
                    else
                        diqid = 2103;
                } else
                { // take away the d
                    qid = 1;
                    diqid = 2203;
                }
            } else
            { // diquark from a neutron
                diqid = 2101;
                mdiq = 0.57933;
                qid = 1;
                mq = 0.33;
                mn = 0.93957;
                if (dist(gen) < 2. / 3.)
                { // take away a d
                    qid = 1;
                    if (dist(gen) < .75)
                        diqid = 2101;
                    else
                        diqid = 2103;
                } else
                { // take away the u
                    qid = 2;
                    diqid = 1103;
                }
            }
            double pabs =
                std::sqrt((mn * mn - std::pow(mq + mdiq, 2)) * (mn * mn - std::pow(mq - mdiq, 2))) /
                (2. * mn);
            double costheta = dist(gen) * 2. - 1.;
            double sintheta = std::sqrt(std::max(1. - costheta * costheta, 1e-9));
            double rphi = 2 * M_PI * dist(gen);
            double Nqz = pabs * costheta, Nqx = pabs * sintheta * std::cos(rphi),
                   Nqy = pabs * sintheta * std::sin(rphi);
            Vec4 pq{Nqx, Nqy, Nqz, 0}, pdiq{-Nqx, -Nqy, -Nqz, 0};
            // decide which object takes the recoil
            if (std::rand() % 2 == 0)
            {
                pq = pq - qmu;
            } else
            {
                pdiq = pdiq - qmu;
            }
            pq.e(std::sqrt(pq.pAbs2() + mq * mq));
            pdiq.e(std::sqrt(pdiq.pAbs2() + mdiq * mdiq));
            Particle recolQ = Particle(qid, 201, i, 0, 0, 0, q_col, 0, pq, mq, 0);
            Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0, 0, q_acol, pdiq, mdiq, 0);

            recoil_remnants.push_back(recolQ);
            x_arr.push_back(ts[j] * HBARC * pq.px() / pq.e());
            y_arr.push_back(ts[j] * HBARC * pq.py() / pq.e());
            z_arr.push_back(ts[j] * HBARC * pq.pz() / pq.e());
            t_arr.push_back(ts[j] * HBARC);

            recoil_remnants.push_back(recoldiQ);
            x_arr.push_back(ts[j] * HBARC * pdiq.px() / pdiq.e());
            y_arr.push_back(ts[j] * HBARC * pdiq.py() / pdiq.e());
            z_arr.push_back(ts[j] * HBARC * pdiq.pz() / pdiq.e());
            t_arr.push_back(ts[j] * HBARC);
        }
        for (auto& p : frag_gluons)
            new_particles.push_back(p);
        for (auto& p : recoil_remnants)
            new_particles.push_back(p);
    }

    int ixyz = 0;
    for (auto& p : new_particles)
    {
        event.append(p.id(), 201, p.col(), p.acol(), p.px(), p.py(), p.pz(), p.e(), p.m());
        p.vProd(p.vProd() + Vec4(x_arr[ixyz] + Rx * HBARC, y_arr[ixyz] + Ry * HBARC,
                                 z_arr[ixyz] + Rz * HBARC, t_arr[ixyz]));
        ixyz++;
    }
}
