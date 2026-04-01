#include "hadronizer.hpp"

#include <cmath>
#include <iostream>

using namespace Pythia8;

namespace {
    void rotate(double px, double py, double pz, double pr[4], int icc) {

        //     input:  (px,py,pz), (wx,wy,wz), argument (i)
        //     output: new (wx,wy,wz)
        //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
        //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

        double E;
        double pt;
        double cosa, sina, cosb, sinb;
        double wx, wy, wz;
        double wx1, wy1, wz1;

        wx = pr[1];
        wy = pr[2];
        wz = pr[3];

        E = std::sqrt(px * px + py * py + pz * pz);
        pt = std::sqrt(px * px + py * py);

        // w = sqrt(wx*wx+wy*wy+wz*wz);

        //  if(pt==0)

        if (pt < 1e-6) {
            cosa = 1;
            sina = 0;
        } else {
            cosa = px / pt;
            sina = py / pt;
        }

        if (E > 1e-6) {
            cosb = pz / E;
            sinb = pt / E;
            if (icc == 1) {
                wx1 = wx * cosb * cosa + wy * cosb * sina - wz * sinb;
                wy1 = -wx * sina + wy * cosa;
                wz1 = wx * sinb * cosa + wy * sinb * sina + wz * cosb;
            } else {
                wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
                wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
                wz1 = -wx * sinb + wz * cosb;
            }
            wx = wx1;
            wy = wy1;
            wz = wz1;
        } else {
            std::cout << "warning: small E in rotation" << std::endl;
        }

        pr[1] = wx;
        pr[2] = wy;
        pr[3] = wz;
    }
} // namespace

Hadronizer::Hadronizer() : pythia(), rd(), gen(rd()), dist(0., 1.) {

    pythia.readString("ProcessLevel:all = off");
    pythia.readString("Print:quiet = on");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    // parton tune and PDF set, please check
    pythia.readString("Tune:pp = 19");
    pythia.readString("PDF:pSet = 12");
    // These two parameters have something to do with the
    // Lund String interative breaking conditions
    // ** We have changed stopMass=0.0 GeV, different from Pythia8 Default
    // too match HERMES FF measurements of pion and Kaon.
    pythia.readString("StringFragmentation:stopMass = 0.0");
    pythia.readString("HadronLevel:mStringMin = 0.5");
    // Set some hadronic & decay specific channels
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("111:mayDecay=off");
    pythia.readString("211:mayDecay=off");
    pythia.readString("321:mayDecay=off");
    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    pythia.readString("PDF:nPDFSetA=0");
    pythia.readString("PDF:nPDFSetB=0");
    // Enabling setting of vertex information.
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("PartonVertex:modeVertex = 2");
    pythia.readString("PartonVertex:ProtonRadius = 0.85");
    pythia.readString("PartonVertex:EmissionWidth = 0.1");
    pythia.readString("Fragmentation:setVertices = on");
    pythia.readString("HadronVertex:mode = 0");
    pythia.readString("HadronVertex:smearOn = on");
    pythia.readString("HadronVertex:xySmear = 0.7");

    pythia.init();
}

/*
This function takes the shower PythiaIn (ep-shower with recoil particles)
and assume it fragments in a nuclear medium with charge Z and atomic mass A
*/
std::vector<Particle> Hadronizer::hadronize(Pythia& pythiaIn,
                                            int Z,
                                            int A,
                                            double Rx,
                                            double Ry,
                                            double Rz) {

    // Define charge fraction
    double ZoverA = Z * 1.0 / A;

    // This vector will store the final hadrons after hadronization, and will be
    // returned to the main function
    std::vector<Particle> FinalStateParticles;
    FinalStateParticles.clear();

    // Get the initial hard parton ID
    int hardid = pythiaIn.event[5].id();

    // Empty the event record of the hadronizer, and prepare to fill it with
    // partons from the shower
    pythia.event.reset();

    // Loop over the partons in the shower PythiaIn, and put them in the hadronizer
    for (int i = 0; i < pythiaIn.event.size(); ++i) {

        auto& particle = pythiaIn.event[i];

        // Only consider final state partons
        // Drop intermediate partons and non-partons
        if (!(particle.isFinal() && particle.isParton())) {            
            continue;
        }

        // Di-quark remnants
        // status() == 63 : outgoing beam remnant from particles produced by beam-remnant treatment
        // idAbs() > 1000 && idAbs() < 3000 : di-quarks
        //   (dd)_1 : 1103
        //   (ud)_0 : 2101
        //   (ud)_1 : 2103
        //   (uu)_1 : 2203
        if (particle.status() == 63 && 1000 < particle.idAbs() && particle.idAbs() < 3000) {

            // valence stuff, the remnants will contain the rest flavor component
            // note that the hard quark has already been sampled according to the
            // the isospin content of the nuclear PDF;
            // *** However, the remanent is generated assuming the rest stuff comes
            // from a proton. Therefore, we need to resample it according to the Z/A
            // ratio this nuclei
            // 1) Decide whether it is from a neutron or proton
            if (dist(gen) < ZoverA) { // From a proton 2212
                if (hardid == 1) {
                    // Produce (uu)_1 : 2203
                    particle.id(2203);
                } else if (hardid == 2) {
                    // Produce (ud)_0 : 2101 and (ud)_1 : 2103 with ratio 3:1
                    if (dist(gen) < 0.75) {
                        particle.id(2101);
                    } else {
                        particle.id(2103);
                    }
                }
            }
            else { // From a neutron 2112
                if (hardid == 1) {
                    // Produce (ud)_0 : 2101 and (ud)_1 : 2103 with ratio 3:1
                    if (dist(gen) < 0.75) {
                        particle.id(2101);
                    } else {
                        particle.id(2103);
                    }
                }
                else if (hardid == 2) {
                    // Produce (dd)_1 : 1103
                    particle.id(1103);
                }
            }
        }

        // For other partons, just put it in the shower
        // status() == 23 : outgoing particles of the hardest subprocess
        pythia.event.append(particle.id(),
                            23,
                            particle.col(),
                            particle.acol(),
                            particle.px(),
                            particle.py(),
                            particle.pz(),
                            particle.e(),
                            particle.m());

        auto& appended = pythia.event.back();
        appended.vProd( particle.vProd() );
    }

    // Assign the space-time to partons
    // Loop over the partons in the hadronizer
    for (int i = 0; i < pythia.event.size(); ++i) {

        auto& particle = pythia.event[i];

        // Only consider partons
        if (particle.isParton()) {

            // Initialize the (x, y, z) position of the parton
            double position[3] = {0.0};

            // Initialize the four-momentum of the parton
            double p0[4] = {0.0};

            // Initialize the four-momentum of the mother
            double p4[4] = {0.0};

            double qt, time_step;
            double timeplus = 0.0;
            int IDmother1, IDmother2;
            int timebreaker = 0;
            int j = i;
            int IDmother0 = j;

            // Compute the formation and spatial information of partons
            while (timebreaker == 0) {

                int IDiii = IDmother0;
                int IDdaughter1, IDdaughter2;

                // 
                if (std::abs(pythia.event[IDiii].status()) == 23 ||
                    std::abs(pythia.event[IDiii].status()) == 21 ||
                    std::abs(pythia.event[IDiii].status()) == 12) {
                    timebreaker = 1;
                }

                // Find the mothers of the parton
                IDmother1 = pythia.event[IDiii].mother1();
                IDmother2 = pythia.event[IDiii].mother2();

                /*
                1) Parton with no mother
                */
                if (IDmother1 == IDmother2 && IDmother1 == 0) {
                    timebreaker = 1;
                }

                /*
                2) Parton is a "carbon copy" of its mother, but with changed momentum
                as a "recoil" effect, e.g. in a shower.
                */
                if (IDmother1 == IDmother2 && IDmother1 > 0) {

                    IDmother0 = IDmother1;

                    // Find the daughters of the mother
                    IDdaughter1 = pythia.event[IDmother0].daughter1();
                    IDdaughter2 = pythia.event[IDmother0].daughter2();

                    if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0) {

                        p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                        p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                        p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                        p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();

                    }

                    if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                        (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
    
                        p4[0] = pythia.event[IDdaughter1].e();
                        p4[1] = pythia.event[IDdaughter1].px();
                        p4[2] = pythia.event[IDdaughter1].py();
                        p4[3] = pythia.event[IDdaughter1].pz();

                    }
                }

                /* 
                3) The "normal" mother case, where it is meaningful to speak of one
                single mother to several products, in a shower or decay.
                */
                if (IDmother1 > 0 && IDmother2 == 0) {

                    IDmother0 = IDmother1;

                    // Find the daughters of the mother
                    IDdaughter1 = pythia.event[IDmother0].daughter1();
                    IDdaughter2 = pythia.event[IDmother0].daughter2();

                    // 
                    if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0) {

                        // Compute the four-momentum of the mother by summing over the two daughters
                        p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                        p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                        p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                        p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();

                    }
                }

                /*
                4,5) For abs(status) = 81 - 86: primary hadrons produced from the
                fragmentation of a string spanning the range from mother1 to mother2,
                so that all partons in this range should be considered mothers; and
                analogously for abs(status) = 101 - 106, the formation of R-hadrons.
                Or, particles with two truly different mothers, in particular the
                particles emerging from a hard 2 → n interaction.
                */
                if (IDmother1 < IDmother2 && IDmother1 > 0 && IDmother2 > 0) {

                    if (pythia.event[IDmother1].e() > pythia.event[IDmother2].e()) {
                        IDmother0 = IDmother1;
                    }

                    if (pythia.event[IDmother1].e() <= pythia.event[IDmother2].e()) {
                        IDmother0 = IDmother2;
                    }

                    IDdaughter1 = pythia.event[IDmother0].daughter1();
                    IDdaughter2 = pythia.event[IDmother0].daughter2();

                    if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0) {

                        p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                        p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                        p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                        p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();

                    }

                    if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                        (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
    
                        p4[0] = pythia.event[IDdaughter1].e();
                        p4[1] = pythia.event[IDdaughter1].px();
                        p4[2] = pythia.event[IDdaughter1].py();
                        p4[3] = pythia.event[IDdaughter1].pz();

                    }
                }

                /*
                6) Particles with two truly different mothers, notably for the
                special case that two nearby partons are joined together into a
                status 73 or 74 new parton, in the g + q → q case the q is made
                first mother to simplify flavour tracing.
                */
                if (IDmother2 < IDmother1 && IDmother1 > 0 && IDmother2 > 0) {

                    if (pythia.event[IDmother1].e() > pythia.event[IDmother2].e()) {
                        IDmother0 = IDmother1;
                    }

                    if (pythia.event[IDmother1].e() <= pythia.event[IDmother2].e()) {
                        IDmother0 = IDmother2;
                    }

                    IDdaughter1 = pythia.event[IDmother0].daughter1();
                    IDdaughter2 = pythia.event[IDmother0].daughter2();

                    if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0) {

                        p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                        p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                        p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                        p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();

                    }

                    if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                        (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {

                        p4[0] = pythia.event[IDdaughter1].e();
                        p4[1] = pythia.event[IDdaughter1].px();
                        p4[2] = pythia.event[IDdaughter1].py();
                        p4[3] = pythia.event[IDdaughter1].pz();

                    }
                }

                double x_split = pythia.event[IDiii].e() / p4[0];

                if (x_split > 1) {
                    x_split = 1.0 / x_split; // revise the mother and daughter
                }

                p0[0] = pythia.event[IDiii].e();
                p0[1] = pythia.event[IDiii].px();
                p0[2] = pythia.event[IDiii].py();
                p0[3] = pythia.event[IDiii].pz();
                
                rotate(p4[1], p4[2], p4[3], p0, 1);
                qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                rotate(p4[1], p4[2], p4[3], p0, -1);
                double kt_daughter = qt;

                if (kt_daughter > 1.e-10) {
                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                pow(kt_daughter, 2) * HBARC;
                    timeplus = timeplus + time_step;
                    position[0] = position[0] + time_step * p4[1] / p4[0];
                    position[1] = position[1] + time_step * p4[2] / p4[0];
                    position[2] = position[2] + time_step * p4[3] / p4[0];

                }
            }

            particle.vProd(particle.vProd() + Vec4(position[0] + Rx * HBARC,
                                                   position[1] + Ry * HBARC,
                                                   position[2] + Rz * HBARC,
                                                   timeplus));
        }
    }

    // Perform hadronization
    const bool ok = pythia.next();
    if (!ok) {
        std::cerr << "Hadronizer: pythia.next() failed\n";
        return {};
    }

    // Return the final state hadrons after hadronization
    for (int i = 0; i < pythia.event.size(); ++i) {

        auto& particle = pythia.event[i];

        if (particle.isFinal()) {
            FinalStateParticles.push_back(particle);

        }
    }

    return FinalStateParticles;
}