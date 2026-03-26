#include "Pythia8/Pythia.h"
#include "cli.hpp"

using namespace Pythia8;
#include <algorithm>
#include <chrono>
#include <cstdlib> // for rand() and srand()
#include <ctime>   // for time()
#include <fstream>
#include <random>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>


std::random_device rrd;  // Non-deterministic generator
std::mt19937 gen(rrd()); // Mersenne Twister engine

// Define a distribution (range 0.0 to 1.0)
std::uniform_real_distribution<> Ran_gen(0.0, 1.0);


namespace fs = std::filesystem;


static std::string zero_pad_int(int64_t value, int width = 8) {
    std::ostringstream os;
    os << std::setw(width) << std::setfill('0') << value;
    return os.str();
}

static fs::path shard_dir_for_event(const fs::path& base_events_dir, int64_t event_id, int64_t chunk_size) {
    const int64_t shard_begin = (event_id / chunk_size) * chunk_size;
    const int64_t shard_end   = shard_begin + chunk_size - 1;

    std::ostringstream name;
    name << "events_" << zero_pad_int(shard_begin) << "-" << zero_pad_int(shard_end);
    return base_events_dir / name.str();
}


void rotate(double px, double py, double pz, double pr[4], int icc)
{
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

    E = sqrt(px * px + py * py + pz * pz);
    pt = sqrt(px * px + py * py);

    // w = sqrt(wx*wx+wy*wy+wz*wz);

    //  if(pt==0)
    if (pt < 1e-6)
    {
        cosa = 1;
        sina = 0;
    } else
    {
        cosa = px / pt;
        sina = py / pt;
    }

    if (E > 1e-6)
    {
        cosb = pz / E;
        sinb = pt / E;
        if (icc == 1)
        {
            wx1 = wx * cosb * cosa + wy * cosb * sina - wz * sinb;
            wy1 = -wx * sina + wy * cosa;
            wz1 = wx * sinb * cosa + wy * sinb * sina + wz * cosb;
        } else
        {
            wx1 = wx * cosa * cosb - wy * sina + wz * cosa * sinb;
            wy1 = wx * sina * cosb + wy * cosa + wz * sina * sinb;
            wz1 = -wx * sinb + wz * cosb;
        }
        wx = wx1;
        wy = wy1;
        wz = wz1;
    } else
    {
        cout << "warning: small E in rotation" << endl;
    }

    pr[1] = wx;
    pr[2] = wy;
    pr[3] = wz;
}

double fermi_distribution(double r, double R_WS, double a_WS)
{
    double f = 1. / (1. + exp((r - R_WS) / a_WS));
    return (f);
}

// S. Chakraborty et al., Physical Review C 107, 064318 (2023)
double sample_r_from_woods_saxon(double R_WS, double a_WS)
{
    double rmaxCut = R_WS + 10. * a_WS;
    double r = 0.;
    do
    {
        r = rmaxCut * pow(Ran_gen(gen), 1.0 / 3.0);
    } while (Ran_gen(gen) > fermi_distribution(r, R_WS, a_WS));
    return (r);
}

double samplePointInSphere(double R)
{
    double u = Ran_gen(gen);
    double r = R * cbrt(u);
    return (r);
}

// ==============================================================

// A separate Pythia instance that only handles hadronization
class Hadronizer
{
  public:
    // Constructor, set Pythia, random generator
    Hadronizer() : pythia(), rd(), gen(rd()), dist(0., 1.)
    {
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
        pythia.readString("StringFragmentation:stopMass = 0.0");         // Wenbin lower limit of string mass
        pythia.readString("HadronLevel:mStringMin = 0.5"); // Wenbin
        // Set some hadronic & decay specific channls
        pythia.readString("HadronLevel:Decay = off");
        pythia.readString("111:mayDecay=off");
        pythia.readString("211:mayDecay=off");
        pythia.readString("321:mayDecay=off");

        // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
        pythia.readString("PDF:nPDFSetA=0");
        pythia.readString("PDF:nPDFSetB=0");
        // Enabling setting of vertex information.
        pythia.readString("PartonVertex:setVertex = on");
        pythia.readString("PartonVertex:protonRadius = 0.7");
        pythia.readString("PartonVertex:emissionWidth = 0.1");
        pythia.readString("Fragmentation:setVertices = on");

        pythia.init();
    }

    // This function takes the shower PythiaIn (ep-shower with recoil particles)
    // and assume it fragments in a nuclear medium with charge Z and atomic mass A
    std::vector<Particle> hadronize(Pythia& pythiaIn, int Z, int A, double Rx, double Ry, double Rz)
    {

        // Define charge fraction
        double ZoverA = Z * 1.0 / A;

        // This vector will store the final hadrons after hadronization, and will be returned to the main function
        std::vector<Particle> FinalStateParticles;
        FinalStateParticles.clear();

        // Get the initial hard parton ID
        int hardid = pythiaIn.event[5].id();

        // Empty the event record of the hadronizer, and prepare to fill it with partons from the shower
        pythia.event.reset();

        // Loop over the partons in the shower PythiaIn, and put them in the hadronizer
        for (int i = 0; i < pythiaIn.event.size(); ++i)
        {

            auto& particle = pythiaIn.event[i];

            // Only consider final state partons
            // Drop intermediate partons and non-partons
            if (!(particle.isFinal() && particle.isParton()))
                continue;
            
            // Di-quark remnants
            // status() == 63 : outgoing beam remnant from particles produced by beam-remnant treatment
            // idAbs() > 1000 && idAbs() < 3000 : di-quarks
            //   (dd)_1 : 1103
            //   (ud)_0 : 2101
            //   (ud)_1 : 2103
            //   (uu)_1 : 2203
            if (particle.status() == 63 && 1000 < particle.idAbs() && particle.idAbs() < 3000)
            {
                // valence stuff, the remnants will contain the rest flavor component
                // note that the hard quark has already been sampled according to the
                // the isospin content of the nuclear PDF;
                // *** However, the remanent is generated assuming the rest stuff comes
                // from a proton. Therefore, we need to resample it according to the Z/A
                // ratio this nuclei
                // 1) Decide whether it is from a neutron or proton
                if (dist(gen) < ZoverA) // From a proton 2212
                { 
                    if (hardid == 1)
                    {
                        // Produce (uu)_1 : 2203
                        particle.id(2203);
                    }
                    if (hardid == 2)
                    { 
                        // Produce (ud)_0 : 2101 and (ud)_1 : 2103 with ratio 3:1
                        if (dist(gen) < 0.75)
                            particle.id(2101);
                        else
                            particle.id(2103);
                    }
                } else // From a neutron 2112
                { 
                    if (hardid == 1)
                    {
                        // Produce (ud)_0 : 2101 and (ud)_1 : 2103 with ratio 3:1
                        if (dist(gen) < 0.75)
                            particle.id(2101);
                        else
                            particle.id(2103);
                    }
                    if (hardid == 2)
                    { 
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
        }

        // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
        // Assign the space-time to partons
        // Loop over the partons in the hadronizer
        for (int i = 0; i < pythia.event.size(); ++i)
        {

            auto& particle = pythia.event[i];

            // Only consider partons
            if (particle.isParton())
            {
                
                // Initialize the (x, y, z) position of the parton
                double position[3] = {0.0};

                // Initialize the four-momentum of the parton
                double p0[4] = {0.0};

                // Initialize the four-momentum of the mother
                double p4[4] = {0.0};


                double qt, time_step;
                double timeplus = 0.0;
                int IDmom1, IDmom2;
                int timebreaker = 0;
                int j = i;
                int IDmom0 = j;

                // Compute the formation and spatial information of partons
                while (timebreaker == 0)
                {
                    int IDiii = IDmom0;
                    if (std::abs(pythia.event[IDiii].status()) == 23 ||
                        std::abs(pythia.event[IDiii].status()) == 21 ||
                        std::abs(pythia.event[IDiii].status()) == 12)
                        timebreaker = 1;

                    // Find the mothers of the parton
                    IDmom1 = pythia.event[IDiii].mother1();
                    IDmom2 = pythia.event[IDiii].mother2();

                    // 1) Parton with no mother
                    if (IDmom1 == IDmom2 && IDmom1 == 0)
                        timebreaker = 1;

                    // 2) Parton is a "carbon copy" of its mother, but with changed momentum as a "recoil" effect, e.g. in a shower
                    if (IDmom1 == IDmom2 && IDmom1 > 0)
                        IDmom0 = IDmom1;

                    // 3) The "normal" mother case, where it is meaningful to speak of one single mother to several products, in a shower or decay;
                    if (IDmom1 > 0 && IDmom2 == 0)
                    {

                        IDmom0 = IDmom1;

                        // Find the daughters of the mother
                        double IDdaughter1 = pythia.event[IDmom0].daughter1();
                        double IDdaughter2 = pythia.event[IDmom0].daughter2();

                        // 
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0)
                        {

                            // Compute the four-momentum of the mother by summing over the two daughters
                            p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                            p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                            p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                            p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();

                            //
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split;

                            // Compute the four-momentum of the parton itself
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();

                            // 
                            rotate(p4[1], p4[2], p4[3], p0, 1);

                            // 
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);

                            // 
                            rotate(p4[1], p4[2], p4[3], p0, -1);

                            //
                            double kt_daughter = qt;

                            //
                            if (x_split < 0.5)
                            {
                                // 
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }
                    }

                    // 4) For abs(status) = 81 - 86: primary hadrons produced from the fragmentation of a
                    // string spanning the range from mother1 to mother2, so that all partons in this
                    // range should be considered mothers; and analogously for abs(status) = 101 - 106,
                    // the formation of R-hadrons. Or, particles with two truly different mothers, in
                    // particular the particles emerging from a hard 2 → n interaction.
                    if (IDmom1 < IDmom2 && IDmom1 > 0 && IDmom2 > 0)
                    {
                        if (pythia.event[IDmom1].e() > pythia.event[IDmom2].e())
                        {
                            IDmom0 = IDmom1;
                        }
                        if (pythia.event[IDmom1].e() <= pythia.event[IDmom2].e())
                        {
                            IDmom0 = IDmom2;
                        }
                        double IDdaughter1 = pythia.event[IDmom0].daughter1();
                        double IDdaughter2 = pythia.event[IDmom0].daughter2();
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0)
                        {
                            p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                            p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                            p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                            p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split;

                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            rotate(p4[1], p4[2], p4[3], p0, 1);
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);
                            double kt_daughter = qt;

                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }
                        if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                            (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1))
                        {
                            p4[0] = pythia.event[IDdaughter1].e();
                            p4[1] = pythia.event[IDdaughter1].px();
                            p4[2] = pythia.event[IDdaughter1].py();
                            p4[3] = pythia.event[IDdaughter1].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split; // revise the mother and daughter
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            rotate(p4[1], p4[2], p4[3], p0, 1); // rotate into the
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);
                            double kt_daughter = qt;

                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }
                    }

                    // 5) Particles with two truly different mothers, notably for the special
                    // case that two nearby partons are joined together into a status 73 or 74
                    // new parton, in the g + q → q case the q is made first mother to
                    // simplify flavour tracing.
                    if (IDmom1 > IDmom2 && IDmom1 > 0 && IDmom2 > 0)
                    {
                        if (pythia.event[IDmom1].e() > pythia.event[IDmom2].e())
                        {
                            IDmom0 = IDmom1;
                        }
                        if (pythia.event[IDmom1].e() <= pythia.event[IDmom2].e())
                        {
                            IDmom0 = IDmom2;
                        }
                        double IDdaughter1 = pythia.event[IDmom0].daughter1();
                        double IDdaughter2 = pythia.event[IDmom0].daughter2();
                        if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0)
                        {
                            p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                            p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                            p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                            p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split;
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            rotate(p4[1], p4[2], p4[3], p0, 1);
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);

                            double kt_daughter = qt;
                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }

                        if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                            (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1))
                        {
                            p4[0] = pythia.event[IDdaughter1].e();
                            p4[1] = pythia.event[IDdaughter1].px();
                            p4[2] = pythia.event[IDdaughter1].py();
                            p4[3] = pythia.event[IDdaughter1].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split; // revise the mother and daughter
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            // double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
                            // double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));
                            rotate(p4[1], p4[2], p4[3], p0, 1); // rotate into the
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);
                            double kt_daughter = qt;
                            // double Q2 = 1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }

                    }


                    if (IDmom1 == IDmom2 && IDmom1 > 0)
                    {
                        IDmom0 = IDmom1;
                        int IDdaughter1 = pythia.event[IDmom0].daughter1();
                        int IDdaughter2 = pythia.event[IDmom0].daughter2();

                        if (IDdaughter1 != IDdaughter2 && IDdaughter1 > 0 && IDdaughter2 > 0)
                        {
                            p4[0] = pythia.event[IDdaughter1].e() + pythia.event[IDdaughter2].e();
                            p4[1] = pythia.event[IDdaughter1].px() + pythia.event[IDdaughter2].px();
                            p4[2] = pythia.event[IDdaughter1].py() + pythia.event[IDdaughter2].py();
                            p4[3] = pythia.event[IDdaughter1].pz() + pythia.event[IDdaughter2].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split; // revise the mother and daughter
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            rotate(p4[1], p4[2], p4[3], p0, 1); // rotate into the
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);
                            double kt_daughter = qt;
                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }

                        if ((IDdaughter1 > 0 && IDdaughter2 == 0) ||
                            (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1))
                        {
                            p4[0] = pythia.event[IDdaughter1].e();
                            p4[1] = pythia.event[IDdaughter1].px();
                            p4[2] = pythia.event[IDdaughter1].py();
                            p4[3] = pythia.event[IDdaughter1].pz();
                            double x_split = pythia.event[IDiii].e() / p4[0];
                            if (x_split > 1)
                                x_split = 1.0 / x_split; // revise the mother and daughter
                            p0[0] = pythia.event[IDiii].e();
                            p0[1] = pythia.event[IDiii].px();
                            p0[2] = pythia.event[IDiii].py();
                            p0[3] = pythia.event[IDiii].pz();
                            rotate(p4[1], p4[2], p4[3], p0, 1); // rotate into the
                            qt = sqrt(p0[1] * p0[1] + p0[2] * p0[2]);
                            rotate(p4[1], p4[2], p4[3], p0, -1);
                            double kt_daughter = qt;
                            if (x_split < 0.5)
                            {
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            } else
                            {
                                x_split = 1. - x_split;
                                if (kt_daughter > 1.e-10)
                                {
                                    time_step = 2.0 * p4[0] * x_split * (1 - x_split) /
                                                pow(kt_daughter, 2) * HBARC;
                                    timeplus = timeplus + time_step;
                                    position[0] = position[0] + time_step * p4[1] / p4[0];
                                    position[1] = position[1] + time_step * p4[2] / p4[0];
                                    position[2] = position[2] + time_step * p4[3] / p4[0];
                                }
                            }
                        }
                    }
                }

                // Compute shift in position
                particle.vProd(particle.vProd() + Vec4(position[0] + Rx * HBARC,
                                                       position[1] + Ry * HBARC,
                                                       position[2] + Rz * HBARC, timeplus));
            }
        }

        // Do Hadronization
        pythia.next();

        // Return only final state particles.
        for (int i = 0; i < pythia.event.size(); ++i)
        {
            auto& particle = pythia.event[i];
            if (particle.isFinal())
            {
                FinalStateParticles.push_back(particle);
            }
        }
        return FinalStateParticles;
    }

  private:
    Pythia pythia;
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen;      // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dist;
};

// Helper function to set Pythia settings
template <class T> void add_arg(Pythia& pythia, std::string name, T value)
{
    std::stringstream ss;
    ss << name << " = " << value;
    std::cout << ss.str() << std::endl;
    pythia.readString(ss.str());
}

// You can put whatever triggers of Pythia events here for e-p.
bool trigger(Pythia& pythia)
{
    Vec4 pProton = pythia.event[1].p(); // four-momentum of proton
    Vec4 peIn = pythia.event[4].p();    // incoming electron
    Vec4 peOut = pythia.event[6].p();   // outgoing electron
    Vec4 pGamma = peIn - peOut;         // virtual boson photon/Z^0/W^+-

    // Q2, W2, Bjorken x, y.
    double ymin = 0.10;
    double ymax = 0.85;
    double xBmin = 0.023;
    double xBmax = 0.6;
    double Q2min = 1.0; // GeV^2
    // double Q2max = 6.0; // GeV^2
    double Wmin = std::sqrt(10.0); // GeV
    // double nuMin = 6.0; // GeV

    double nu = (pProton * pGamma) / std::sqrt(pProton * pProton);  // photon energy in target rest frame
    double Q2 = -pGamma.m2Calc();                                   // hard scale square
    double W2 = (pProton + pGamma).m2Calc();                        // invariant mass square of hadronic final state
    double W = std::sqrt(W2);
    double xBj = Q2 / (2. * pProton * pGamma);                       // Bjorken x
    double y = (pProton * pGamma) / (pProton * peIn);               // inelasticity

    // Apply trigger conditions
    return (ymin < y) && (y < ymax) && (xBmin < xBj) && (xBj < xBmax) && (Wmin < std::sqrt(W2)) && (Q2min < Q2);
}

// Output function, we need the final-particle list and kinematics of the original event
// to compute Q, x, etc.
void Output(int32_t eventNumber, int Z, int A, Pythia& pythia,
            std::vector<Particle>& plist,
            std::ofstream& F,
            std::ofstream& M)
{
    // int32_t eventNumber = Ntriggered - 1; // Start event numbering from 0
    // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = pythia.event[1].p();
    Vec4 peIn = pythia.event[4].p();
    Vec4 peOut = pythia.event[6].p();
    Vec4 pGamma = peIn - peOut;
    // Q2, W2, Bjorken x, y.
    double nu = (pProton * pGamma) / std::sqrt(pProton * pProton);  // photon energy in target rest frame
    double Q2 = -pGamma.m2Calc();                                   // hard scale square
    double W2 = (pProton + pGamma).m2Calc();                        // invariant mass square of hadronic final state
    double W = std::sqrt(W2);
    double xB = Q2 / (2. * pProton * pGamma);                       // Bjorken x
    double y = (pProton * pGamma) / (pProton * peIn);               // inelasticity

    // ===============================
    // Write per-event DIS metadata
    // ===============================
    M << std::setprecision(16);
    M << "{\n";
    M << "  \"event\": " << eventNumber << ",\n";
    M << "  \"Z\": " << Z << ", \"A\": " << A << ",\n";
    M << "  \"xB\": " << xB << ",\n";
    M << "  \"Q2\": " << Q2 << ",\n";
    M << "  \"y\": "  << y  << ",\n";
    M << "  \"nu\": " << nu << ",\n";
    M << "  \"P4\": [" << pProton.px() << ", " << pProton.py() << ", "
                    << pProton.pz() << ", " << pProton.e() << "],\n";
    M << "  \"q4\": [" << pGamma.px() << ", " << pGamma.py() << ", "
                    << pGamma.pz() << ", " << pGamma.e() << "]\n";
    M << "}\n";


    double theta = -pGamma.theta();
    double phi = -pGamma.phi();

    Vec4 pProtonB = pProton;
    pProtonB.rot(0, phi);
    pProtonB.rot(theta + M_PI, 0);
    Vec4 pGamma_test = pGamma;
    pGamma_test.rot(0, phi);
    pGamma_test.rot(theta + M_PI, 0);
    double beta = -pGamma_test.e() / pGamma_test.pz();
    double gamma = 1. / (std::sqrt(1. - beta * beta));
    double P0B = gamma * pProtonB.e() + gamma * beta * pProtonB.pz();

    Vec4 pCoM = pGamma + pProton;
    Vec4 pPhoton_in_com = pGamma;

    pPhoton_in_com.bstback(pCoM);

    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    double gamma_com_theta = pPhoton_in_com.theta();
    double gamma_com_phi = pPhoton_in_com.phi();
    double vz = pPhoton_in_com.e() / std::sqrt(pPhoton_in_com.e() * pPhoton_in_com.e() + Q2);

    // Count number of hadrons in the shower
    int Count = 0;
    for (auto& p : plist)
    {
        if (p.isFinal())
            Count++;
    }

    // Output the event header (OSCAR-style)
    F << "# event " << eventNumber << " out " << Count << "\n";

    int32_t particleIndex = 0; // Internal ID for particles
    for (auto& p : plist)
    {
        // assign 1fm for hadronization process
        Vec4 local_pos = {p.xProd(), p.yProd(), p.zProd(), p.tProd()};
        Vec4 pCoM = p.p();
        Vec4 phadron = p.p();
        phadron.bstback(pCoM); // JORDI, try this off
        local_pos.bstback(pCoM);
        // Asign the formation time 1 fm in its local rest frame
        Vec4 formation_4 = {0.00001 * phadron.px() / phadron.e() + local_pos.px(),
                            0.00001 * phadron.py() / phadron.e() + local_pos.py(),
                            0.00001 * phadron.pz() / phadron.e() + local_pos.pz(),
                            0.00001 + local_pos.e()}; // this is time
        formation_4.bst(pCoM);
        double t_hadron = formation_4.e();
        double x_hadron = formation_4.px();
        double y_hadron = formation_4.py();
        double z_hadron = formation_4.pz();
        // Particle properties
        double MASS = p.m();     // Mass of the hadron in GeV
        double p0 = p.e();       // Energy of the hadron in GeV
        int PDGID = p.id();      // PDGID of the hadron
        int PID = particleIndex; // Internal ID
        int CHARGE = p.charge(); // Charge of the hadron
        if (p.isFinal())
        {
            F << std::fixed << std::setprecision(2) << t_hadron << " " << std::fixed
              << std::setprecision(5) << x_hadron << " " << y_hadron << " " << z_hadron << " "
              << MASS << " " << p0 << " " << p.px() << " " << p.py() << " " << p.pz() << " "
              << PDGID << " " << PID << " " << CHARGE << std::endl;
            particleIndex++;
        }
    }

    // Output nucleons
    double RA = 1.2 * pow(A, 1. / 3.);                   // fm
    // std::srand(static_cast<unsigned int>(std::time(0))); // Seed the random number generator
    int randomBinary = std::rand() % 2;                  // Generate 0 or 1

    // Output proton data
    int totalProtons = Z - randomBinary;
    for (auto i = 0; i < totalProtons; ++i)
    {
        double rr = samplePointInSphere(RA);
        double phi = 2. * M_PI * Ran_gen(gen);
        double costheta = 1. - 2. * Ran_gen(gen);
        double sintheta = std::sqrt(1. - costheta * costheta);
        double rx = rr * sintheta * std::cos(phi);
        double ry = rr * sintheta * std::sin(phi);
        double rz = rr * costheta;
        double t = 0.0;  // Formation time
        double px = 0.0; // Momentum in x   // TODO: JORDI ASK WENBIN !!! (FERMI MOMENTUM)
        double py = 0.0; // Momentum in y
        double pz = 0.0; // Momentum in z
        // Particle properties
        double MASS = 0.938272;                                           // Proton mass in GeV
        double p0 = std::sqrt(MASS * MASS + px * px + py * py + pz * pz); // Energy of proton in GeV
        int PDGID = 2212;                                                 // Proton PDGID
        int PID = particleIndex;                                          // Internal ID
        int CHARGE = 1;                                                   // Charge of proton
        F << std::fixed << std::setprecision(2) << t << " " << std::fixed << std::setprecision(5)
          << rx << " " << ry << " " << rz << " " << MASS << " " << p0 << " " << px << " " << py
          << " " << pz << " " << PDGID << " " << PID << " " << CHARGE << std::endl;
        particleIndex++;
    }

    // Output neutron data
    int totalNeutrons = A - Z + randomBinary;
    for (auto i = 0; i < totalNeutrons; ++i)
    {
        // for (auto i =0; i<A-Z-1 + randomBinary; ++i) {
        double rr = samplePointInSphere(RA);
        double phi = 2. * M_PI * Ran_gen(gen);
        double costheta = 1. - 2. * Ran_gen(gen);
        double sintheta = std::sqrt(1. - costheta * costheta);
        double rx = rr * sintheta * std::cos(phi);
        double ry = rr * sintheta * std::sin(phi);
        double rz = rr * costheta;
        double t = 0.0;  // Formation time
        double px = 0.0; // Momentum in x
        double py = 0.0; // Momentum in y
        double pz = 0.0; // Momentum in z
        // Particle properties
        double MASS = 0.939565; // Neutron mass in GeV
        double p0 =
            std::sqrt(MASS * MASS + px * px + py * py + pz * pz); // Energy of neutron in GeV
        int PDGID = 2112;                                         // Neutron PDGID
        int PID = particleIndex;                                  // Internal ID
        int CHARGE = 0;                                           // Charge of neutron
        F << std::fixed << std::setprecision(2) << t << " " << std::fixed << std::setprecision(5)
          << rx << " " << ry << " " << rz << " " << MASS << " " << p0 << " " << px << " " << py
          << " " << pz << " " << PDGID << " " << PID << " " << CHARGE << std::endl;
        particleIndex++;
    }

    F << "# event " << eventNumber << " end 0" << std::endl;
}

// low-Q2 medium correction (a Monte Carlo version of the the modified FF model)
class Modified_FF
{
  public:
    Modified_FF(int mode_, int Z_, int A_, double K, double n, double lambda, std::string TablePath)
        : mode(mode_), Z(Z_), A(A_), ZoverA(Z * 1. / A), Coll(K, n, lambda), eHIJING_Geometry(A, Z),
          rd(), gen(rd()), dist(0., 1.)
    {
        Coll.Tabulate(TablePath);
    };
    void sample_FF_partons(Event& event, double& Rx, double& Ry, double& R);

  private:
    int mode, Z, A;
    double ZoverA;
    EHIJING::MultipleCollision Coll;
    EHIJING::NuclearGeometry eHIJING_Geometry;
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen;      // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dist;
};

int main(int argc, char* argv[])
{

    const RunConfig cfg = parse_args(argc, argv);

    const int nEvents = cfg.nEvents;
    const int Z = cfg.Z;
    const int A = cfg.A;
    const int mode = cfg.mode;
    const double K = cfg.K;
    const std::string& tableDir = cfg.tableDir;
    const std::string& outDir = cfg.runDir;
    const std::string& configFile = cfg.configFile;
    const int64_t first_event_id = cfg.firstEventId;
    const int64_t chunk_size = cfg.chunkSize;
    const uint32_t seed = cfg.seed;

    // Ensure output directories exist
    try // Attempt to create the directories if they do not exist
    {
        fs::create_directories(fs::path(outDir));
        fs::create_directories(fs::path(tableDir));
    } catch (const fs::filesystem_error& e) // Handle any errors that occur during directory creation
    {
        std::cerr << "ERROR: cannot create output/table directories:\n  " << e.what() << "\n";
        return 3;
    }

    // Build the nucleus ID used by Pythia & the PDF (the isospin effect)
    const int iNuclei = 100000000 + Z * 10000 + A * 10;

    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    // Shadowing effect:
    //   nPDFset=0: only isospin
    //   nPDFset = 1: EPS09 LO
    //   nPDFset = 2: EPS09 NLO
    //   nPDFset = 2: EPPS16 NLO
    // We will use only isospin for deuteron,
    // and EPPS16 NLO for heavier nucleus
    int nPDFset = 0; // (A>2)?3:0;

    // Initialize the Pythia instance for hadronization
    Hadronizer hadronizer;

    // Initialize the eHIJING-Pythia for high-Q parton shower in medium
    Pythia pythia;

    // Read settings from the config file and apply them to the Pythia instance
    pythia.readFile(configFile);
    
    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    // Make Pythia deterministic
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    // ALSO seed your own RNGs deterministically (important!)
    gen.seed(seed);
    std::srand(seed);

    // Set Pythia and eHIJING settings from command line arguments
    // Note: Command line arguments will override settings in the config file if there are conflicts
    add_arg<int>(pythia, "PDF:nPDFSetA", nPDFset);
    add_arg<int>(pythia, "PDF:nPDFBeamA", iNuclei);
    add_arg<int>(pythia, "eHIJING:Mode", mode);
    add_arg<int>(pythia, "eHIJING:AtomicNumber", A);
    add_arg<int>(pythia, "eHIJING:ChargeNumber", Z);
    add_arg<double>(pythia, "eHIJING:Kfactor", K);
    add_arg<std::string>(pythia, "eHIJING:TablePath", tableDir);

    // Prepare Pythia for event generation
    pythia.init();

    // Initialize the modified FF module for medium corrections to the parton shower
    Modified_FF MFF(mode, Z, A, K, pythia.settings.parm("eHIJING:xG-n"),
                    pythia.settings.parm("eHIJING:xG-lambda"), tableDir);

    // Define counter for triggered events
    int nTriggered = 0;

    // Define counter for failed events
    int nFailed = 0;

    // Define counter for total events
    int nTotal = 0;
    
    // Begin event loop
    while (nTriggered < nEvents)
    {

        // Add to total event count
        nTotal++;

        // Skip event if it fails to generate
        if (!pythia.next())
        {
            // Add to failed event count
            nFailed++;
            continue;
        };

        // Skip event if it doesn't satisfy the trigger conditions
        if (!trigger(pythia))
            continue;

        // Event passed the trigger, add to triggered event count
        nTriggered++;

        // Initialize the (x, y, z) size
        double Rx, Ry, Rz;

        // Modify the final shower with low-Q^2 medium corrections
        MFF.sample_FF_partons(pythia.event, Rx, Ry, Rz);

        // Put the parton-level event into the separate hadronizer
        auto hadronizerEvent = hadronizer.hadronize(pythia, Z, A, Rx, Ry, Rz);

        // Set event ID
        const int64_t event_id = first_event_id + (nTriggered - 1);

        // Set output paths for this event, organized by shard directories
        const fs::path base_events_dir(outDir);
        const fs::path shard_dir = shard_dir_for_event(base_events_dir, event_id, chunk_size);

        // Create the shard directory if it does not exist
        try 
        {
            fs::create_directories(shard_dir);
        } catch (const fs::filesystem_error& e) 
        {
            std::cerr << "ERROR: cannot create shard directory:\n  " << shard_dir << "\n  " << e.what() << std::endl;
            return 3;
        }

        // Set the zero-padded event ID string for file naming
        const std::string event_str = zero_pad_int(event_id);

        // Set output file paths for this event
        const fs::path event_path = shard_dir / ("event_" + event_str + ".oscar");
        const fs::path meta_path  = shard_dir / ("event_" + event_str + ".meta.json");

        // Open OSCAR output file
        std::ofstream fout_event(event_path);
        if (!fout_event) 
        {
            std::cerr << "ERROR: cannot open output file: " << event_path << std::endl;
            return 1;
        }

        // Open metadata output file
        std::ofstream fout_meta(meta_path);
        if (!fout_meta) 
        {
            std::cerr << "ERROR: cannot open metadata file: " << meta_path << std::endl;
            return 1;
        }

        // Write headers for the OSCAR event file
        fout_event << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n";
        fout_event << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none\n";

        // Write the event data and metadata to the respective files
        Output(event_id, Z, A, pythia, hadronizerEvent, fout_event, fout_meta);

        // Close the output files
        fout_event.close();
        fout_meta.close();

    }

    // Done
    return 0;
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
