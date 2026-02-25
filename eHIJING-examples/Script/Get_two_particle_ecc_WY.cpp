#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
int main(int argc, char* argv[]) {
        // Load all data
    std::vector<double> xB, Q2, pid, chg, pT, eta, phi, theta, zh;
    std::vector<int> breaks;
    char** stemp1;
    double xBd, Q2d, pidd, chgd, pTd, etad, phid, thetad, zhd;
    int Count;
    int Neves = 0;
    char inputfile[128];
    sprintf(inputfile, "%s", argv[1]);
    FILE* infile;
    infile = fopen(inputfile,"r");
    breaks.push_back(0);
    int conut_particle = 0;
    while (!feof(infile)) {
        fscanf(infile,"%s %lf %lf %d\n",&stemp1, &Q2d, &xBd, &Count);
        for (int ip = 0; ip < Count; ip++) {
            fscanf(infile,"%lf %lf %lf %lf %lf %lf %lf\n", &pidd, &chgd, &pTd, &etad, &phid, &thetad, &zhd);
            xB.push_back(xBd);
            Q2.push_back(Q2d);
            pid.push_back(pidd);
            chg.push_back(chgd);
            eta.push_back(etad);
            phi.push_back(phid);
            theta.push_back(thetad);
            zh.push_back(zhd);
        }
        conut_particle = conut_particle + Count;
        breaks.push_back(conut_particle);
        Neves++;
    }
    fclose(infile);
    
    /*
    // Define the angular variable
    double u[41] = {-9.   , -8.675, -8.35 , -8.025, -7.7  , -7.375, -7.05 , -6.725,
       -6.4  , -6.075, -5.75 , -5.425, -5.1  , -4.775, -4.45 , -4.125,
       -3.8  , -3.475, -3.15 , -2.825, -2.5  , -2.175, -1.85 , -1.525,
       -1.2  , -0.875, -0.55 , -0.225,  0.1  ,  0.425,  0.75 ,  1.075,
        1.4  ,  1.725,  2.05 ,  2.375,  2.7  ,  3.025,  3.35 ,  3.675,
        4.};
    std::vector<double> taubins;
    for (double val : u) {
        taubins.push_back(1.0 / (1.0 + std::exp(-val)));
    }
    std::vector<double> dtau;
    for (size_t i = 0; i < taubins.size() - 1; i++) {
        dtau.push_back(taubins[i + 1] - taubins[i]);
    }
    std::vector<double> tau;
    for (size_t i = 0; i < taubins.size() - 1; i++) {
        tau.push_back((taubins[i + 1] + taubins[i]) / 2.0);
    }
    */
    std::vector<double> tau;
    std::vector<double> taubins;
    std::vector<double> dtau;
    double u_start = -15.0;
    double u_end = 4.0;
    int numPoints = 70;
    double du = (u_end - u_start) / (numPoints - 1);

    for (int i = 0; i < numPoints; i++) {
        double u = u_start + i * du;
        double taubin = 1.0 / (1.0 + exp(-u));
        taubins.push_back(taubin);
    }

    for (int i = 0; i < numPoints - 1; i++) {
        double dtaubin = taubins[i + 1] - taubins[i];
        dtau.push_back(dtaubin);
        double tau_val = (taubins[i + 1] + taubins[i]) / 2.0;
        tau.push_back(tau_val);
    }
    
    std::vector<double> Res(tau.size(), 0.0);
    
    // Look over each event
    std::cout << "Nevent " <<  breaks.size() << std::endl;
    for (size_t idx = 0; idx < breaks.size()-1; idx++) {
        std::vector<double> thetas, phis, zs;
        for (size_t i = breaks[idx]; i < breaks[idx+1]; i++) {
            thetas.push_back(theta[i]);
            phis.push_back(phi[i]);
            zs.push_back(zh[i]);
        }
        std::vector<double> tau_list;
        std::vector<double> weight_list;
        for (size_t i = 0; i < zs.size(); i++) {
            for (size_t j = i + 1; j < zs.size(); j++) {
                double thetaI = thetas[i];
                double thetaJ = thetas[j];
                double phiI = phis[i];
                double phiJ = phis[j];
                double cos_val = std::sin(thetaI) * std::sin(thetaJ) * (std::cos(phiI - phiJ) - 1)
                    + std::cos(thetaI - thetaJ);
                tau_list.push_back((1 - cos_val) / 2.0);
                weight_list.push_back(zs[i] * zs[j]);
            }
        }
        /*
        for (int i = 0; i < weight_list.size(); i++) {
            cout << weight_list[i] << " ";
        }
        cout << endl;
        */
        std::vector<double> H(taubins.size(), 0.0);
        for (size_t i = 0; i < tau_list.size(); i++) {
            for (size_t j = 0; j < taubins.size() - 1; j++) {
                if (tau_list[i] >= taubins[j] && tau_list[i] < taubins[j + 1]) {
                    H[j] += weight_list[i];
                    break;
                }
            }
        }
        for (size_t i = 0; i < tau.size(); i++) {
            Res[i] += H[i] / dtau[i];
        }
    }
    
    for (size_t i = 0; i < Res.size(); i++) {
        Res[i] /= Neves;
    }
    
    // Return tau and Res vectors or perform further operations
    
    // Print tau and Res vectors
    char output_filename1a[128];
    sprintf(output_filename1a,"results/output_%s.dat", argv[1]);
    std::ofstream output1a(output_filename1a);
    output1a <<  breaks.size() << "  " <<  breaks.size()  << std::endl;
    for (size_t i = 0; i < tau.size(); i++) {
       output1a  << tau[i] << "  " << tau[i]* Res[i] << std::endl;
    }
    output1a.close();
    
    return 0;
}

