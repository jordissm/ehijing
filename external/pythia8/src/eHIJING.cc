#include "eHIJING/eHIJING.h"
#include <iostream>
#include <fstream>
#include "eHIJING/integrator.h"
#include <cmath>
#include <thread>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <string>
#include <sstream>
#include <filesystem>
#include <atomic>

namespace EHIJING {

// Unit conversion constants
// https://pdg.lbl.gov/2025/reviews/rpp2025-rev-phys-constants.pdf
const double GeVfm = 1. / 0.197'326'980'4; // ≈ 5.068
const double GeV2fm2 = GeVfm * GeVfm;
const double GeV3fm3 = GeVfm * GeV2fm2;

// Color algebra constants
const double CA = 3;        // Casimir for gluon emission from a gluon
const double CF = 4./3.;    // Casimir for gluon emission from a quark
const double dA = 8;        // Dimension of the adjoint representation, i.e., N_c^2 - 1 for SU(N_c)

// Lambda_{QCD} and Lambda_{QCD}^2
const double mu = 0.25;             // [GeV]
const double mu2 = std::pow(mu, 2); // [GeV]^2

// Minimum and maximum values of nuclear thickness function TA
const double TAmin = 0.1 / GeV2fm2; // [GeV]^2
const double TAmax = 2.8 / GeV2fm2; // [GeV]^2

// 1-loop QCD beta function coefficient
// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-qcd.pdf#page=2.63
// b0 = (11 * N_c - 2 * n_f) / (6 * 2 * pi)
//    = (9 / 2) / (2 * pi) for N_c = 3 and n_f = 3
const double b0 = 9. / 2.;

// Handy constants
const double twoPioverb0 = 2. * M_PI / b0;
const double piCAoverdA = M_PI * CA / dA;
const double Mproton = 0.938;                   // [GeV]
const double rho0 = 0.17 / GeV3fm3;             //
const double r0 = 1.12 * GeVfm;                 // [GeV]^{-1}

// Strong coupling constant
double alphas(double Q2) {
    double scale_ratio_squared = Q2 / mu2;
    if (scale_ratio_squared < 2.71828) {
        return twoPioverb0;
    } else {
        return twoPioverb0 / std::log(scale_ratio_squared);
    }
}

// Gluon TMD that reduces to the Weizsäcker-Williams distribution at large momentum-squared k_T^2
// Eq. (19) in https://doi.org/10.1103/PhysRevD.110.034001
double alphas_PhiG(double x, double kT2, double Qs2,
		   double powerG, double lambdaG, double K){
    return K * std::pow(1. - x, powerG) * std::pow(x, lambdaG) / (kT2 + Qs2);
}

// integrate du (1-cos(1/u)) / u
double inte_C(double u){
    return gsl_sf_Ci(1./u) + std::log(u);
}

double CHT_F1(double x){
    return 1.0 - sin(x)/x;
}

///////// Class: Multiple Collision /////////////////////
// initializer
MultipleCollision::MultipleCollision(double Kfactor, double powerG, double lambdaG):
Kfactor_(Kfactor),
powerG_(powerG),
lambdaG_(lambdaG),
rd(),
gen(rd()),
flat_gen(0.,1.),
Qs2Table(3, {31,31,31}, // TA, ln(x), ln(Q2) --> size and grid for Qs table
           {TAmin, std::log(1e-6), std::log(1.0)},
           {TAmax, std::log(.99), std::log(1e5)}
       )
{
}
// Tabulate Qs
void MultipleCollision::Tabulate(std::filesystem::path table_path){
    std::filesystem::path fname = table_path/std::filesystem::path("Qs.dat");
    if (std::filesystem::exists(fname)) {
        std::cout << "Loading Qs Table" << std::endl;
        std::ifstream f(fname.c_str());
        int count = 0;
        double entry;
        std::string line;
        while(getline(f, line)){
            std::istringstream in(line);
            in >> entry;
            if (count >= Qs2Table.size()){
                std::cerr << "Loading table Qs: mismatched size - 1" << std::endl;
                exit(-1);
            }
            Qs2Table.set_with_linear_index(count, entry);
	    count ++;
        }
        if (count<Qs2Table.size()){
            std::cerr << "Loading table Qs: mismatched size - 2" << std::endl;
            exit(-1);
        }
    } else {
        try {
            const std::filesystem::path out_path(fname);
            const std::filesystem::path dir_path = out_path.parent_path();

            // Create parent directory tree if needed
            if (!dir_path.empty()) {
                std::error_code ec;
                std::filesystem::create_directories(dir_path, ec);
                if (ec) {
                    throw std::runtime_error(
                        "Failed to create directory '" + dir_path.string() +
                        "': " + ec.message());
                }
            }

            // Write to temporary file first
            const std::filesystem::path tmp_path = out_path.string() + ".tmp";

            std::ofstream f(tmp_path, std::ios::out | std::ios::trunc);
            if (!f.is_open()) {
                throw std::runtime_error(
                    "Failed to open temporary output file '" + tmp_path.string() + "'");
            }

            std::cout << "Generating Qs^2 table: " << out_path << std::endl;

            // Table Qs as a function of lnx, lnQ2, TA
            for (int c = 0; c < Qs2Table.size(); ++c) {

                // 
                auto index = Qs2Table.LinearIndex2ArrayIndex(c);
                auto xvals = Qs2Table.ArrayIndex2Xvalues(index);

                const double TA = xvals[0];
                const double xB = std::exp(xvals[1]);
                const double Q2 = std::exp(xvals[2]);

                const double aQs2 = compute_Qs2(TA, xB, Q2);

                // Validate value before storing/writing
                if (!std::isfinite(aQs2) || aQs2 < 0.0) {
                    throw std::runtime_error(
                        "Invalid Qs^2 value at linear index " + std::to_string(c) +
                        " (TA=" + std::to_string(TA) +
                        ", xB=" + std::to_string(xB) +
                        ", Q2=" + std::to_string(Q2) +
                        ", Qs2=" + std::to_string(aQs2) + ")");
                }

                Qs2Table.set_with_linear_index(c, aQs2);

                f << aQs2 << '\n';
                if (!f) {
                    throw std::runtime_error(
                        "Write failure while writing temporary file '" + tmp_path.string() + "'");
                }
            }

            f.flush();
            if (!f) {
                throw std::runtime_error(
                    "Flush failure for temporary file '" + tmp_path.string() + "'");
            }

            f.close();
            if (!f) {
                throw std::runtime_error(
                    "Close failure for temporary file '" + tmp_path.string() + "'");
            }

            // Atomically replace final file
            std::error_code ec;
            std::filesystem::rename(tmp_path, out_path, ec);
            if (ec) {
                // On some systems rename may fail if target exists
                std::filesystem::remove(out_path, ec);
                ec.clear();
                std::filesystem::rename(tmp_path, out_path, ec);
                if (ec) {
                    throw std::runtime_error(
                        "Failed to move temporary file '" + tmp_path.string() +
                        "' to final file '" + out_path.string() +
                        "': " + ec.message());
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error generating Qs^2 table '" << fname
                    << "': " << e.what() << std::endl;
            throw;
        }
    }
}

/*
 * Self-consistent equation for Qs^2
 * Eq. (20) in https://doi.org/10.1103/PhysRevD.110.034001
 * Before solving, a change of variable is performed to obtain
 * a more well-behaved integrand for numerical integration.
 * The change of variable is
 * 
 *                  u = ln(1 + q^2 / Q_s^2),
 * 
 * which implies
 * 
 *           Q_s^2 = (pi * C_A * T_A / d_A) 
 *                   * int du alpha_s 
 *                   * phi_g(x_g(u), q^2(u), Q_s^2) 
 *                   * (Q_s^2 + q^2(u)).
 */
double MultipleCollision::Qs2_self_consistent_eq(double Qs2, double TA, double xB, double Q2) {

    // Prefactor constants
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
}

/*
 * Compute Qs^2 using the bisection method to solve its self-consistent equation.
 */
double MultipleCollision::compute_Qs2(double TA, double xB, double Q2){
    const double EPS = 1e-4;
    const int MAX_BRACKET_ITERS = 100;
    // double xleft = mu2, xright = Q2*100;

    // Starting points for bracketing
    double xleft  = mu2;
    double xright = std::max(Q2, mu2) * 100.0;

    double yleft  = Qs2_self_consistent_eq(xleft,  TA, xB, Q2);
    double yright = Qs2_self_consistent_eq(xright, TA, xB, Q2);

    // If already exactly on the root
    if (std::abs(yleft) < EPS) {
        return xleft;
    }

    int iter = 0;
    while (yleft * yright > 0.0 && iter < MAX_BRACKET_ITERS) {
        xleft *= 0.5;
        xright *= 2.0;
        yleft  = Qs2_self_consistent_eq(xleft, TA, xB, Q2);
        yright = Qs2_self_consistent_eq(xright, TA, xB, Q2);
        ++iter;
    }

    // If still not bracketed, give up safely
    if (yleft * yright > 0.0) {
        std::cerr << "eHIJING warning: could not bracket Qs2 root; "
                  << "TA=" << TA
                  << ", xB=" << xB
                  << ", Q2=" << Q2
                  << ", yleft=" << yleft
                  << ", yright=" << yright
                  << ", using Qs2 = mu2"
                  << std::endl;
        return mu2;
    }

    // Standard bisection
    while ((xright - xleft) > EPS) {
        const double xmid = 0.5 * (xleft + xright);
        const double ymid = Qs2_self_consistent_eq(xmid, TA, xB, Q2);

        if (std::abs(ymid) < EPS) {
            return xmid;
        }

        if (yleft * ymid < 0.0) {
            xright = xmid;
            yright = ymid;
        } else {
            xleft = xmid;
            yleft = ymid;
        }
    }

    return 0.5 * (xleft + xright);
}

// Sample all elastic collisions, without radiation, ordered from high to low scales
int MultipleCollision::sample_all_qt2(int pid, double E, double L, double Thickness, double xB, double Q2,
                             std::vector<double> & q2_list, std::vector<double> & t_list,
                             std::vector<double> & phi_list) {
    q2_list.clear();
    t_list.clear();
    phi_list.clear();
    double TA = rho0 * L;
    double CR = (pid==21)? CA : CF;
    double qs2 = TA / Thickness * Qs2(Thickness, xB, Q2);
    double tildeTA = Kfactor_ * M_PI * CR * TA / dA;
    if (tildeTA < 1e-6) return q2_list.size();
    double qs2overCTA = qs2 / tildeTA;
    double q2max = Q2 / xB;
    double q2min = 0.01 * EHIJING::mu2;
    double q2 = q2max;
    double t0 = EHIJING::r0; // exclude double scattering on the same nucleon
    if (L < t0) return q2_list.size();
    while (q2 > q2min) {
        // Sample the next hard multiple collision
        double lnr = std::log(flat_gen(gen));
        double xg = q2 / q2max;
	    q2 = q2 * std::pow(1.0 + (lambdaG_ - 1.0) * lnr * q2 / tildeTA * std::pow(xg, -lambdaG_), 1.0 / (lambdaG_ - 1.0));
        double t = t0 + flat_gen(gen) * (L - t0);
        if ( flat_gen(gen) < std::pow(1.0 - xg, powerG_) * (q2 / (q2+qs2)) ) {
            // Correct for the gluon distribtuion at large x_g and the screening effect
            q2_list.push_back(q2);
            t_list.push_back(t);
            phi_list.push_back(2.0 * M_PI * flat_gen(gen));
        }
    }
    return q2_list.size();
}

/////////// Class: eHIJING
// initializer
eHIJING::eHIJING(int mode, double Kfactor,
                 double powerG, double lambdaG):
MultipleCollision(Kfactor, powerG, lambdaG),
mode_(mode),
Kfactor_(Kfactor),
powerG_(powerG),
lambdaG_(lambdaG),
rd(),
gen(rd()),
flat_gen(0.,1.),
GHT_Angular_Table(2, {51, 51}, // X = delta = 2kq/(k^2+q^2),
                               // ln(1+Y) = ln(1+ (|k|-|q|)^2*t/(2*z*(1-z)*E) )
                               // 0.5<delta<1, ln(1)<ln(1+Y)<ln(11)
           {0.5, 1e-3},
           {0.99, std::log(11.0)}
       )
{
}

// Tabulate the Qs and (if necessary) the GHT / GLV table (collinear H-T do not need a separate table)
void eHIJING::Tabulate(std::filesystem::path table_path) {
    MultipleCollision::Tabulate(table_path);
    // If mode = 1, use GHT
    if (mode_ == 1) {
        // GHT table is at least 4-dimensional with each entry a 2D integral
        // the following routine will use the max number of hardware concurrency of your computer
        // to parallel the computation of the table
        std::filesystem::path fname = table_path/std::filesystem::path("GHT.dat");
        if (std::filesystem::exists(fname)) {
            std::cout << "Loading GHT Table" << std::endl;
            std::ifstream f(fname.c_str());
            int count = 0;
            double entry1, entry2;
            std::string line;
            while(getline(f, line)){
                std::istringstream in(line);
                in >> entry1;
                if (count>=GHT_Angular_Table.size()){
                    std::cerr << "Loading table GHT: mismatched size - 1" << std::endl;
                    exit(-1);
                }
                GHT_Angular_Table.set_with_linear_index(count, entry1);
                count ++;
            }
            if (count<GHT_Angular_Table.size()){
                std::cerr << "Loading table GHT: mismatched size - 2" << std::endl;
                exit(-1);
            }
        } else {
            std::filesystem::create_directory(table_path);
            std::ofstream f(fname.c_str());
            // Table Qs as a function of lnx, lnQ2, TA
            std::atomic_int counter =  0;
            int percentbatch = int(GHT_Angular_Table.size()/100.);
            auto code = [this, percentbatch](int start, int end) {
                static std::atomic_int counter;
                for (int c=start; c<end; c++) {
                    counter ++;

                    if (counter%percentbatch==0) {
                      std::cout <<std::flush << "\r" << counter/percentbatch << "% done";
                    }
                    auto index = GHT_Angular_Table.LinearIndex2ArrayIndex(c);
                    auto xvals = GHT_Angular_Table.ArrayIndex2Xvalues(index);
                    double X = xvals[0];
                    double ln1Y = xvals[1];
                    double Y = std::exp(ln1Y)-1.;
                    double entry1 = 0.;
                    entry1 = compute_GHT_Angular_Table(X, Y);
                    GHT_Angular_Table.set_with_linear_index(c, entry1);
                }
            };
            std::vector<std::thread> threads;
            int nthreads = std::thread::hardware_concurrency();
            int padding = int(std::ceil(GHT_Angular_Table.size()*1./nthreads));
            std::cout << "Generating GHT angular tables with " << nthreads << " thread" << std::endl;
            for(auto i=0; i<nthreads; ++i) {
                int start = i*padding;
                int end = std::min(padding*(i+1), GHT_Angular_Table.size());
                threads.push_back( std::thread(code, start, end) );
            }
            for(auto& t : threads) t.join();
            for (int c=0; c<GHT_Angular_Table.size(); c++) {
                f << GHT_Angular_Table.get_with_linear_index(c) << std::endl;
            }
        }
        std::cout << "... done" << std::endl;
    }
}

/*
 * Compute the GHT angular table
 */
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

bool eHIJING::next_kt2_stochastic(double & kt2, 
	                int pid,
                        double E,
                        double kt2min,
                        std::vector<double> qt2s,
                        std::vector<double> ts) {
    double CR = (pid==21)? CA : CF;
    double CAoverCR = CA/CR;
    double CR_2overb0 = CR*2.0/b0;
    double zmin = std::min(.4/E, .4);
    double zmax = 1. - zmin;
    double logvac = std::log(zmax/zmin);
    int Ncolls = ts.size();
    double acceptance = 0.;
    // If mode = 0, use HT
    if (mode_ == 0){
        while (acceptance<flat_gen(gen) && kt2>kt2min) {
            double maxlogmed = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                if (q2>kt2) continue;
                double phasemax = inte_C((2*zmax*E)/(t*kt2min)) - inte_C((2*zmin*E)/(t*kt2));
                maxlogmed += q2 * phasemax ;
            }
            maxlogmed *= 2. / kt2min * CAoverCR;
            double Crad = CR_2overb0 * (logvac + maxlogmed);
            double r = flat_gen(gen);
            kt2 = mu2 * std::pow(kt2/mu2, std::pow(r, 1.0/Crad) );
            double logmed = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                if (q2>kt2) continue;
                double phase = inte_C((2*zmax*E)/(t*kt2)) - inte_C((2*zmin*E)/(t*kt2));
                logmed += q2 * phase;
            }
            logmed *= 2./kt2*CAoverCR;
            acceptance = (logvac + logmed) / (logvac + maxlogmed);
        }
    // If mode = 1, use GHT
    } else {
        double maxdiffz = 1./zmin - 1./zmax + 2.*logvac;
        while (acceptance<flat_gen(gen) && kt2>kt2min) {
            double maxmedcoeff = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                maxmedcoeff += t*(kt2 + q2);
            }
            maxmedcoeff *= CAoverCR/(2.*E)*maxdiffz;
            double Crad = CR_2overb0 * (logvac + maxmedcoeff);
            double r = flat_gen(gen);
            kt2 = mu2 * std::pow(kt2/mu2, std::pow(r, 1.0/Crad) );
            // compute acceptance
            double medcoeff = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                medcoeff += t*(kt2 + q2);
            }
            medcoeff *= CAoverCR/(2.*E)*maxdiffz;
            acceptance = (logvac + medcoeff) / (logvac + maxmedcoeff);
        }
    }
    return (kt2>kt2min);
}

double eHIJING::sample_z_stochastic(double & z, int pid,
                        double E,
                        double kt2,
                        std::vector<double> qt2s,
                        std::vector<double> ts,
                        std::vector<double> phis) {
    double CR = (pid==21)? CA : CF;
    double zmin = std::min(.4/E, .4);
    double zmax = 1. - zmin;

    int Ncolls = ts.size();
    double acceptance = 0.;
    double weight = 1.0;
    // If mode = 0, use HT
    if (mode_ == 0) {
        while (acceptance<flat_gen(gen)) {
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
        }
    // If mode = 1, use GHT
    } else {
        double a = 0.;
        for (int i=0; i<Ncolls; i++) {
            double q2 = qt2s[i], t = ts[i];
            a += (kt2 +q2)*t;
        }
        a *= CA/CR/(2.*E);
        double Norm = (1. + 2.*a)*std::log(zmax/zmin) + a*(1./zmin-1/zmax);
        // sample z ~ 1/z + a/[z^2(1-z)]
      	double left = zmin;
      	double right = zmax;
      	double mid = (left+right)/2.;
      	double r = flat_gen(gen);
      	while(right-left > 1e-3) {
        		mid = (left+right)/2.;
        		double fmid = ( (1+a)*std::log(mid/zmin)
        			    + a*(1./zmin-1/mid)
        			    + a*std::log(zmax/(1-mid)) ) / Norm;
        		if (fmid<r) left = mid;
        		else right = mid;
      	}
        z = mid;
        if (Ncolls==0) weight=1.;
        else {
           double wmax = CR/CA;
           for (int i=0; i<Ncolls; i++) {
               double q2 = qt2s[i], t = ts[i];
               wmax += (kt2+q2)*t;
           }
           wmax /= (2.*z*(1-z)*E);

           double w = CR/CA;
           for (int i=0; i<Ncolls; i++) {
               double dw;
               double q2 = qt2s[i], t = ts[i], phi = phis[i];
	       double A = (kt2+q2)*t/(2.*z*(1-z)*E),
                      B = 2.*std::sqrt(kt2*q2)*t/(2.*z*(1-z)*E);
               double X = B/A;
               double Y = A*(1.-X);
               if (X<.5){
                    double jv0 = gsl_sf_bessel_J0(B),
                           jv1 = gsl_sf_bessel_J1(B);
                    dw = - X*std::sin(A)*jv1
                         + X*X * ( .5 + X*std::cos(A) * (jv1/B - jv0) );
               } else {
                   if (Y>10.){
                       dw = 1./std::sqrt(1.-X*X)-1;
                   }
                   else {
                       dw = GHT_Angular_Table.interpolate({X, std::log(1.+Y)});
                   }
               }
               w += dw;
           }
           weight = w/wmax;
        }
    }
    return std::min(std::max(weight,0.),1.);
}
} // namespace EHIJING
