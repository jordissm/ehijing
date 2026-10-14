#include "eHIJING/NuclearTarget.h"
#include "eHIJING/eHIJING.h"
#include "nuclear_recoil.hpp"
#include "Pythia8/Event.h"
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
void close(double a,double b,double tolerance,const char* message) {
    require(std::abs(a-b) < tolerance,message);
}
void elastic(Pythia8::Vec4& p,Pythia8::Vec4& k,double q2,double phi) {
    const auto sum = p+k; const double pm=p.m2Calc(),km=k.m2Calc();
    require(ehijing::recoil_two_body(p,k,q2,phi),"Physical collision rejected");
    const auto difference=p+k-sum;
    require(std::abs(difference.e())+difference.pAbs() < 1e-10,"Four-momentum not conserved");
    close(p.m2Calc(),pm,1e-9,"Jet mass changed");
    close(k.m2Calc(),km,1e-9,"Target mass changed");
}
}
int main() {
    try {
        std::mt19937 engine(17); std::uniform_real_distribution<double> uniform(0.,1.);
        auto random = [&]() { return uniform(engine); };
        EHIJING::NuclearTarget proton(1,1,.5); proton.sample(random);
        require(proton.primary()==0 && proton.nucleons()[0].hits==1,"Single-nucleon primary");
        require(proton.thickness({0.,0.,0.},{0.,0.,1.})==0.,"Primary included in medium");
        EHIJING::MultipleCollision no_medium(0.,4.,-.2);
        std::vector<double> q2,times,phis;std::vector<int> ids;
        require(no_medium.sample_all_qt2_target(2,10.,proton,{0.,0.,0.},{0.,0.,1.},.1,2.,q2,times,phis,ids)==0,
                "Single nucleon unexpectedly scatters on primary");
        EHIJING::table_nd cube(3,{2,2,2},{0.,0.,0.},{1.,1.,1.});
        for(int i=0;i<8;++i) cube.set_with_linear_index(i,7.);
        close(cube.interpolate({1.,1.,1.}),7.,1e-14,"Unsafe upper interpolation endpoint");
        EHIJING::NuclearTarget deuteron(2,1,.5); deuteron.sample(random);
        const auto& d = deuteron.nucleons();
        for (int j=0;j<3;++j) close(d[0].position[j],-d[1].position[j],1e-12,"Hulthen COM");
        const int spectator=1-deuteron.primary();
        deuteron.nucleons()[spectator].position={0.,0.,0.};
        const double w=deuteron.width();
        const double expected=.5/(2.*3.141592653589793*w*w);
        close(deuteron.thickness({0.,0.,0.},{0.,0.,1.}),expected,1e-12,"Gaussian half-ray integral");
        // Independent numerical integration of the three-dimensional profile.
        double integral=0.; const int intervals=20000; const double step=12.*w/intervals;
        for(int i=0;i<=intervals;++i) integral += (i==0 || i==intervals ? .5 : 1.)*
            deuteron.density({0.,0.,i*step})*step;
        close(integral,expected,1e-9,"Density/path normalization mismatch");
        int cdf_count=0;
        for(int i=0;i<20000;++i) {
            double time; require(deuteron.sample_collision({0.,0.,0.},{0.,0.,.5},random,time)==spectator,
                "Collision selected primary");
            if(time*.5<w) ++cdf_count;
        }
        close(cdf_count/20000.,std::erf(1./std::sqrt(2.)),.015,"Joint time/profile sample");
        EHIJING::NuclearTarget lead(208,82,.5); int selected_protons=0;
        for(int i=0;i<2000;++i) {
            lead.sample(random); int protons=0,hits=0;
            for(const auto& n:lead.nucleons()) { protons+=n.proton; hits+=n.hits; }
            require(protons==82 && hits==1,"Target composition/primary not exact");
            selected_protons+=lead.nucleons()[lead.primary()].proton;
        }
        close(selected_protons/2000.,82./208.,.04,"Z/A primary selection");
        bool invalid=false; try { EHIJING::NuclearTarget unsupported(12,6,.5); }
        catch(const std::invalid_argument&) {invalid=true;}
        require(invalid,"Unsupported Woods-Saxon parameters silently invented");
        EHIJING::NuclearTarget custom(12,6,.5,2.4,.5); custom.sample(random);
        Pythia8::Vec4 p(0.,0.,10.,10.),k(0.,0.,0.,.938);
        elastic(p,k,.1,.3); elastic(p,k,.1,2.); elastic(p,k,.03,4.);
        p=Pythia8::Vec4(1.,2.,3.,std::sqrt(14.+.33*.33)); k=Pythia8::Vec4(.1,-.2,.3,std::sqrt(.14+.938*.938));
        elastic(p,k,.01,.8);
        p=Pythia8::Vec4(0.,0.,1.,1.);k=Pythia8::Vec4(0.,0.,0.,.938);auto before=p;
        require(!ehijing::recoil_two_body(p,k,100.,0.),"Impossible transfer accepted");
        close(p.e(),before.e(),1e-14,"Rejected transfer mutated jet");
        Pythia8::Particle original;original.AddMultipleCollision({2.},{.1},{.4},{7});original.nuclearTargetId(7);
        Pythia8::Particle copied=original,assigned;assigned=original;
        require(copied.coll_targets()==original.coll_targets() && assigned.nuclearTargetId()==7,"Identity lost on copy");
        copied.ResetMultipleCollision();require(copied.coll_targets().empty(),"Stale collision identity after reset");
        Pythia8::Event event;event.nuclearTarget(std::make_shared<EHIJING::NuclearTarget>(2,1,.5));
        Pythia8::Event saved=event;require(saved.nuclearTarget()==event.nuclearTarget(),"Shared target lost on event copy");
        event.clear();require(!event.nuclearTarget(),"Previous target leaked into next event");
        std::cout << "Geometry, identity and exact recoil checks passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
