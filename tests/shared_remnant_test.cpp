#include "modified_ff.hpp"
#include "ehijing_constants.hpp"
#include "spacetime_history.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
struct FixtureTable {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("electra-synthetic-qs-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    FixtureTable() {
        std::filesystem::create_directory(path);
        std::ofstream file(path/"Qs.dat");
        for(int i=0;i<31*31*31;++i) file << .0625 << '\n';
    }
    ~FixtureTable() {std::filesystem::remove_all(path);}
};
}
int main() {
    try {
        FixtureTable table;
        Modified_FF ff(0,1,2,0.,4.,-.2,table.path.string());
        Pythia8::Pythia pythia("",false);
        pythia.event.init("synthetic",&pythia.particleData);pythia.event.reset();
        auto target=std::make_shared<EHIJING::NuclearTarget>(2,1,.5);
        target->sample([](){return .25;});pythia.event.nuclearTarget(target);
        pythia.event.setRx(0.);pythia.event.setRy(0.);pythia.event.setRz(0.);
        pythia.event.setSeparationScale(.0625);
        // Soft massless jets avoid radiation and isolate exact elastic handling.
        int first=pythia.event.append(2,23,101,0,0.,0.,.35,.35,0.);
        int second=pythia.event.append(-2,23,0,101,0.,0.,.35,.35,0.);
        pythia.event[first].AddMultipleCollision({3.,1.},{.001,.001},{.3,2.},{1,1});
        pythia.event[second].AddMultipleCollision({2.},{.001},{4.},{1});
        const auto initial=pythia.event[first].p()+pythia.event[second].p()+
            Pythia8::Vec4(0.,0.,0.,ehijing::constants::mass::neutron_gev);
        const int before=pythia.event.size();
        DISKinematics kin;kin.pProton=Pythia8::Vec4(0.,0.,0.,.938);kin.Q2=2.;kin.bjorken_x=.1;
        double x=0.,y=0.,z=0.;SpacetimeHistory history(0,1);
        require(ff.sample_ff_partons(pythia,kin,x,y,z,&history),"Valid shared recoil rejected");
        require(pythia.event.size()==before+2,"Repeated hits created extra quark-diquark pairs");
        require(target->nucleons()[1].hits==3 && target->nucleons()[0].hits==1,"Hit identities lost");
        Pythia8::Vec4 total;
        for(int i=first;i<pythia.event.size();++i) total+=pythia.event[i].p();
        auto delta=total-initial;require(std::abs(delta.e())+delta.pAbs()<1e-9,"Shared recoil/split lost momentum");
        for(int i=before;i<pythia.event.size();++i) require(pythia.event[i].nuclearTargetId()==1,"Remnant ID lost");
        require(std::abs(pythia.event[before].tProd()/ehijing::constants::units::fm_to_mm/
            ehijing::constants::units::hbarc_gev_fm-3.)<1e-10,"Split vertex not at last chronological hit");
        std::ostringstream ascii;history.write(ascii);
        require(ascii.str().find("nuclear_target_id 1")!=std::string::npos,"Collision history lacks target ID");
        // Repeat with a boosted target/output frame: the rest-frame hit time
        // is transformed together with spatial displacement, once.
        pythia.event.reset();pythia.event.nuclearTarget(target);
        pythia.event.setSeparationScale(.0625);
        kin.pProton=Pythia8::Vec4(0.,0.,2.,std::sqrt(4.+.938*.938));
        Pythia8::Vec4 boosted(0.,0.,.35,.35);boosted.bst(kin.pProton);
        int boosted_jet=pythia.event.append(2,23,101,0,boosted.px(),boosted.py(),boosted.pz(),boosted.e(),0.);
        pythia.event[boosted_jet].AddMultipleCollision({2.},{.001},{.3},{1});
        require(ff.sample_ff_partons(pythia,kin,x,y,z),"Boosted exact recoil rejected");
        auto vertex=pythia.event[pythia.event.size()-1].vProd();vertex.bstback(kin.pProton);
        require(std::abs(vertex.e()/ehijing::constants::units::fm_to_mm/
            ehijing::constants::units::hbarc_gev_fm-2.)<1e-10,"Collision vertex boost/conversion inconsistent");
        // Unsatisfiable transfers return rejection, never append invented remnants.
        pythia.event.reset();pythia.event.nuclearTarget(target);pythia.event.setSeparationScale(.0625);
        kin.pProton=Pythia8::Vec4(0.,0.,0.,.938);
        int jet=pythia.event.append(2,23,101,0,0.,0.,.35,.35,0.);
        pythia.event[jet].AddMultipleCollision({1.},{100.},{.3},{1});
        require(!ff.sample_ff_partons(pythia,kin,x,y,z),"Impossible exact recoil not rejected");
        std::cout << "Shared-remnant and rejection checks passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
