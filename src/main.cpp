#include "allvars.hpp"
#include "jet.hpp"

using namespace std;
using namespace jets;
int main() {
    ClusterSequence cs;
    TimeInfo ti;
    ti.func = __func__;
    MyGetTimeStart(ti, __LINE__);
    cs._jets.resize(100);
    std::default_random_engine generator;
    std::normal_distribution<double> normaldistrib(0.0,1.0);
    std::uniform_real_distribution<double> uniformdistrib(0.0,1.0);
    for (auto &j:cs._jets) {
        auto px = normaldistrib(generator);
        auto py = normaldistrib(generator);
        auto pz = normaldistrib(generator);
        auto mass = uniformdistrib(generator);
        auto E = sqrt(px*px+py*py+pz*pz+mass*mass);
        j = PseudoJet(px,py,pz,E);
    }
    cs.fill_initial_history();
    MyGetTimeEnd(ti, __LINE__);
    ReportElapsedTime(ti);
    cs.simple_N2_cluster<BriefJet>();
    return 0;
    
}