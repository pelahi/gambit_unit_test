#include "allvars.hpp"
#include "jet.hpp"

using namespace std;
using namespace jets;
typedef float real;
// typedef double real;
int main() {
    TimeInfo ti;
#ifndef VECTORIZATIONTEST
    ClusterSequence cs;
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
#endif

#ifdef VECTORIZATIONTEST
    int N=1000000;
    vector<real> x, y, z;
    // real *x, *y, *z;
    x.resize(N);
    y.resize(N);
    z.resize(N);
    // x=new real[N];
    // y=new real[N];
    // z=new real[N];
    MyGetTimeStart(ti, __LINE__);
    #pragma omp simd
    for (auto i=0;i<N;i++) 
        z[i] = (x[i] * y[i])+(x[i] * y[i]);
    MyGetTimeEnd(ti, __LINE__);
    ReportElapsedTime(ti);
#endif

    return 0;
}