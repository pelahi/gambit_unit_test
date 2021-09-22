#ifndef _ALLVARS_H
#define _ALLVARS_H

// some constants 
const double pi = 3.141592653589793238462643383279502884197;
const double twopi = 6.283185307179586476925286766559005768394;
const double pisq  = 9.869604401089358618834490999876151135314;
const double zeta2 = 1.644934066848226436472415166646025189219;
const double zeta3 = 1.202056903159594285399738161511449990765;
const double eulergamma = 0.577215664901532860606512090082402431042;
const double ln2   = 0.693147180559945309417232121458176568076;
const int Invalid = -1;

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <valarray>
#include <random>
#include <cassert>
#include <chrono>
#include <map>
#include <sys/stat.h>
#include <unistd.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#ifdef USEMPI
#include <mpi.h>
#endif

void GetMemUsage(std::string funcname, bool printreport);
struct TimeInfo{
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int start_line, end_line;
    double elapsed;
    std::string func;
};

std::chrono::time_point<std::chrono::high_resolution_clock> MyGetTime();
double MyElapsedTime(std::chrono::time_point<std::chrono::high_resolution_clock>  before);
void MyGetTimeStart(TimeInfo &, int line);
void MyGetTimeEnd(TimeInfo &, int line);
void ReportElapsedTime(TimeInfo &ti);

class Options
{
    public:
    int njets;
    unsigned long long seed;
    Options(){
        njets = 100;
        seed = 1234567890987654321ull;
    };
};

void GetArgs(int argc, char *argv[], Options &opt);
void Usage();
void Status(Options &opt);
#endif