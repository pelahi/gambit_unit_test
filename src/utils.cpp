#include "allvars.hpp"

/// get the memory use looking at the task
void GetMemUsage(std::string funcname, bool printreport){
#ifndef USEMPI
    int ThisTask=0;
#endif
    std::ifstream f;
    size_t pos1, pos2;
    unsigned long long size, resident, shared, text, data, library, dirty, peak;
    std::map<std::string,float> memuse;
    bool iflag = true;
    char buffer[2500];
    std::string memreport, temp, delimiter1("VmPeak:"), delimiter2("kB");
    // Open the file storing memory information associated with the process;
    f.open("/proc/self/statm");
    if (f.is_open()) {
        f >> size >> resident >> shared >> text >> library >> data >>dirty;
        // nscan = fscanf(file, "%lld %lld %lld %lld %lld %lld %lld",
        //     &size, &resident, &shared, &text, &library, &data, &dirty);
        f.close();
    }
    else {
        iflag = false;
    }
    f.open("/proc/self/status");
    if (f.is_open()) {
        while (f.getline(buffer, 2500)){
            temp = std::string(buffer);
            if ((pos1 = temp.find(delimiter1)) != std::string::npos) {
                pos2 = temp.find(delimiter2);
                temp = temp.substr(pos1+delimiter1.size(), pos2);
                peak = stol(temp)*1024;
                break;
            }
        }
        f.close();
    }
    else {
        iflag = false;
    }

    memreport += std::string("Memory report, func = ")+funcname+std::string(" task = ")+std::to_string(ThisTask)+std::string(" : ");

    //having scanned data for memory footprint in pages, report
    //memory footprint in GB
    if (iflag) {
        // Convert pages into bytes. Usually 4096, but could be 512 on some
        // systems so take care in conversion to KB. */
        uint64_t sz = sysconf(_SC_PAGESIZE);
        size *= sz ;
        resident *= sz ;
        shared *= sz ;
        text *= sz ;
        library *= sz ;
        data *=  sz ;
        dirty *= sz ;
        float bytestoGB;
        bytestoGB = 1.0/(1024.0*1024.*1024.);
        // float bytestoGB = 1.0;
        // if (opt.memuse_peak < peak) opt.memuse_peak = peak;
        // opt.memuse_nsamples++;
        // opt.memuse_ave += size;
        memuse["Size"] = size*bytestoGB;
        memuse["Resident"] = resident*bytestoGB;
        memuse["Shared"] = shared*bytestoGB;
        memuse["Text"] = text*bytestoGB;
        memuse["Data"] = data*bytestoGB;
        memuse["Peak"] = peak*bytestoGB;
        for (std::map<std::string,float>::iterator it=memuse.begin(); it!=memuse.end(); ++it) {
            memreport += it->first + std::string(" = ") + std::to_string(it->second) + std::string(" GB, ");
        }
    }
    else{
        memreport+= std::string(" unable to open or scane system file storing memory use");
    }
    if (printreport) std::cout<<memreport<<std::endl;
}

std::chrono::time_point<std::chrono::high_resolution_clock> MyGetTime(){
    auto now = std::chrono::high_resolution_clock::now();
    return now;
}
double MyElapsedTime(std::chrono::time_point<std::chrono::high_resolution_clock> before)
{
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(now - before);
    return elapsed.count()*1e-9;
}
void MyGetTimeStart(TimeInfo &ti, int line){
    ti.start = MyGetTime();
    ti.start_line = line;
}
void MyGetTimeEnd(TimeInfo &ti, int line){
    ti.end = MyGetTime();
    ti.end_line = line;
    ti.elapsed = MyElapsedTime(ti.start);
}

void ReportElapsedTime(TimeInfo &ti) {
    std::cout<<" Time taken for "<<ti.func<<" from lines "<<ti.start_line<<" - "<<ti.end_line<<" : "<<ti.elapsed<<std::endl;
}

