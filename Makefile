CPUTYPE ?= GNU
# CPUTYPE = CLANG
# GPUTYPE = AMD
GPUTYPE ?= NVIDIA
OPTLEVEL ?= 2

CC = gcc-11
CXX = g++-11
CLANG = clang
CLANGCXX = clang++
AOMP = aompcc
AOMPCXX = aompcc
AOMPACCC = /opt/sourceryg++-2020.05/bin/x86_64-none-linux-gnu-gcc
AOMPACCCXX = /opt/sourceryg++-2020.05/bin/x86_64-none-linux-gnu-g++

NCC = nvcc
NCXX = nvcc
HCC = hipcc
HCXX = hipcc
HIPIFY = hipify-perl -inplace

OMPCC = $(CC)
OMPCXX = $(CXX)
OACCCC = $(CC)
OACCCXX = $(CXX)
OCLC = $(CC)
OCLCXX = $(CXX)

ifeq ($(CPUTYPE), CLANG)
	OMPCC = $(CLANG)
	OMPCXX = $(CLANGCXX)
	OACCCC = $(CLANG)
	OACCCXX = $(CLANGCXX)
endif
ifeq ($(CPUTYPE), AOMP)
	OMPCC = $(AOMP)
	OMPCXX = $(AOMPCXX)
	OACCCC = $(AOMPACCC)
	OACCCXX = $(AOMPACCCXX)
endif
ifeq ($(CPUTYPE), CLANG-M1)
	OMPCC = $(CLANG)
	OMPCXX = $(CLANGCXX)
	OACCCC = $(CLANG)
	OACCCXX = $(CLANGCXX)
endif

ifeq ($(GPUTYPE), AMD)
	OCLC = $(HCC)
	OCLCXX = $(HCXX)
endif
ifeq ($(GPUTYPE), NVIDIA)
	OCLC = $(NCC)
	OCLCXX = $(NCXX)
endif
ifeq ($(GPUTYPE), M1)
	OCLC = $(CC)
	OCLCXX = $(CXX)
endif

OMP_FLAGS = -fopenmp -DUSEOPENMP
OMPTARGET_FLAGS = -fopenmp -DUSEOPENMPTARGET
OACC_FLAGS = -fopenacc -fopt-info-optimized-omp -DUSEOPENACC

CUDA_FLAGS = -DUSECUDA
HIPIFY_FLAGS = -DUSEHIPIFY
HIP_FLAGS = -DUSEHIP
OCL_FLAGS = -lOpenCL -DUSEOPENCL
ifeq ($(CPUTYPE), CLANG)
    OMPTARGET_FLAGS += --offload-arch=gfx906
endif
ifeq ($(CPUTYPE), AOMP)
    OMP_FLAGS = -DUSEOPENMP
    OMPTARGET_FLAGS = -DUSEOPENMPTARGET
endif
ifeq ($(CPUTYPE), CLANG-M1)
    OMP_FLAGS = -arch arm64 -Xclang -fopenmp -DUSEOPENMP -I/opt/homebrew/include/ -L/opt/homebrew/lib/
    OMPTARGET_FLAGS = -arch arm64 -Xclang -fopenmp -DUSEOPENMPTARGET -I/opt/homebrew/include/
    OACC_FLAGS = -arch arm64 -Xclang -fopenacc -fopt-info-optimized-omp -DUSEOPENACC -I/opt/homebrew/include/
endif


ifeq ($(GPUTYPE), AMD)
    OCL_FLAGS += -I/opt/rocm/opencl/include/ -L/opt/rocm/opencl/lib/
endif
ifeq ($(GPUTYPE), NVIDIA)
    OMPTARGET_FLAGS+="-foffload=-lm -fno-fast-math -fno-associative-math"
endif
ifeq ($(GPUTYPE), NVIDIA)
    OACC_FLAGS+="-foffload=-lm -fno-fast-math -fno-associative-math"
endif
ifeq ($(GPUTYPE), M1)
    OCL_FLAGS+=-I/opt/homebrew/include/ -L/opt/homebrew/lib/
endif

CFLAGS = -O$(OPTLEVEL) -D$(CALCTYPE) -D$(FUNCTYPE)
CXXFLAGS = -std=c++17 -O$(OPTLEVEL) -Iinclude/ -fopenmp

BINDIR = $(shell pwd)/bin/
SRCDIR = $(shell pwd)/src/
INCLUDE = $(SRCDIR)/*.h

all: gambit_unit_test

# explicit compilation

gambit_unit_test: obj/main.o obj/jet.o 
	$(CXX) $(CXXFLAGS) -o bin/gambit_unit_test obj/main.o obj/jet.o 


obj/main.o : src/main.cpp include/allvars.hpp include/jet.hpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o obj/main.o

obj/jet.o : src/jet.cpp include/allvars.hpp include/jet.hpp
	$(CXX) $(CXXFLAGS) -c src/jet.cpp -o obj/jet.o


.PHONY: clean all gambit_unit_test

clean:
	rm bin/gambit_unit_test
	rm obj/*.o
