Settings for `Rcpp` and `MKL`.

Include the following in `~/.my.bashrch`:
```
if [ "${PS1:+set}" = set ]; then

	use UGER
	use OpenblasR
	use R-3.3
	use GCC-5.2
	use Python-3.4

	use .openblas-0.2.8 
	use .boost-1.60.0	
	use .zlib-1.2.8
	use .git-2.5.0
	use .openssl-1.0.2g 
	use .libssh2-1.3.0 

	use GSL
	use .icc-2015
	use .setenv-ldflags++
	use .htslib-1.3.2

	use .parallel-20140722
	use .automake-1.12.5
	use .autoconf-2.69
	use .imagemagick-6.9.2-3 

	export MKROOT="/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/intel/mkl"
	export MKLFLAGS="-m64 -I${MKLROOT}/include"
	export MKLLINK="-L${MKLROOT}/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl"

	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HOME/work/common/lib
	export LD_PRELOAD=${MKLROOT}/lib/intel64/libmkl_core.so:${MKLROOT}/lib/intel64/libmkl_sequential.so
	export LIBRARY_PATH=${LIBRARY_PATH}:$HOME/work/common/lib
	export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:$HOME/work/common/include

	export MKL="$MKLLINK $MKLFLAGS"
	export CPPFLAGS="${CPPFLAGS} -std=c++14 -O3 ${MKLFLAGS} -msse2 -DEIGEN_USE_MKL_ALL"
	export LDFLAGS="${LDFLAGS} ${MKLLINK} -lopenblas -Wl,--no-as-needed"
fi
```

Include the following in `~/.R/Makevars`:
```
CC = gcc
CXX = g++

MKROOT = /broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/intel/mkl
MKLFLAGS = -m64 -I$(MKLROOT)/include
MKLLINK = -L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl

CFLAGS = -O3 -I. -I.. -I../inst -I/broad/software/free/Linux/redhat_6_x86_64/pkgs/boost_1.60.0/include/
CXXFLAGS = -O3 -I. -I.. -I../inst -I/broad/software/free/Linux/redhat_6_x86_64/pkgs/boost_1.60.0/include/
LDFLAGS = -L.

PKG_CFLAGS = -O3 $(CFLAGS)
PKG_CXXFLAGS = $(CXXFLAGS) -std=c++14 $(MKLFLAGS) -msse2 -DEIGEN_USE_MKL_ALL
PKG_LDFLAGS = $(LDFLAGS) $(MKLFLAGS) $(MKLLINK) $(LDFLAGS)
PKG_LIBS = $(PKG_LDFLAGS)
```
