# bayesecopop

In Windows: add the following lines to Documents/.R/Makevars.win 
for compilation:

CXX14=$(BINPREF)g++ -O2 -march=native -mtune=native -m$(WIN) -std=c++1y
CXX14FLAGS=-O3 -march=native -mtune=native
CXX11FLAGS=-O3 -march=native -mtune=native
PKG_LIBS = $(SHLIB_OPENMP_CFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
