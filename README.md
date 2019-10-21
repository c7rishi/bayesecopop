# bayesecopop

## Installation

In Windows: add the following lines to Documents/.R/Makevars.win 
for compilation:

```
CXX14=$(BINPREF)g++ -O2 -march=native -mtune=native -m$(WIN) -std=c++1y
CXX14FLAGS=-O3 -march=native -mtune=native
CXX11FLAGS=-O3 -march=native -mtune=native
PKG_LIBS = $(SHLIB_OPENMP_CFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

Then run
```
devtools::install_github(
  "c7rishi/bayesecopop",
  auth_token = "2e3ff8f12251e40ab2c6cab4942745ef0d52363d"
)
```
