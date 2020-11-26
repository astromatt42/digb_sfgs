INCLFLAGS = -I$(HOME)/Software/gsl-2.6/interpolation -I$(HOME)/Software/gsl-2.6/integration -I$(HOME)/Software/gsl-2.6/specfunc\
            -I$(HOME)/Software/cubature -I$(HOME)/Software/gsl_lib/lib -L$(HOME)/Software/gsl_lib/lib -L$(HOME)/Software/cubature

LDFLAGS = -Wl,-rpath,$(HOME)/Software/gsl_lib/lib:$(HOME)/Software/cubature

OPTFLAGS = -O3

CFLAGSTC = -I/home/mar/Software/gsl-2.6/integration

GDBFLAGS = -ggdb -Wall -pg

CPYFLAGS = -fpic -shared

OMPFLAGS = -fopenmp

spectra:
	gcc $(INCLFLAGS) $(OMPFLAGS) $(LDFLAGS) $(OPTFLAGS) spectra.c spectra_funcs.h -lgsl -lgslcblas -lhcubature -lm -o spectra

spectra_gas:
	gcc $(INCLFLAGS) $(OMPFLAGS) $(LDFLAGS) $(OPTFLAGS) spectra_gas.c spectra_funcs.h -lgsl -lgslcblas -lhcubature -lm -o spectra

spectra_so:
	gcc $(INCLFLAGS) $(CPYFLAGS) spectra.c spectra_funcs.h -lgsl -lgslcblas -lhcubature -lm -o spectra.so

cosmo_so:
	gcc $(INCLFLAGS) $(CPYFLAGS) cosmo_funcs.c -lgsl -lgslcblas -lm -o cosmo.so

clean:
	rm -f spectra spectra.so cosmo.so


