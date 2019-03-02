test : unitTests.f90 util.o basis.o neighbours.o bispectrum.o
	gfortran unitTests.f90 util.o basis.o neighbours.o bispectrum.o -o test -Wall -Wno-tabs -fPIC -lblas -llapack -O2


bispectrum.o : bispectrum.f90 util.o basis.o neighbours.o
	gfortran -c bispectrum.f90 util.o basis.o neighbours.o -fPIC -Wall -Wno-tabs -lblas -llapack

neighbours.o : neighbours.f90 util.o
	gfortran -c neighbours.f90 util.o -fPIC -lblas -llapack -Wall -Wno-tabs -O2

basis.o : basis.f90 util.o
	gfortran -c basis.f90 util.o -o basis.o -fPIC -Wall -Wno-tabs -O2

util.o : util.f90
	gfortran -c util.f90 -fPIC -lblas -llapack -Wall -Wno-tabs -O2
