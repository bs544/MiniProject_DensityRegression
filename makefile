bispectrum.o : bispectrum.f90 util.o basisFns.o neighbours.o
	gfortran -c bispectrum.f90 util.o basisFns.o neighbours.o -lblas -llapack


test : unitTests.f90 util.o basisFns.o neighbours.o
	gfortran unitTests.f90 util.o basisFns.o neighbours.o -o test -lblas -llapack

neighbours.o : neighbours.f90 util.o
	gfortran -c neighbours.f90 util.o -lblas -llapack

basisFns.o : basis.f90 util.o
	gfortran -c basis.f90 util.o -o basisFns.o

util.o : util.f90
	gfortran -c util.f90 -lblas -llapack
