
test : unitTests.f90 util.o basisFns.o
	gfortran unitTests.f90 util.o basisFns.o -o test -lblas

#basisFns.o : sphericalharm.f90 util.o
#	gfortran -c sphericalharm.f90 util.o -o basisFns.o

util.o : util.f90
	gfortran -c util.f90 -lblas -llapack
