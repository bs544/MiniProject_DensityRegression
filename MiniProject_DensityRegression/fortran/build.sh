#!/bin/bash
#kind_map taken from Andrew Fowler's github page:
#https://github.com/andrew31416/densityregression
#in features

modname="f90_descriptor"

file1="neighbours"
file2="util"
file3="basis"
file4="unitTests"
file5="bispectrum"

#gfortran_opts="-W -Wall -fPIC -llapack -lblas -O2"

rm .f2py_f2cmap
rm f90_descriptor.cpython-36m-x86_64-linux-gnu.so
for file in {$file1,$file2,$file3,$file4,$file5}
do
	for suffix in {"o","mod"}
	do
		rm $file"."$suffix
	done
	rm "f90wrap_"$file".f90"
	
done
make

python_facing_fns="get_CG_tensor"
python_facing_fns=$python_facing_fns" getLocalBispectrum"
python_facing_fns=$python_facing_fns" CG"


python_facing_files="bispectrum.f90"




export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3

f90wrap -m $modname $python_facing_files -k kind_map -S 12 --only $python_facing_fns
f2py -c -m $modname -lblas f90wrap*.f90 *.o --f90flags="-fPIC -llapack -lblas -O2" --fcompiler=gfortran

#rm $modname".py"

