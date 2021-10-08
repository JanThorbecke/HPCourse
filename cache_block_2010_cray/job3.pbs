#!/bin/bash
#
#PBS -N test3
#PBS -j oe
#PBS -q istanbul
#PBS -o cache_block3_mc6
#PBS -l mppwidth=12
#PBS -S /bin/bash
#PBS -V
#PBS -l walltime=05:05:00

cd /cray/css/users/jan/cache_block_2010
source ${MODULESHOME}/init/bash
module swap PrgEnv-pgi PrgEnv-cray
rm tensor
cp tensor.f90.org tensor.f90
make tensor
mv tensor tensor.org
#cp tensor.f90.opt tensor.f90
#make tensor

aprun -n 1 ./tensor.org
#aprun -n 1 ./tensor

#exit

for (( l=0; l<101; l+=4 )) do
	for (( k=0; k<121; k+=4 )) do
		sed -e 's/555/'"$l"'/g' -e 's/444/'"$k"'/g' < nep3.f90 > tensor.f90
		make tensor
		aprun -n 1 ./tensor > out3.txt
 		ltime=`grep Elapsed < out3.txt | awk '{print $6}'`
 		echo "Blocksize $l x $k $ltime"
 		echo "$l $k $ltime" >> data3.txt
	done
 	echo " " >> data3.txt
done

