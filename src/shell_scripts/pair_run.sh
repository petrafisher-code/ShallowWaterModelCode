#! /bin/bash
# Runs a batch while altering two variables

NAME1="Re"
ARRAY1=(0.8 1.0 1.2) 

NAME2="grav"
ARRAY2=(1366 1466 1566) 

ELEMENTS1=${#ARRAY1[@]} # elements in first array
ELEMENTS2=${#ARRAY2[@]} # elements in second array

cd ..
# make

for (( i=0;i<$ELEMENTS1;i++)); do
	for (( j=0;j<$ELEMENTS2;j++)); do
			# Runs with the hm process switched on:
	 		echo "${NAME1} = ${ARRAY1[${i}]},    ${NAME2} = ${ARRAY1[${i}]}"
            sed -e "s|output.nc|/output_${i}_${j}.nc|" ../config/namelist.in > ../config/namelist.tmp
            sed -e "s/${NAME1}=1.0/${NAME1}=${ARRAY1[${i}]}/" ../config/namelist.tmp > ../config/namelist.tmp2	
            sed -e "s/${NAME2}=1566./${NAME2}=${ARRAY2[${j}]}/" ../config/namelist.tmp2 > ../config/namelist.run
 						
            mpiexec -n 8 ./main.exe ../config/namelist.run > /tmp/std.out
            rm "../config/namelist.tmp"
            mv "/tmp/output_${i}_${j}.nc" "../tests/output.nc"

            # prepare to run python files
            cd python 

            # create height and streamlines figure
            echo "Creating height and streamlines figure"
            python3 height_and_streamlines.py
            mv "../../output/frames/frame.png" "../../output/streamlines/streamlines_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.png" 

            # rename output.nc file for storage
            mv "../../tests/output.nc" "../../tests/output_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${i}]}.nc"
            cd ..
 			
	done
done