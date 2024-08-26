#! /bin/bash
# Runs a batch while altering two variables

# NAME1="Re"
# ARRAY1=(0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1 1.02 1.04 1.06 1.08) 

# NAME2="grav"
# ARRAY2=(1344 1368 1392 1416 1440 1464 1488 1512 1536 1560 1584 1608 1632 1656 1680 1704 1728 1752 1776 1800 1824) 

NAME1="nlat_thresh"
ARRAY1=(86.5 86.6 86.7 86.8 86.9 87.0 87.1 87.2 87.3 87.4 87.5 87.6 87.7 87.8 87.9 88.0)

NAME2="u_vortex"
ARRAY2=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)

ELEMENTS1=${#ARRAY1[@]} # elements in first array
ELEMENTS2=${#ARRAY2[@]} # elements in second array

cd .. && make

# rm ../output/animations/animation.gif
# rm ../output/animations/animation.mp4

for (( i=0; i<ELEMENTS1; i++)); do
	for (( j=0; j<ELEMENTS2; j++)); do
			# Runs with the hm process switched on:
	 		echo "${NAME1} = ${ARRAY1[${i}]},    ${NAME2} = ${ARRAY2[${j}]}"
            sed -e "s|output.nc|/output_${i}_${j}.nc|" ../config/namelist.in > ../config/namelist.tmp
            # sed -e "s/${NAME1}=1.0/${NAME1}=${ARRAY1[${i}]}/" ../config/namelist.tmp > ../config/namelist.tmp2	
            # sed -e "s/${NAME2}=1566./${NAME2}=${ARRAY2[${j}]}/" ../config/namelist.tmp2 > ../config/namelist.run
            sed -e "s/${NAME1}=86.5/${NAME1}=${ARRAY1[${i}]}/" ../config/namelist.tmp > ../config/namelist.tmp2	
            sed -e "s/${NAME2}=7./${NAME2}=${ARRAY2[${j}]}/" ../config/namelist.tmp2 > ../config/namelist.run
 						
            mpiexec -n 8 ./main.exe ../config/namelist.run > /tmp/std.out
            rm "../config/namelist.tmp"
            mv "/tmp/output_${i}_${j}.nc" "../tests/output.nc"

            # prepare to run python files
            cd python 

            # # clear old animations and create new ones
            # echo "Creating animations"
            # python3 animate_output.py -y
            # mv "../../output/animations/animation.gif" "../../output/animations/animation_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.gif" 
            # mv "../../output/animations/animation.mp4" "../../output/animations/animation_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.mp4" 

            # create height and streamlines figure
            echo "Creating height and streamlines figure"
            python3 height_and_streamlines.py
            mv "../../output/frames/frame.png" "../../output/streamlines/streamlines_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.png" 

            # create height and streamlines figure
            echo "Creating velocity latitude figure"
            python3 velocity_latitudes.py
            mv "../../output/velocities/velocity_latitudes.png" "../../output/velocities/velocity_latitudes_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.png"

            # rename output.nc file for storage
            mv "../../tests/output.nc" "../../tests/output_${NAME1}_${ARRAY1[${i}]}_${NAME2}_${ARRAY2[${j}]}.nc"
            cd ..
 			
	done
done