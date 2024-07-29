#! /bin/bash
# Batch running

VEL_ARRAY=(8. 10. 12. 15.) 

VEL_ELEMENTS=${#VEL_ARRAY[@]} # elements in velocity array

cd ..
# make

rm ../output/animations/animation.gif
rm ../output/animations/animation.mp4

for (( i=0;i<$VEL_ELEMENTS;i++)); do
	# Runs with the hm process switched on:
	echo "ujet = ${VEL_ARRAY[${i}]}"
    sed -e "s|output.nc|/output_ujet_${VEL_ARRAY[${i}]}.nc|" ../config/namelist.in > ../config/namelist.tmp
	sed -e "s/u_jet=12./u_jet=${VEL_ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
 			
    if [ -z "$1" ]
    then
        mpiexec -n 1 ./main.exe ../config/namelist.run > /tmp/std.out
    else
                mpiexec -n $1 ./main.exe ../config/namelist.run > /tmp/std.out
    fi
    rm ../config/namelist.tmp
    mv /tmp/output_ujet_${VEL_ARRAY[${i}]}.nc ../tests/output.nc

    # clear old animations and create new ones
    cd python && python3 animate_output.py -y
    mv ../../output/animations/animation.gif ../../output/animations/animation_ujet_${VEL_ARRAY[${i}]}.gif 
    mv ../../output/animations/animation.mp4 ../../output/animations/animation_ujet_${VEL_ARRAY[${i}]}.mp4 

    # create height and streamlines figure
    python3 height_and_streamlines.py
    mv ../../output/frames/frame.png ../../output/streamlines/streamlines_ujet_${VEL_ARRAY[${i}]}.gif 

    # rename output.nc file for storage
    mv ../../tests/output.nc ../../tests/output_ujet_${VEL_ARRAY[${i}]}.nc
    cd ..
done
