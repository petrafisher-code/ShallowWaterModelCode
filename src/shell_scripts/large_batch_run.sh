#! /bin/bash
# Batch running - large


ALL_NAMES=("u_jet" "h_jet" "theta_jet" "jet_noise" "perturb_strength")
ELEMENTS_ALL_NAMES=${#ALL_NAMES[@]} # elements in array

# cd .. && make
# rm ../output/animations/animation.gif
# rm ../output/animations/animation.mp4

for((k=0; k<ELEMENTS_ALL_NAMES; k++)); do
    cd ~/Documents/Petra/ShallowWaterModelCode/src

    NAME=${ALL_NAMES[${k}]}
    if test "${NAME}" = "u_jet"; then
        ARRAY=(6. 6.5 7. 7.5 8. 8.5 9. 9.5 10. 10.5 11. 11.5 12. 12.5 13. 13.5 14. 14.4 15. 15.5 16. 16.5 17. 17.5 18. 18.5 19. 19.5 20. 21. 22. 23. 24. 25. 26. 27. 28. 29. 30. 31. 32.)
    elif test "${NAME}" = "h_jet"; then
        ARRAY=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95)
    elif test "${NAME}" = "theta_jet"; then
        ARRAY=(70. 70.5 71. 71.5 72. 72.5 73. 73.5 74. 74.5 75. 75.5 76. 76.5 77. 77.5 78. 78.5 79. 79.5 80. 70.25 70.75 71.25 71.75 72.25 72.75 73.25 73.75 74.25 74.75 75.25 75.75 76.25 76.75 77.25 77.75 78.25 78.75 79.25 
            79.75 80.25)
    elif test "${NAME}" = "jet_noise"; then
        ARRAY=(0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.000000001 0.0000000001 0.00000000001 0.000000000001 0.0000000000001 0.00000000000001 0.000000000000001 0.0000000000000001 0.00000000000000001 
            0.000000000000000001 0.0000000000000000001 0.00000000000000000001 0.05 0.005 0.0005 0.00005 0.000005 0.0000005 0.00000005 0.000000005 0.0000000005 0.00000000005 0.000000000005 0.0000000000005 0.00000000000005 
            0.000000000000005 0.0000000000000005 0.00000000000000005 0.000000000000000005 0.0000000000000000005 0.00000000000000000005 0.000000000000000000005)
    elif test "${NAME}" = "perturb_strength"; then
        ARRAY=(0.1 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25 8.75 9.25 9.75)
    fi
    ELEMENTS=${#ARRAY[@]} # elements in array

    for (( i=0; i<ELEMENTS; i++)); do
        # Runs with the hm process switched on:
        echo "${NAME} = ${ARRAY[${i}]}"
        sed -e "s|output.nc|/output_${NAME}_${ARRAY[${i}]}.nc|" ../config/namelist.in > ../config/namelist.tmp
        if test "${NAME}" = "u_jet"; then
            sed -e "s/${NAME}=15./${NAME}=${ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
        elif test "${NAME}" = "h_jet"; then
            sed -e "s/${NAME}=1.0/${NAME}=${ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
        elif test "${NAME}" = "theta_jet"; then
            sed -e "s/${NAME}=75.9/${NAME}=${ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
        elif test "${NAME}" = "jet_noise"; then
            sed -e "s/${NAME}=0.01/${NAME}=${ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
        elif test "${NAME}" = "perturb_strength"; then
            sed -e "s/${NAME}=0.5/${NAME}=${ARRAY[${i}]}/" ../config/namelist.tmp > ../config/namelist.run	
        fi
                
        mpiexec -n 8 ./main.exe ../config/namelist.run > /tmp/std.out
        rm "../config/namelist.tmp"
        mv "/tmp/output_${NAME}_${ARRAY[${i}]}.nc" "../tests/output.nc"

        # prepare to run python files
        cd python 

        # create height and streamlines figure
        echo "Creating height and streamlines figure"
        python3 height_and_streamlines.py
        mv "../../output/frames/frame.png" "../../output/streamlines/streamlines_${NAME}_${ARRAY[${i}]}.png" 

        # rename output.nc file for storage
        mv "../../tests/output.nc" "../../tests/output_${NAME}_${ARRAY[${i}]}.nc"
        cd ..
    done
done