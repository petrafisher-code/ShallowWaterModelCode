#! /bin/bash
# Single run

cd ..
# make

mpiexec -n 8 ./main.exe ../config/namelist.in
mv /tmp/output.nc ../tests/output.nc

cd python && python3 animate_output.py
rm ../../output/animations/animation.gif
rm ../../output/animations/animation.mp4

python3 height_and_streamlines.py