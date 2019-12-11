#!/bin/bash

#================================================================================
# ET04
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
export PATH="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex11_forward3D_het/"
export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export SENSORS_LOW="sensor_low.dat" 
export SENSORS_MID="sensor_mid.dat" 
export SENSORS_HIGH="sensor_high.dat" 
export FORWARD_SIGNAL="forwardSignal_RT.dat"

# Mode
export MODE="-f"
export GPU_INDEX=0
# Generate dimensions file
Nx=256 dx=0.0001
Ny=128 dy=0.0001
Nz=128 dz=0.0001
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Parameters
nRaysPhi=2
nRaysTheta=30
dt=4e-8
tMax=2e-5

zPos=$(echo "scale=6;($dz*($Nz-2))/2" | bc)
yPos=$(echo "scale=6;($dy*($Ny-2))/2" | bc)

# SENSORS
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.4 0.4 1.1 1.8" >> $SENSORS_LOW

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL


