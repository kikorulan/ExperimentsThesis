#!/bin/bash

#====================
# QSUB
#====================

#$ -P gpu
#$ -l gpu=1
#$ -l h_rt=10:50:0
#$ -l tmem=3G
#$ -N RTsolver
#$ -wd /home/frullan/C++
#$ -S /bin/bash

#$ -o RTsolver.txt
#$ -j y

#================================================================================
# Real Data 06 : finger
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex17_real_recon3D_DS2/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ] || [ "$HOSTNAME" = "ember.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
elif [ "$HOSTNAME" = "hannover" ]; then
    export HOST_FOLDER="/home/wontek/sharedWK/ExperimentsThesis/"
else
    export HOST_FOLDER="/home/frullan/HighFreqCode/ExperimentsThesis/"
fi
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="pixelPressure_0.dat"
export SENSORS="sensors_subsampled_3600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_3600sensors_490timesteps.dat"
export STDOUT="stdout-adjoint.txt"

# Mode
MODE='-a'
export GPU_INDEX=0
# Generate dimensions file
Nx=80  dx=0.000053
Ny=240 dy=0.000053
Nz=240 dz=0.000053
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

#==============================
# SENSORS
#==============================
sensors_y=60
sensors_z=60
nRaysPhi=1024 
nRaysTheta=1024
dt=1.6667e-8
tMax=8.16e-6

# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" > $INPUT_FOLDER$SENSORS

subsampleFactor=4
sensors_y=$(echo "scale=0;($Ny/$subsampleFactor)" | bc)
sensors_z=$(echo "scale=0;($Nz/$subsampleFactor)" | bc)
echo $sensors_y $sensors_z
# YZ
for ((k=0; k<$sensors_z; k++)); do
    zPos=$(echo "scale=6;($k*$dz*$subsampleFactor)" | bc)
    for ((i=0; i<$sensors_y; i++)); do
        yPos=$(echo "scale=6;($i*$dy*$subsampleFactor)" | bc)
        echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done 
done 

#==============================
# SENSORS
#==============================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
# > $OUTPUT_FOLDER$STDOUT
