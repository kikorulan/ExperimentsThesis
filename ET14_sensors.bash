#!/bin/bash

#================================================================================
# ET14
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex14_synth_recon3D_homo/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ] || [ "$HOSTNAME" = "blaze.cs.ucl.ac.uk" ] || [ "$HOSTNAME" = "ember.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
else
    export HOST_FOLDER="/home/frullan/HighFreqCode/ExperimentsThesis/"
fi
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SENSORS="sensors_subsampled_14400.dat" 

#==============================
# DIMENSIONS
#==============================
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
nRaysPhi=1024 
nRaysTheta=1024
dt=1.667e-8
tMax=8.0836e-06

# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" > $INPUT_FOLDER$SENSORS

subsampleFactor=2
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
