#!/bin/bash

#====================
# QSUB
#====================

#$ -l gpu=true
#$ -l h_rt=100:00:00
#$ -l tmem=3G
#$ -N RTsolver
#$ -S /bin/bash

#$ -o RTsolver.txt
#$ -j y

#================================================================================
# EXAMPLE 85
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
if [ "$HOSTNAME" = "miller.local" ] || [ "$HOSTNAME" = "armstrong.local" ]; then
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build_miller/bin:$PATH"
else
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
fi
export EXAMPLE="Ex04_simulation3D_homo/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ] || [ "$HOSTNAME" = "blaze.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
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
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240_smooth.dat"
export SENSORS="sensors_subsampled_180k.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_noisy5_14400sensors.dat"
export STDOUT="stdout-forward.txt"

# Mode
export MODE="-f"
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
nRaysPhi=1024 
nRaysTheta=1024
dt=1.6667e-8
tMax=8.0836e-06

# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS

# XY Bottom
zPos=$(echo "scale=6; 0" | bc)
for ((j=0; j<Ny; j++)); do
    for ((i=0; i<Nx; i++)); do
        xPos=$(echo "scale=6;($i*$dx)" | bc)
        yPos=$(echo "scale=6;($j*$dy)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta -3.14 3.14 0.04 1.57" >> $INPUT_FOLDER$SENSORS
    done
done

# XY Bottom
zPos=$(echo "scale=6; 0" | bc)
for ((j=0; j<Ny; j++)); do
    for ((i=0; i<Nx; i++)); do
        xPos=$(echo "scale=6;($i*$dx)" | bc)
        yPos=$(echo "scale=6;($j*$dy)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta -3.14 3.14 0.04 1.57" >> $INPUT_FOLDER$SENSORS
    done
done

# XZ - YZ
for ((k=1; k<Nz; k++)); do
    zPos=$(echo "scale=6;$k*$dz" | bc)
    # XZ
    yPos=$(echo "scale=6; 0" | bc)
    for ((i=0; i<Nx; i++)); do
        xPos=$(echo "scale=4;($i*$dx)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 0 3.14 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done
    # XY 
    xPos=$(echo "scale=6;($dx*($Nx-1))" | bc)
    for ((j=1; j<Ny-1; j++)); do
        yPos=$(echo "scale=6;($j*$dy)" | bc)
        echo "0     $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta  1.57 4.71 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done 
    # XZ
    yPos=$(echo "scale=6;($dy*($Ny-1))" | bc)
    for ((i=0; i<Nx; i++)); do
        xPos=$(echo "scale=6;($i*$dx)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta 3.14 6.28 0.04 3.1" >> $INPUT_FOLDER$SENSORS
    done
done 

# XY Top
zPos=$(echo "scale=6;($dz*($Nz-1))" | bc)
for ((j=0; j<Ny; j++)); do
    for ((i=0; i<Nx; i++)); do
        xPos=$(echo "scale=6;($i*$dx)" | bc)
        yPos=$(echo "scale=6;($j*$dy)" | bc)
        echo "$xPos $yPos $zPos $nRaysPhi $nRaysTheta -3.14 3.14 1.57 3.1" >> $INPUT_FOLDER$SENSORS
    done
done

#====================
# RUN 
#====================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL > $OUTPUT_FOLDER$STDOUT
