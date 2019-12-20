#!/bin/bash

#================================================================================
# ET04
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
export PATH="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex10_forward3D/"
export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SENSORS_LOW="sensor_low.dat" 
export SENSORS_MID="sensor_mid.dat" 
export SENSORS_HIGH="sensor_high.dat" 
export FORWARD_SIGNAL="forwardSignal_RT.dat"

# Mode
export MODE="-f"
export GPU_INDEX=0
export RESOLUTION=3
# Generate dimensions file
case $RESOLUTION in
    1)
        echo "RESOLUTION 1024x512x512"
        Nx=1024 dx=0.000025
        Ny=512 dy=0.000025
        Nz=512 dz=0.000025
        ;;
    2)
        echo "RESOLUTION 512x256x256"
        Nx=512 dx=0.00005
        Ny=256 dy=0.00005
        Nz=256 dz=0.00005
        ;;
    3)
        echo "RESOLUTION 256x128x128"
        Nx=256 dx=0.0001
        Ny=128 dy=0.0001
        Nz=128 dz=0.0001
        ;;
esac
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Sound speed 
export SOUND_SPEED="sound_speed_1500_"$Nx$".dat"

# Parameters
nRaysPhi=1024
nRaysTheta=1024
tMax=2e-5

#zPos=$(echo "scale=6;($dz*($Nz-2))/2" | bc)
#yPos=$(echo "scale=6;($dy*($Ny-2))/2" | bc)
zPos=0.0063
yPos=0.0063
#==============================
# SENSORS
#==============================
dt=4e-8
echo "==================================== dt = 4e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5 0.5 1.0 2.2" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6 0.6 1.0 2.2" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0005 2.2005" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0005 2.2005" >> $SENSORS_HIGH

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_low_$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_mid_$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_high_$dt.dat
  
#==============================
# SENSORS
#==============================
dt=2e-8
echo "==================================== dt = 2e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5 0.5 1.0 2.2" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6 0.6 1.0 2.2" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0005 2.2005" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0005 2.2005" >> $SENSORS_HIGH

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_low_$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_mid_$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_high_$dt.dat

#==============================
# SENSORS
#==============================
dt=1e-8
echo "==================================== dt = 1e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5 0.5 1.0 2.2" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6 0.6 1.0 2.2" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0000 2.2000" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5995 0.6000 1.0005 2.2005" >> $SENSORS_HIGH
#echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.6000 0.6005 1.0005 2.2005" >> $SENSORS_HIGH
#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_low_$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_mid_$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high_"$Nx".dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL
mv $OUTPUT_FOLDER"/ForwardSignal.dat" forwardSignal_RT_high_$dt.dat
