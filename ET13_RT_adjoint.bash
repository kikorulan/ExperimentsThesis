#!/bin/bash

#================================================================================
# ET04
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
export PATH="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex13_adjoint3D_het/"
export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SENSORS_LOW="sensor_low.dat" 
export SENSORS_MID="sensor_mid.dat" 
export SENSORS_HIGH="sensor_high.dat" 

# Mode
export MODE="-a"
export GPU_INDEX=0
# Generate dimensions file
Nx=256 dx=0.0001
Ny=128 dy=0.0001
Nz=128 dz=0.0001
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

# Sound speed
export SOUND_SPEED="sound_speed.dat"

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
export FORWARD_SIGNAL_LOW="forwardSignal_kWave_low_"$dt".dat"
export FORWARD_SIGNAL_MID="forwardSignal_kWave_mid_"$dt".dat"
export FORWARD_SIGNAL_HIGH="forwardSignal_kWave_high_"$dt".dat"
echo "==================================== dt = 4e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -.5 .5 1.17 1.97" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1 1 0.6 2.54" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.3 2.84" >> $SENSORS_HIGH

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL_LOW
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_low_"$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL_MID
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_mid_"$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL_HIGH
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_high_"$dt.dat
  
#==============================
# SENSORS
#==============================
dt=2e-8
export FORWARD_SIGNAL_LOW="forwardSignal_kWave_low_"$dt".dat"
export FORWARD_SIGNAL_MID="forwardSignal_kWave_mid_"$dt".dat"
export FORWARD_SIGNAL_HIGH="forwardSignal_kWave_high_"$dt".dat"
echo "==================================== dt = 2e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -.5 .5 1.17 1.97" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1 1 0.6 2.54" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.3 2.84" >> $SENSORS_HIGH

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL_LOW
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_low_"$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL_MID
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_mid_"$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL_HIGH
mv $OUTPUT_FOLDER"/PixelPressure.dat"  "adjoint_pressure_RT_high_"$dt.dat

#==============================
# SENSORS
#==============================
dt=1e-8
export FORWARD_SIGNAL_LOW="forwardSignal_kWave_low_"$dt".dat"
export FORWARD_SIGNAL_MID="forwardSignal_kWave_mid_"$dt".dat"
export FORWARD_SIGNAL_HIGH="forwardSignal_kWave_high_"$dt".dat"
echo "==================================== dt = 1e-8 ========================================="
# SENSORS LOW
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_LOW
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -0.5 0.5 1.17 1.97" >> $SENSORS_LOW
# SENSORS MID
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_MID
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1 1 0.6 2.54" >> $SENSORS_MID
# SENSORS HIGH
echo "$dt $tMax 0 0 0 0 0 0 0" > $SENSORS_HIGH
echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.3 2.84" >> $SENSORS_HIGH

#====================
# RUN 
#====================
# Low pressure
export INITIAL_PRESSURE="u0_low.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_LOW $OUTPUT_FOLDER $FORWARD_SIGNAL_LOW
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_low_"$dt.dat
# Mid pressure
export INITIAL_PRESSURE="u0_mid.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_MID $OUTPUT_FOLDER $FORWARD_SIGNAL_MID
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_mid_"$dt.dat
# High pressure
export INITIAL_PRESSURE="u0_high.dat"
RTsolver_GPU $MODE $DIMENSIONS $SOUND_SPEED $INITIAL_PRESSURE $SENSORS_HIGH $OUTPUT_FOLDER $FORWARD_SIGNAL_HIGH
mv $OUTPUT_FOLDER"/PixelPressure.dat" "adjoint_pressure_RT_high_"$dt.dat
