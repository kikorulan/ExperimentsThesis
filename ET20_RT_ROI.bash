#!/bin/bash
#================================================================================
# EXAMPLE 84
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex20_ROI/"
# Output folder
export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER


# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_ROI_3600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_noisy5_3600sensors.dat"
export STDOUT="stdout-forward.txt"

# Mode
export MODE="-a"
export GPU_INDEX=0

#====================
# RUN 
#====================
RTroi_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
          $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL
