#!/bin/bash

#====================
# QSUB
#====================

#$ -l gpu=true
#$ -l h_rt=100:00:00
#$ -l tmem=3G
#$ -N ET15_forward
#$ -S /bin/bash

#$ -o RTsolver.txt
#$ -j y

#================================================================================
# ET15
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
export EXAMPLE="Ex15_synth_recon3D_het/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ] || [ "$HOSTNAME" = "blaze.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentsThesis/"
elif [ "$HOSTNAME" = "ember.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/scratch0/NOT_BACKED_UP/frullan/ExperimentsThesis/"
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
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled_1.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_1sensor.dat"
export STDOUT="stdout-forward.txt"

# Mode
export MODE="-a"
export GPU_INDEX=0

#====================
# RUN 
#====================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
# > $OUTPUT_FOLDER$STDOUT
if [ "$MODE" = "-f" ]; then
    mv $OUTPUT_FOLDER"ForwardSignal.dat" $INPUT_FOLDER"forwardSignal_reference_1sensor.dat"
fi
#mv $OUTPUT_FOLDER"PixelPressure.dat" $OUTPUT_FOLDER"pressure_adjoint_RT.dat"
