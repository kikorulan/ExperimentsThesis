#!/bin/bash
#================================================================================
# EXAMPLE for STOCHASTIC PDHG
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
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
export SENSORS="sensors_ROI_3600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_noisy5_3600sensors.dat"
export PIXEL_PRESSURE="pixelPressure_0.dat"

# Choose GPU
export GPU_INDEX=0
# Choose mode
export MODE='-G'

# Regularization parameters
TAU=5
LAMBDA=1e-10
NITER=30
#=======   GRADIENT DESCENT
if [ "$MODE" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="

    # Output
    export STDOUT="stdout_GD_tau"$TAU"_lambda"$LAMBDA$"_iter"$NITER".txt"
    RTroiIterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                       $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER
# > $OUTPUT_FOLDER$STDOUT
else
    echo "Non supported mode"
fi

