# Script to run Uintah-Wasatch on a specified
cd /data

# Define variables
result_dir="uintah_wasatch_varden_2D" # Name of the result directory
input_file="varden-projection-2d-oscillating-periodic-mms-xy.ups" # Path to the input .fastq file


# Welcome message
echo "Running Uintah-Wasatch Simulation on $input_file..."

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist."
    exit 1
fi

sus $input_file 

# Check if Uintah-Wasatch ran successfully
if [ $? -eq 0 ]; then
    echo "Uintah-Wasatch simulation completed successfully."
    echo "Check the $result_dir directory for the results."
else
    echo "Uintah-Wasatch failed!"
    exit 1
fi
