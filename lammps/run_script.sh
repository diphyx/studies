#!/bin/bash

# Ensure the script stops on any error
set -e

# Navigate to the data directory
cd /data

# Unzip the amoeba.zip file into the amoeba directory
echo "Unzipping amoeba.zip..."
unzip -o amoeba.zip -d amoeba

# Change directory to /data/amoeba
echo "Changing directory to /data/amoeba..."
cd /data/amoeba

# water_dimer conversions
echo "Converting water dimer data..."
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_dimer.xyz -amoeba amoeba_water.prm -data data.water_dimer.amoeba
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_dimer.xyz -hippo hippo_water.prm -data data.water_dimer.hippo

# water_hexamer conversions
echo "Converting water hexamer data..."
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_hexamer.xyz -amoeba amoeba_water.prm -data data.water_hexamer.amoeba
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_hexamer.xyz -hippo hippo_water.prm -data data.water_hexamer.hippo

# water_box conversions
echo "Converting water box data..."
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_box.xyz -amoeba amoeba_water.prm -data data.water_box.amoeba -pbc 18.643 18.643 18.643
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz water_box.xyz -hippo hippo_water.prm -data data.water_box.hippo -pbc 18.643 18.643 18.643

# ubiquitin conversions
echo "Converting ubiquitin data..."
python $LAMMPS_TOOLS_PATH/tinker/tinker2lmp.py -xyz ubiquitin.xyz -amoeba amoeba_ubiquitin.prm -data data.ubiquitin -pbc 54.99 41.91 41.91 -bitorsion bitorsion.ubiquitin.data

echo "All conversions completed successfully."
