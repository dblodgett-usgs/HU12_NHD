#!/bin/sh

# Sample submission script for pbdMPI Hello World

#SBATCH --account=geoint
#SBATCH -p normal 
#SBATCH --time=03:00:00
#SBATCH --job-name=hello_world
#SBATCH -o hello.out
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dblodgett@usgs.gov

# Add modules
module purge
module load gdal/2.4.1-gcc71 gcc/7.1.0 gis/geos-3.5.0 
module load proj/4.9.2-gcc6.1.0 tools/udunits-2.2.16-gnu

# Run R
srun Rscript hr_wbd_runner.R