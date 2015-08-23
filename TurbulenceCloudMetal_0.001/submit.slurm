#!/bin/bash
#SBATCH -J TurbulenceCloudMetas0_001              # job name
#SBATCH -o TurbulenceCloudMetas0_001.o%j          # output and error file name (%j expands to jobID)
#SBATCH -n 128                                    # total number of mpi tasks requested
#SBATCH -p normal                                 # queue -- not sure what this means
#SBATCH -t 48:00:00                               # run time
#SBATCH --mail-user=rfernandez@astro.columbia.edu # email
#SBATCH --mail-type=begin                         # email me when the job starts
#SBATCH --mail-type=end                           # email the MPI executable

module load impi
module load hdf5
#ibrun ./enzo.exe -d CloudCollisionMetalsTurbulence.enzo
ibrun ./enzo.exe -d -r DD0061/cloud_collision_0061
