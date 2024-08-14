#!/bin/bash
#SBATCH -n 20
#SBATCH -t 5-00:00:00		# wall time up until 20 hours (5-00:00:00: wall time up until 5 days)
#SBATCH -p fat			# thin (256GiB) or fat (1TiB)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=huazhi.li@vu.nl

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hydromt

num_procs=20 			# number of parallel processes
basin='SEA'

python -W ignore Scripts/source_station_basin.py $basin

wait

