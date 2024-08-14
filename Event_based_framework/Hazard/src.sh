#!/bin/bash
#SBATCH -n 20
#SBATCH -t 5-00:00:00		# wall time up until 20 hours (5-00:00:00: wall time up until 5 days)
#SBATCH -p fat			# thin (256GiB) or fat (1TiB)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=huazhi.li@vu.nl

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hydromt

num_procs=20 			# number of parallel processes
basin='NEA'

#Parallel simulation
declare -A pids=( )

for cluster in $(seq 1 1 26)
do
  while (( ${#pids[@]} >= num_procs )); do
    sleep 0.2
    for pid in "${!pids[@]}"; do
      kill -0 "$pid" &>/dev/null || unset "pids[$pid]"
    done
  done
#  echo $rp
  python -W ignore Scripts/source_station.py $basin "$cluster" & pids["$!"]=1
done
wait

