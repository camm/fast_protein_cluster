#!/bin/sh
../fast_protein_cluster -i traj -o output_binpos --binary_coords --rmsd --nclusters 10 --nthreads 4 --hcomplete --sse3
