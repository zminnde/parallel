#!/bin/bash
# pin each MPI rank (and the Python child it forks, via inherited affinity) to one core
exec taskset -c $((OMPI_COMM_WORLD_RANK * 2)) "$@"
