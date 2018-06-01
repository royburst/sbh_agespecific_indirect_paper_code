#!/bin/sh
SINGULARITYENV_OMP_NUM_THREADS=10 /share/local/singularity-2.4.2/bin/singularity exec --cleanenv /share/singularity-images/lbd/test_deploy/r_pkgs3.4.3gcc7mkl.simg /usr/local/bin/R <$1 --no-save $@
