#!/bin/bash

BASE="EM_matrices_3D_unsym"

# Homogeneous
# for MATRIX in \
    # "homogeneous/250x250x10/cfl-0.25/time-020" \
    # "homogeneous/250x250x10/cfl-0.25/time-100" \
    # "homogeneous/250x250x10/cfl-1.00/time-020" \
    # "homogeneous/250x250x10/cfl-1.00/time-100" \
    # "homogeneous/250x250x10/cfl-4.00/time-020" \
    # "homogeneous/250x250x10/cfl-4.00/time-100"; do
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    # -p multi_split -s 0.8 -b 5 \
    # --dir "results/$BASE/$MATRIX/Msplit(0.8,5)/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    # -p multi_split -s 0.6 -b 5 \
    # --dir "results/$BASE/$MATRIX/Msplit(0.6,5)/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    # -p amg \
    # --dir "results/$BASE/$MATRIX/amg/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    # -p diag \
    # --dir "results/$BASE/$MATRIX/diag/"
# done

# Heterogeneous 3D
# for MATRIX in \
    # "heterogeneous/3D/cfl-1/time-020" \
    # "heterogeneous/3D/cfl-4/time-020" \
    # "heterogeneous/3D/cfl-8/time-020" \
    # "heterogeneous/3D/cfl-8/time-100"; do
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 1 \
    # -p multi_split -s 0.8 -b 5 \
    # --dir "results/$BASE/$MATRIX/Msplit(0.8,5)/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 1 \
    # -p multi_split -s 0.6 -b 5 \
    # --dir "results/$BASE/$MATRIX/Msplit(0.6,5)/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 1 \
    # -p amg \
    # --dir "results/$BASE/$MATRIX/amg/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 1 \
    # -p diag \
    # --dir "results/$BASE/$MATRIX/diag/"
# done

# Heterogeneous 2D
for MATRIX in \
    "heterogeneous/2D/layer_64/cfl-1/time-020"; do 
./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    -p multi_split -s 0.96 -b 5 \
    --dir "results/$BASE/$MATRIX/Msplit(0.96,5)/"
./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    -p multi_split -s 0.9 -b 5 \
    --dir "results/$BASE/$MATRIX/Msplit(0.9,5)/"
./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    -p multi_split -s 0.8 -b 5 \
    --dir "results/$BASE/$MATRIX/Msplit(0.8,5)/"
# ./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    # -p multi_split -s 0.6 -b 5 \
    # --dir "results/$BASE/$MATRIX/Msplit(0.6,5)/"
./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    -p amg \
    --dir "results/$BASE/$MATRIX/amg/"
./spe_prec -u -m ${BASE}/${MATRIX}/matrix.crs -v ${BASE}/${MATRIX}/rhs.crs -a 5 \
    -p diag \
    --dir "results/$BASE/$MATRIX/diag/"
done
