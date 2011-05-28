#!/bin/bash

# Homogeneous
# DIR="results/EM_matrices_3D_unsym/homogeneous/250x250x10/cfl-0.25/time-020"

# Heterogeneous 3D
# DIR="results/EM_matrices_3D_unsym/heterogeneous/3D/cfl-1/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/3D/cfl-4/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/3D/cfl-8/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/3D/cfl-8/time-100"

# Heterogeneous 2D
# DIR="results/EM_matrices_3D_unsym/heterogeneous/2D/layer_1/cfl-1/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/2D/layer_22/cfl-1/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/2D/layer_43/cfl-1/time-020"
DIR="results/EM_matrices_3D_unsym/heterogeneous/2D/layer_64/cfl-1/time-020"
# DIR="results/EM_matrices_3D_unsym/heterogeneous/2D/layer_85/cfl-1/time-020"

pushd .
cd $DIR

for prec in *; do
    echo "    $prec &" `cat $prec/t_const` "&" `cat $prec/t_prec` "&" `cat $prec/niter` "&" `cat $prec/t_sol` \
	"& {\\bf " `cat $prec/t_total` "} \\\\" >> tabular.tex
    echo "    \hhline{------}" >> tabular.tex
done

echo '\begin{table}[h!]' > table.tex
echo "\\caption{$DIR}" >> table.tex
echo '\begin{center}' >> table.tex
echo '    \begin{tabular}{c||r|r|r|r|r}' >> table.tex
echo '    Prec & $t_{const}$ & $t_{prec}$ & \# iter & $t_{sol}$ & $t_{total}$ \\' >> table.tex
echo '    \hhline{======}' >> table.tex
cat tabular.tex >> table.tex
rm tabular.tex
echo '    \end{tabular}' >> table.tex
echo '\end{center}' >> table.tex
echo '\end{table}' >> table.tex

popd
cp -f $DIR/table.tex .
