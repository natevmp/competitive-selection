#!/bin/bash
#$ -pe smp 10 # Request CPU cores
#$ -l h_vmem=2G
#$ -l h_rt=16:0:0
#$ -wd ~/CompetitiveSelection
#$ -j y
#$ -N sizeDist
#$ -o ~/CompetitiveSelection
#$ -m beas

export PATH="$PATH:~/julia-1.8.5/bin"
export JULIA_NUM_THREADS=10

SRCDIR=$HOME/CompetitiveSelection/
DATADIR=$HOME/CompetitiveSelection/SimResults/
mkdir -p $DATADIR


rsync -rt $SRCDIR $TMPDIR/
# rsync -rltv $SRCDIR/src $TMPDIR/
# rsync -rltv $SRCDIR/HPC $TMPDIR/
# rsync -rltv $SRCDIR/Data $TMPDIR/

cd $TMPDIR

echo "starting julia script..."
# julia HPC/addPackages.jl
julia HPC/sizeDistABC.jl "./params/paramsChisquare6575.jl"

mv abcResult*.jld2 $DATADIR/

echo "success"