#!/bin/bash
#$ -pe smp 1
#$ -pe smp 10  # Request CPU cores
#$ -l h_vmem=12G
#$ -l h_rt=12:0:0
#$ -wd ~/CompetitiveSelection
#$ -j y
#$ -N cloneNumbers
#$ -o ~/CompetitiveSelection
#$ -m beas

module load julia

SRCDIR=$HOME/CompetitiveSelection/
DATADIR=$HOME/CompetitiveSelection/Data/
# mkdir -p $DATADIR

rsync -rltv $SRCDIR/src $TMPDIR/
rsync -rltv $SRCDIR/HPC $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

# julia HPC/addPackages.jl
julia HPC/cloneNumbersAlphaMuSpace.jl

# rsync -rltv $TMPDIR/ $DATADIR/
mv sizeDistABCResult.jld2 $DATADIR/

echo "success"