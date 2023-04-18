qlogin -l h_vmem=16G

wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
tar zxvf julia-1.8.5-linux-x86_64.tar.gz

export PATH="$PATH:~/julia-1.8.5/bin"
export JULIA_NUM_THREADS=10