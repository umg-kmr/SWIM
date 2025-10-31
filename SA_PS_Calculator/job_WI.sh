#! /bin/bash
#PBS -N SWIM_test
#PBS -o /home/umang.kumar/SWIM_repo/SA_PS_Calculator/out.log
#PBS -e /home/umang.kumar/SWIM_repo/SA_PS_Calculator/err.log
#PBS -l nodes=1:ppn=64
#PBS -q cpu

#export NUMEXPR_MAX_THREADS=104
export OMP_NUM_THREADS=8

module load compiler/anaconda3
cd /home/umang.kumar/SWIM_repo/SA_PS_Calculator
source activate pyjulia_cobaya
#g++ -shared -fPIC model_calc.cpp -I /home/umang.kumar/boost_1_87_0 -lm -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize -o libmodel.so
mpirun -n 8 cobaya-run /home/umang.kumar/SWIM_repo/SA_PS_Calculator/Input.yaml >> /home/umang.kumar/SWIM_repo/SA_PS_Calculator/output.txt 

