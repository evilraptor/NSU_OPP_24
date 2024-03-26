mpicc parallel_slau.c -o parallel_slau -lm

echo 1
mpirun -np 1 ./parallel_slau 1

echo 2
mpirun -np 2 ./parallel_slau 1

echo 8
mpirun -np 8 ./parallel_slau 1
