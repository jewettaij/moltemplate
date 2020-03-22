
mpiexec -np 36 lmp_mpi -i run.in.min
mpiexec -np 36 lmp_mpi run_increase_Udd=U0cd.in
