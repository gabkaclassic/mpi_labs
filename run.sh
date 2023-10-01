mpicxx main.cpp -o main
sudo mpirun -np $1 ./main
sudo rm *.btr