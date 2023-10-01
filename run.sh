mpicxx main.cpp -o main
sudo mpirun -np $1 ./main $2
sudo rm *.btr