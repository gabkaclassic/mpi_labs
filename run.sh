mpicxx main.cpp -o main
sudo mpirun -np $1 ./main $1
sudo rm *.btr