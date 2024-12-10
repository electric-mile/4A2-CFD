cd Code
make clean
make
#cd ..
#cd Cases
../Code/solver.x < input_bump.txt
cd ..
cd Code
python3 plot_guess.py bump
python3 plot_conv.py bump
python3 plot_contours.py bump

