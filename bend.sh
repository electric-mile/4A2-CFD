
cd Code
make clean
make
../Code/solver.x < input_bend.txt
cd ..
cd Code
python3 plot_guess.py bend
python3 plot_conv.py bend
python3 plot_contours.py bend


