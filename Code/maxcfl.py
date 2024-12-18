import os
import subprocess
import numpy as np
from routines import *

def modify_input_parameters(filename, new_cfl=None, new_sfac=None):
    """
    Modify the input file to update CFL and sfac values.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    if len(lines) >= 3:  # Ensure there are at least 3 lines
        line = lines[2]  # Line indexing starts at 0
        parts = line.split()  # Split the line into parts

        if new_cfl is not None:
            parts[0] = f"{new_cfl:.3f}"  # Update CFL

        if new_sfac is not None:
            if len(parts) > 1:  # Ensure the second part exists
                parts[1] = f"{new_sfac:.3f}"  # Update sfac

        lines[2] = ' '.join(parts) + '\n'  # Reconstruct the line

    with open(filename, 'w') as file:
        file.writelines(lines)

def run_solver(input_file, log_file):
    with open(input_file, 'r') as inp, open(log_file, 'w') as log:
        subprocess.run(['../Code/solver.x'], stdin=inp, stdout=log, stderr=subprocess.STDOUT)

def analyze_convergence(file_path, tolerance=1e-3):
    """
    Analyze convergence history from a space-separated file.
    Returns True if stable (residual decreases consistently below tolerance), False otherwise.
    """
    try:
        # Load the data using space as the delimiter
        data = np.loadtxt(file_path, delimiter=None)  # Automatically handles space-separated values
        # Check if data is 2-dimensional
        if data.ndim == 1:
            print(f"Data in {file_path} is not in the expected format.")
            return "format issue"
        residuals = data[:, 1]  # Assuming the second column contains residuals
        
        # Check if residuals consistently drop below the tolerance
        if np.all(residuals[-10:] < tolerance):  # Check the last 5 iterations
            return True
    except Exception as e:
        print(f"Error analyzing {file_path}: {e}")
    return False

def main():
    input_file = 'input_bend.txt'
    settings_file = 'input_bend.txt'
    conv_file = 'conv_bend.csv'

    # Initialize with CFL = 0.4 and sfac = 0.3
    initial_cfl = 0.4
    initial_sfac = 0.3
    modify_input_parameters(settings_file, new_cfl=initial_cfl, new_sfac=initial_sfac)

    # Test increasing CFL values
    cfl_values = np.arange(0.4, 10, 0.2)  # Increments of 0.02 from 0.4 to 1.0
    sfac_values = np.arange(0.3, 0.00, -0.02)  # Decrementing by 0.02 from 0.3 to 0.1

    max_stable_cfl = None
    min_stable_sfac = None
 
    # Step 1: Test for maximum stable CFL
    for cfl in cfl_values:
        print(f'Testing CFL = {cfl:.2f}')
        modify_input_parameters(settings_file, new_cfl=cfl)
        log_file = f'log_bend_CFL_{cfl:.2f}.txt'
        count  = 0.0
        for i in range(3):
            run_solver(input_file, log_file)
            result = analyze_convergence(conv_file)
            if result:
                print(f'CFL = {cfl:.2f} is stable.')
                count += 1
            elif result == "format issue":
                print("Retrying due to format issue.")
                i -= 1
            else:
                print(f'CFL = {cfl:.2f} is not stable.')
        if count < 2:
            print('Stopping further tests.')
            break
        else:
            max_stable_cfl = cfl

    if not max_stable_cfl:
        print("No stable CFL found.")
        return

    print(f"Maximum stable CFL: {max_stable_cfl:.2f}")

    #modify_input_parameters(settings_file, new_cfl=initial_cfl, new_sfac=initial_sfac)

    # Step 2: Test for minimum stable sfac
    for sfac in sfac_values:
        print(f'Testing sfac = {sfac:.2f} with CFL = {max_stable_cfl:.2f}')
        modify_input_parameters(settings_file, new_cfl=max_stable_cfl, new_sfac=sfac)
        log_file = f'log_bend_CFL_{max_stable_cfl:.2f}_sfac_{sfac:.2f}.txt'
        count = 0.0
        for i in range(3):

            run_solver(input_file, log_file)
            result = analyze_convergence(conv_file)
            if result:
                print(f'sfac = {sfac:.2f} is stable.')
                count += 1
            elif result == "format issue":
                print("Retrying due to format issue.")
                i -= 1
            else:
                print(f'sfac = {sfac:.2f} is not stable.')
        if count < 2:
            print('Stopping further tests.')
            break
        else:
            min_stable_sfac = sfac

    if min_stable_sfac:
        print(f"Minimum stable sfac: {min_stable_sfac:.2f}")
    else:
        print("No stable sfac found.")


if __name__ == '__main__':
    main()