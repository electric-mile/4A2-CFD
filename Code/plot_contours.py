#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
    cut = cut_i(g,0)
    pstagref, mass = mass_av(cut, 'pstag')
    pref, l = area_av(cut, 'p')
    pstag_in, l = area_av(cut, 'pstag')
    g['cp'] = (g['p'] - pref) / (pstagref - pref)
    g['pstag_coeff'] = (g['pstag'] - pstag_in) / (pstag_in - pref)
    tstag_ratio = []
    mass_ratio = []
    for i in range(g['ni']):
        mass_av_tstag, mass = mass_av(cut_i(g,i), 'tstag_ratio')
        tstag_ratio.append(mass_av_tstag)
        mass_ratio.append(np.sum(mass))

    
    plt.plot(np.linspace(1,g['ni'],g['ni']),tstag_ratio)
    plt.show()
    plt.plot(np.linspace(1,g['ni'],g['ni']),mass_ratio/mass_ratio[0])
    plt.show()
        

    # Specify the parameters to plot
    fieldnames = ['cp', 'mach', 'pstag_coeff']
    colnames = ['Static pressure coefficient', 'Mach number', 'Pstag coefficient']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],np.linspace(-1.0,1.0,40),colors='w',
                linewidths=0.5)
        if name == 'pstag_coeff':
            ax.contour(g['x'],g['y'],g['pstag_coeff'],np.linspace(-1.0,1.0,40),colors='w',
                linewidths=0.5)

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.show()

    
main()


