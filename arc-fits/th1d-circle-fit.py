import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

fits_file = "/Users/will/Dropbox/JorgeBowshocks/MIR/SmithBally/wcs.mos11jy.fits"

for delta_theta in 65, 70, 75, 80:
    plotfile = f"th1D-{delta_theta:02d}.pdf"
    print('#### '*10)
    print("Creating", plotfile)
    circle_fit.plot_solution(
        f"th1d-ori-smith-2005-forma.reg",
        fits_file,
        plotfile,
        delta_theta=delta_theta,
        vmin=0.0, vmax=0.01,
    )
