import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

fits_file = "../data/LP-Ori-HST/lp-ori-f435w.fits"

for delta_theta in 110, 130, 150:
    plotfile = f"lp-ori-{delta_theta:03d}.pdf"
    print('#### '*10)
    print("Creating", plotfile)
    circle_fit.plot_solution(
        f"lp-ori-acs-forma.reg",
        fits_file,
        plotfile,
        delta_theta=delta_theta,
        vmin=0.0, vmax=8.0,
    )
