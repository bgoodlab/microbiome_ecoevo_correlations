import os
import sys

sys.stderr.write("Plotting Figure 1 & supplements...\n")
os.system("python plot_figure_1.py")
sys.stderr.write("Done!")

sys.stderr.write("Plotting Figure 2 & supplements...\n")
os.system("python plot_figure_2.py")
sys.stderr.write("Done!")

sys.stderr.write("Plotting Figure 3 & supplements...\n")
os.system("python plot_figure_3.py")
sys.stderr.write("Done!")

sys.stderr.write("Plotting supplemental private marker fig...\n")
os.system("python plot_supplemental_private_marker.py")
sys.stderr.write("Done!")




