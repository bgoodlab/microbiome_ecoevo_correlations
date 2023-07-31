import os
import sys

os.system("mkdir -p figures")
os.system("mkdir -p supplemental_data")

sys.stderr.write("Checking that data is there...\n")
os.system("python parse_data.py")
sys.stderr.write("Done!\n")

sys.stderr.write("Plotting Figure 1 & supplements...\n")
os.system("python plot_figure_1.py")
sys.stderr.write("Done!\n")

sys.stderr.write("Plotting Figure 2 & supplements...\n")
os.system("python plot_figure_2.py")
sys.stderr.write("Done!\n")

sys.stderr.write("Plotting Figure 3 & supplements...\n")
os.system("python plot_figure_3.py")
sys.stderr.write("Done!\n")

sys.stderr.write("Plotting supplemental private marker fig...\n")
os.system("python plot_supplemental_private_marker.py")
sys.stderr.write("Done!\n")

sys.stderr.write("Making Supplemental Tables...\n")
os.system("cp HMP1-2_ids_order.txt supplemental_data/TableS1.txt")
os.system("cp within_host_changes_consecutive.txt supplemental_data/TableS2.txt")
sys.stderr.write("Done!\n")