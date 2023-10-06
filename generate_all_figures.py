import os
import sys
import config

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

sys.stderr.write("Making Supplementary Data Files...\n")
os.system("cp %s supplemental_data/Supplementary_Data_1.txt" % config.default_metadata_filename)
os.system("cp %s supplemental_data/Supplementary_Data_2.txt" % config.default_within_host_filename)
os.system("gzcat within_host_changes_all_genes.txt.gz > supplemental_data/Supplementary_Data_3.txt")
os.system("cp %s supplemental_data/Supplementary_Data_4.txt" % config.default_coverage_filename)

sys.stderr.write("Done!\n")