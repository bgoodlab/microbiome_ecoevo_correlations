import numpy
import pylab
import pylab as plt
import parse_data
import analysis
import config
from numpy.random import choice,shuffle,normal,uniform,binomial
from scipy.stats import gaussian_kde, linregress
import sys
import colorbrewer
import operator

#####
#
#     Set up plotting stuff 
#
#####
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
#mpl.rcParams['legend.fontsize']  = 'small'

maxX=1
cNorm  = colors.Normalize(vmin=-maxX, vmax=maxX)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pylab.get_cmap('coolwarm') )

null_color = scalarMap.to_rgba(0)
modification_color = scalarMap.to_rgba(-1)
replacement_color = scalarMap.to_rgba(1)

#####
#
# Set up figure
#
#########


main_fig = plt.figure(figsize=(2.5,1.7))
outer_grid = gridspec.GridSpec(1,1)


axis = plt.Subplot(main_fig, outer_grid[0])
main_fig.add_subplot(axis)

axis.set_xlabel('Num SNV differences')
axis.set_ylabel('Fraction private\nmarkers preserved')

#####
#
#  Load data and perform initial calculations 
#
#####

num_bootstraps = config.default_num_bootstraps

sys.stderr.write("Loading data...\t")
sample_metadata_map = parse_data.parse_sample_metadata_map()
within_host_changes = parse_data.parse_within_host_changes()
abundance_matrix,speciess,samples = parse_data.parse_abundances()
sys.stderr.write("Done!\n")

for record_idx in range(0,len(within_host_changes)):
    
    cohort,subject,t0,t1,species,Lsnp,ksnp,Lprivate,kprivate = within_host_changes[record_idx]

    if ksnp < 0.5:
        continue
        
    if Lprivate<9.5:
        continue
            
    fraction_preserved = 1-kprivate*1.0/Lprivate
    
    if ksnp > config.default_replacement_threshold:
        color = replacement_color
    else:
        color = modification_color
        
    axis.semilogx(ksnp, fraction_preserved,'o',markersize=4,color=color,alpha=0.5)
    
axis.plot([1,1],[2,3],'-',color=replacement_color,label='Rep')
axis.plot([1,1],[2,3],'-',color=modification_color, label='Mod')
axis.set_ylim([-0.05,1.05])
axis.legend(loc='upper right',frameon=False)

main_fig.savefig('supplemental_private_marker.pdf',bbox_inches='tight')