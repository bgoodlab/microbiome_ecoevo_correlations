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

output_file = open("figure_2_output.txt","w")
import bacterial_phylogeny_utils

genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

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
mpl.rcParams['legend.fontsize']  = 'small'

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


main_fig = plt.figure(figsize=(8.6,1.8))

outer_grid = gridspec.GridSpec(1,2, width_ratios=[2.6,4.4], wspace=0.25) 


right_grid = gridspec.GridSpecFromSubplotSpec(1,2,width_ratios=[1,1.2],subplot_spec=outer_grid[1],wspace=0.35)

right_right_grid = gridspec.GridSpecFromSubplotSpec(1,2,width_ratios=[1.4,1.1],subplot_spec=right_grid[1],wspace=0.1)

explained_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[1,1.8],subplot_spec=right_right_grid[0],hspace=0.2)

js_distribution_axis = plt.Subplot(main_fig, right_grid[0])
main_fig.add_subplot(js_distribution_axis)

bias_axis = plt.Subplot(main_fig, explained_grid[0])
main_fig.add_subplot(bias_axis)

explained_axis = plt.Subplot(main_fig, explained_grid[1])
main_fig.add_subplot(explained_axis)

extinct_axis = plt.Subplot(main_fig, right_right_grid[1])
main_fig.add_subplot(extinct_axis)


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

genus_abundance_matrix, genera = analysis.coarse_grain_abandance_matrix(abundance_matrix,speciess,samples,"genus")

family_abundance_matrix, families = analysis.coarse_grain_abandance_matrix(abundance_matrix,speciess,samples,"family")

phylum_abundance_matrix, phyla = analysis.coarse_grain_abandance_matrix(abundance_matrix,speciess,samples,"phylum")
    
sys.stderr.write("Performing preliminary calculations...\t")
within_host_events = analysis.calculate_within_host_events(within_host_changes)
    
host_classification_map = analysis.calculate_host_classifications(within_host_events)
host_records = host_classification_map.keys()
host_event_map = analysis.collate_within_host_events(within_host_events)
host_event_by_species_map = analysis.collate_within_host_events_by_species(within_host_events)
host_ecological_map = analysis.calculate_host_ecological_changes(host_records,abundance_matrix,speciess,samples)

# Calculate same thing at higher levels of taxonomic organization. 
host_genus_ecological_map = analysis.calculate_host_ecological_changes(host_records,genus_abundance_matrix,genera,samples)
host_family_ecological_map = analysis.calculate_host_ecological_changes(host_records,family_abundance_matrix,families,samples)
host_phylum_ecological_map = analysis.calculate_host_ecological_changes(host_records,phylum_abundance_matrix,phyla,samples)


sys.stderr.write("Done!\n")  

species_event_map = analysis.collate_species_events(within_host_events)


###
#
# Prepare new data structures for statistics 
# (will be numpy arrays w/ rows = # of timepoint comparisons)
#
###

jsds = [] # Jensen-Shannon (JS) distances between timepoints
profile_jsdss = [] # JS distances w/ different cutoffs on delta and favg
alpha0s = [] # shannon entropies at timepoint 0
dalphas = [] # difference in shannon entropy between timepoint 1 and 0
nextincts = [] # num species that go "extinct" between timepoint 0 and 1
ninvades = [] # num species that go "invade" between timepoint 0 and 1
neffs = [] # num species that "contribute" to the Jensen Shannon distance
delta_effs = [] # a measure of effect size contributing to JSD, <delta^2>/<delta> 
x_resolved_jsdss = []

genus_jsds = []
family_jsds = []
phylum_jsds = []

null_idxs = [] # idxs addressing hosts w/ no evolutionary change
modification_idxs = [] # w/ modification
replacement_idxs = [] # w/ replacement


# Between host versions of above (for visual comparison)
between_jsds = []
between_alpha0s = []
between_dalphas = []
between_nextincts = []
between_delta_effs= []
num_between_comparisons = 1000

sys.stderr.write("Collating calculations...\t")
for host_idx in xrange(0,len(host_records)):
    
    host_record = host_records[host_idx]
    
    alpha0,alpha1,js,profile_jsd_tuple,x_resolved_jsds_tuple, neff_delta, neff_avg, delta_eff,extinctions,invasions = host_ecological_map[host_record]
    
    profile_jsds, delta_mins, favg_mins = profile_jsd_tuple
    #x_resolved_jsds, x_bins = x_resolved_jsds_tuple
    dalpha = alpha1-alpha0
    nextinct = len(extinctions)
    ninvade = len(invasions)
    neff = neff_delta #*1.0/neff_avg
    
    jsds.append(js)
    profile_jsdss.append(profile_jsds)
    alpha0s.append(alpha0)
    dalphas.append(dalpha)
    neffs.append(neff)
    delta_effs.append(delta_eff)
    nextincts.append(nextinct)
    ninvades.append(ninvade)
    #x_resolved_jsdss.append(x_resolved_jsds)
    
    genus_jsds.append(host_genus_ecological_map[host_record][2])
    family_jsds.append(host_family_ecological_map[host_record][2])
    phylum_jsds.append(host_phylum_ecological_map[host_record][2])
    
    if host_classification_map[host_record]==0:
        # not a genetic change
        null_idxs.append(True)
        modification_idxs.append(False)
        replacement_idxs.append(False)
    elif host_classification_map[host_record]==1:
        # at least one modification event (no replacements)
        null_idxs.append(False)
        modification_idxs.append(True)
        replacement_idxs.append(False)
                 
    elif host_classification_map[host_record]==2:
        null_idxs.append(False)
        modification_idxs.append(False)
        replacement_idxs.append(True)
        
jsds = numpy.array(jsds)
profile_jsdss = numpy.array(profile_jsdss)
alpha0s = numpy.array(alpha0s)
dalphas = numpy.array(dalphas)
nextincts = numpy.array(nextincts)
ninvades = numpy.array(ninvades)
neffs = numpy.array(neffs)
delta_effs = numpy.array(delta_effs)

genus_jsds = numpy.array(genus_jsds)
family_jsds = numpy.array(family_jsds)
phylum_jsds = numpy.array(phylum_jsds)

#x_resolved_jsdss = numpy.array(x_resolved_jsdss)

null_idxs = numpy.array(null_idxs)
modification_idxs = numpy.array(modification_idxs)
replacement_idxs = numpy.array(replacement_idxs)

# Implement a screen
#null_idxs = (null_idxs)*(alpha0s>=(alpha0s[replacement_idxs].min()))
#modification_idxs = (modification_idxs*(jsds<=0.8))
#null_idxs = (null_idxs*(alpha0s<=3.1))
#replacement_idxs = (replacement_idxs*(alpha0s<=3.1))

sys.stderr.write("Done!\n")

###
#
# Do same calculations for samples from different hosts
#
###

sys.stderr.write("Calculating between-host distribution...\t")
# First construct fake host timepoint pairs
fake_host_classification_map = {}
shuffled_host_records = [host_record for host_record in host_records]
for i in xrange(0,num_between_comparisons):
    
    shuffle(shuffled_host_records)
    while shuffled_host_records[0][1]==shuffled_host_records[1][1]:
        shuffle(shuffled_host_records)
    
    fake_host_classification_map[(shuffled_host_records[0][0], shuffled_host_records[0][1],shuffled_host_records[0][2],shuffled_host_records[1][3])] = -1

# Calculate ecological changes for them
fake_host_records = fake_host_classification_map.keys()

fake_host_ecological_map = analysis.calculate_host_ecological_changes(fake_host_records, abundance_matrix, speciess,samples)   

for host_record in fake_host_records:
    
    cohort,subject,t0,t1 = host_record
    
    alpha0,alpha1,js,profile_jsds,dummy,neff_delta, neff_avg, delta_eff, extinctions,invasions = fake_host_ecological_map[host_record]
    dalpha=alpha1-alpha0
    nextinct = len(extinctions)
    ninvade = len(invasions)
    
    between_jsds.append(js)
    between_nextincts.append(nextinct)
    between_alpha0s.append(alpha0)
    between_dalphas.append(dalpha)
    between_delta_effs.append(delta_eff)
    
between_jsds = numpy.array(between_jsds)
between_alpha0s = numpy.array(between_alpha0s)
between_dalphas = numpy.array(between_dalphas)
between_nextincts = numpy.array(between_nextincts)
between_delta_effs = numpy.array(between_delta_effs)
sys.stderr.write("Done!\n")

output_file.write("Sample sizes: total=%d, none=%d, modification=%d, replacement=%d, fake_between=%d\n" % ( len(null_idxs), null_idxs.sum(),modification_idxs.sum(),replacement_idxs.sum(),len(between_jsds)))

#######
#
# Calculate KS statistics
#
######
sys.stderr.write("Calculating KS statistics...\t")

# Comparing JS distance distributions
observed_modification_js_ks, observed_modification_jsstar = analysis.ks_2samp_greater(jsds[modification_idxs], jsds[null_idxs])
observed_replacement_js_ks, observed_replacement_jsstar = analysis.ks_2samp_greater(jsds[replacement_idxs], jsds[null_idxs])

# Comparing Profile JS distance distributions
modification_profile_js_kss = numpy.zeros_like(profile_jsdss[0,0])
replacement_profile_js_kss = numpy.zeros_like(profile_jsdss[0,0])
for i in xrange(0,profile_jsdss[0].shape[1]):
    for j in xrange(0,profile_jsdss[0].shape[2]):
        ks, jsstar = analysis.ks_2samp_greater(profile_jsdss[modification_idxs,0,i,j], profile_jsdss[null_idxs,0, i,j])
        modification_profile_js_kss[i,j] = ks
        
        ks, jsstar = analysis.ks_2samp_greater(profile_jsdss[replacement_idxs,0,i,j], profile_jsdss[null_idxs,0,i,j])
        replacement_profile_js_kss[i,j] = ks
    
    
# Comparing entropy distributions
observed_modification_alpha0_ks, observed_modification_alpha0star = analysis.ks_2samp_greater(alpha0s[modification_idxs], alpha0s[null_idxs])
observed_replacement_alpha0_ks, observed_replacement_alpha0star = analysis.ks_2samp_greater(alpha0s[replacement_idxs], alpha0s[null_idxs])

# Comparing change in entropy distributions 
# (multiply by negative 1 because alternative hypothesis is that
#  evolutionary events lead to reduction in entropy)
observed_modification_dalpha_ks, observed_modification_dalphastar = analysis.ks_2samp_greater(-1*dalphas[modification_idxs], -1*dalphas[null_idxs])
observed_replacement_dalpha_ks, observed_replacement_dalphastar = analysis.ks_2samp_greater(-1*dalphas[replacement_idxs], -1*dalphas[null_idxs])
# Same but for absolute value
observed_modification_absdalpha_ks, observed_modification_absdalphastar = analysis.ks_2samp_greater(numpy.fabs(dalphas[modification_idxs]), numpy.fabs(dalphas[null_idxs]))
observed_replacement_absdalpha_ks, observed_replacement_absdalphastar = analysis.ks_2samp_greater(numpy.fabs(dalphas[replacement_idxs]), numpy.fabs(dalphas[null_idxs]))


# Comparing num extinct distributions
observed_modification_nextinct_ks,dummy = analysis.ks_2samp_greater(nextincts[modification_idxs], nextincts[null_idxs])
observed_modification_pextinct = (nextincts[modification_idxs]>0).sum()*1.0/modification_idxs.sum()

observed_replacement_nextinct_ks,dummy = analysis.ks_2samp_greater(nextincts[replacement_idxs], nextincts[null_idxs])
observed_replacement_pextinct = (nextincts[replacement_idxs]>0).sum()*1.0/replacement_idxs.sum()

# Same for invasions
observed_modification_ninvade_ks,dummy = analysis.ks_2samp_greater(ninvades[modification_idxs], ninvades[null_idxs])

observed_replacement_ninvade_ks,dummy = analysis.ks_2samp_greater(ninvades[replacement_idxs], ninvades[null_idxs])


# Comparing Coarse JS distance distributions
observed_modification_genus_js_ks, observed_modification_genus_jsstar = analysis.ks_2samp_greater(genus_jsds[modification_idxs], genus_jsds[null_idxs])
observed_replacement_genus_js_ks, observed_replacement_genus_jsstar = analysis.ks_2samp_greater(genus_jsds[replacement_idxs], genus_jsds[null_idxs])

observed_modification_family_js_ks, observed_modification_family_jsstar = analysis.ks_2samp_greater(family_jsds[modification_idxs], family_jsds[null_idxs])
observed_replacement_family_js_ks, observed_replacement_family_jsstar = analysis.ks_2samp_greater(family_jsds[replacement_idxs], family_jsds[null_idxs])

observed_modification_phylum_js_ks, observed_modification_phylum_jsstar = analysis.ks_2samp_greater(phylum_jsds[modification_idxs], phylum_jsds[null_idxs])
observed_replacement_phylum_js_ks, observed_replacement_phylum_jsstar = analysis.ks_2samp_greater(phylum_jsds[replacement_idxs], phylum_jsds[null_idxs])


sys.stderr.write("Done!\n")

######
#
# Make plots!
#
######
sys.stderr.write("Making plots...\n")

######
#
# Stacked bar charts w/ a few examples
#
######

# Used to plot specific examples
example_host_record_1 = ('hmp','160603188', '700034926', '700112755') # example of big 
example_host_record_2 = ('hmp', '706846339', '700106170', '700163772') # stercoris 
example_host_record_3 = ('hmp','884419868', '700113975', '700163811') # no change
example_host_records = [example_host_record_1, example_host_record_2, example_host_record_3]
example_host_classifications = ['Rep','Mod','None']

example_grid = gridspec.GridSpecFromSubplotSpec(1,len(example_host_records)+1,width_ratios=([2]*len(example_host_records)+[5]),wspace=0.3,subplot_spec=outer_grid[0])

legend_axis = plt.Subplot(main_fig, example_grid[-1])
main_fig.add_subplot(legend_axis)
legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])      



# Build a joint abundance matrix for the samples we're looking at
fs = []
for example_idx in xrange(0,len(example_host_records)):
    
    host_record = example_host_records[example_idx]

    example_fs = numpy.vstack([abundance_matrix[:,samples.index(host_record[2])], 
abundance_matrix[:,samples.index(host_record[3])]]).T
    fs.append(example_fs)
    
fs = numpy.array(fs)
fs = numpy.swapaxes(fs,0,1)

# Now do all the sorting and colormapping

species_abundances_by_family = {}
family_abundance_map = {}


for species_idx in xrange(0,len(speciess)):

    species = speciess[species_idx]
    family = bacterial_phylogeny_utils.get_family_name(species, genus_family_map=genus_family_map)
    
    # Only include things that hit a given size
    if not (fs[species_idx,:,:].flatten()>2e-02).any():
        continue

    #if species=='Bacteroides_stercoris_56735':
    #    print fs[species_idx]    
    # DEBUG GOOD TO HERE
    
    if family not in family_abundance_map:
        family_abundance_map[family] = 0
        species_abundances_by_family[family] = []
                
    family_abundance_map[family] += fs[species_idx,:,:].sum()
    species_abundances_by_family[family].append((species, fs[species_idx,:,0],fs[species_idx,:,1])) 
    
# sort families by abundance
# This will be a problem because dict is not sorted.
family_items = [(family,family_abundance_map[family]) for family in family_abundance_map]

sorted_family_items = list(sorted(family_items, key=operator.itemgetter(1),reverse=True))

sorted_species_idxs_by_family = {}
species_color_map = {}

for family_idx in xrange(0,len(sorted_family_items)):
    family = sorted_family_items[family_idx][0]
    
    colors = colorbrewer.get_colors(family_idx,len(species_abundances_by_family[family])+1)
    
    species_items = [(species_idx, species_abundances_by_family[family][species_idx][1].sum()+species_abundances_by_family[family][species_idx][2].sum()) for species_idx in xrange(0,len(species_abundances_by_family[family]))]
    
    sorted_species_items = list(sorted(species_items, key=operator.itemgetter(1),reverse=True))
    
    sorted_species_idxs_by_family[family] = []
    for species_idx,dummy in sorted_species_items:
        sorted_species_idxs_by_family[family].append(species_idx)
        color = colors[species_idx]
        species_name = species_abundances_by_family[family][species_idx][0]
        species_color_map[species_name] = color
        
        #if species_name=='Bacteroides_stercoris_56735':
        #    print species_abundances_by_family[family][species_idx][1]     
        # DEBUG GOOD TO HERE
        
        #legend_axis.plot([-2,-1], [-2,-1],'s',markersize=3,markeredgewidth=0.0,color=color,label=bacterial_phylogeny_utils.get_pretty_name(species_name))
    
    legend_axis.plot([-2,-1], [-2,-1],'s',markersize=3,markeredgewidth=0.0,color=colors[0],label=family)

legend_axis.plot([-2,-1], [-2,-1],'s', markersize=3,markeredgewidth=0.0,color='#F2C1FB',label='Other')

handles, labels = legend_axis.get_legend_handles_labels()
legend_axis.legend(handles[::-1], labels[::-1],loc=(-0.08,-0.05),frameon=False,fontsize=6,numpoints=1,ncol=1,handlelength=1,handletextpad=0.1)   

        
for example_idx in xrange(0,len(example_host_records)):

    #print "Example", example_idx+1
    host_record = example_host_records[example_idx]
    host_classification = example_host_classifications[example_idx]
    
    output_file.write("Host %d example: %s, classification=%s\n" % (example_idx+1,str(host_record),host_classification))
    
    example_axis = plt.Subplot(main_fig, example_grid[example_idx])
    main_fig.add_subplot(example_axis)
    example_axis.set_xlabel('Host %d\n(%s)' % (example_idx+1,host_classification))  
    xs = numpy.array([0,1,2,3])
    lowers = numpy.zeros_like(xs)*1.0
    uppers = numpy.array(lowers,copy=True)
    for family_idx in xrange(0,len(sorted_family_items)):
        
        family = sorted_family_items[family_idx][0]
    
        for species_idx in sorted_species_idxs_by_family[family]:
            species_name = species_abundances_by_family[family][species_idx][0]
            f0 = species_abundances_by_family[family][species_idx][1][example_idx]
            f1 = species_abundances_by_family[family][species_idx][2][example_idx]
            
            #if species_name=='Bacteroides_stercoris_56735':
            #    print f0,f1     
            # DEBUG GOOD TO HERE
            
            if f0<2e-02 and f1<2e-02:
                continue # don't waste space to plot if too small...
                
            #print species_name,f0,f1
            
            uppers = lowers+numpy.array([f0,f0,f1,f1])
            
            edgecolor='0.3'
            linewidth=0.25
            
            
            if species_name in host_event_by_species_map[host_record]:
                event = host_event_by_species_map[host_record][species_name]
                
                if event>0:
                
                    print species_name, event, uppers[0]-lowers[0],uppers[-1]-lowers[-1]
                    
                    if event==2:
                        event_color=replacement_color
                        #edgecolor = replacement_color
                        #linewidth=2
                    elif event==1:
                        event_color=modification_color
                        #edgecolor=modification_color
                        #linewidth=2
                
                    example_axis.text(0.43, (uppers[0]+lowers[0])/2,'$\\blacktriangleright$',color=event_color,va='center',fontsize=7,fontweight='bold')
                
            color = species_color_map[species_name]
            example_axis.fill_between(xs, lowers,uppers,facecolor=color,edgecolor=edgecolor,linewidth=linewidth)
            lowers = numpy.array(uppers,copy=True)
            #print lowers[0]
    uppers = numpy.ones_like(lowers)
    example_axis.fill_between(xs, lowers,uppers,facecolor='#F2C1FB',edgecolor='0.3',linewidth=0.25)

    example_axis.fill_between([xs[1],xs[2]],[0,0],[1,1],color='w',linewidth=0,alpha=0.5)
    example_axis.set_xlim([0.75,2.25])
    example_axis.set_ylim([0,1])

    example_axis.spines['top'].set_visible(False)
    example_axis.spines['right'].set_visible(False)
    example_axis.spines['left'].set_visible(False)
    example_axis.spines['bottom'].set_visible(False)

    example_axis.set_xticks([])
    example_axis.set_yticks([])
           

##### 
#
# Distributions of JS distances
#
####

js_distribution_axis.set_ylabel('Jensen-Shannon distance')
js_distribution_axis.set_ylim([0.1,0.9])

Xs = [0,-1,1,2]
for ys,colorVal,X,width,light_colorVal in zip([jsds[modification_idxs],jsds[null_idxs],jsds[replacement_idxs],between_jsds],[modification_color,'0.7',replacement_color,'0.3'],Xs,[0.3,0.3,0.1,0.3],['#CED2ED','#ECECEC','#E8C8CA','#D2D2D2']):

    kernel = gaussian_kde(ys)

    
    theory_ys = numpy.linspace(ys.min(),ys.max(),100)
    theory_pdf = kernel(theory_ys)
    
    scale = width/theory_pdf.max()
    
    
    xs = uniform(-1,1,size=len(ys))*kernel(ys)*scale

    q25 = numpy.quantile(ys,0.25)
    q50 = numpy.quantile(ys,0.5)
    q75 = numpy.quantile(ys,0.75)

    
    
    js_distribution_axis.fill_betweenx(theory_ys, X-theory_pdf*scale,X+theory_pdf*scale,linewidth=0.25,facecolor=light_colorVal,edgecolor=colorVal)
    
    
    
    #pylab.plot([X-scale*kernel(q25), X+scale*kernel(q25)],[q25,q25],'-',color=colorVal,linewidth=0.25)
    #pylab.plot([X-scale*kernel(q50), X+scale*kernel(q50)],[q50,q50],'-',color=colorVal,linewidth=1)
    #pylab.plot([X-scale*kernel(q75), X+scale*kernel(q75)],[q75,q75],'-',color=colorVal,linewidth=0.25)
    
    other_width = width+0.1
    js_distribution_axis.plot([X-other_width,X+other_width],[q25,q25],'-',color=colorVal,linewidth=1)
    js_distribution_axis.plot([X-other_width,X+other_width],[q50,q50],'-',color=colorVal,linewidth=1)
    js_distribution_axis.plot([X-other_width,X+other_width],[q75,q75],'-',color=colorVal,linewidth=1)
    js_distribution_axis.plot([X-other_width,X-other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    js_distribution_axis.plot([X+other_width,X+other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    
    if len(ys)<900:
        js_distribution_axis.plot(X+xs,ys,'.',color=colorVal,alpha=0.5,markersize=5,markeredgewidth=0.0)    
    #pylab.fill_between([X-0.05,X+0.05],[q25,q25],[q75,q75],facecolor='0.3')
    #pylab.fill_between([X-0.2,X+0.2], [q50-0.01,q50-0.01],[q50+0.01,q50+0.01],facecolor='0.3')
    #pylab.plot([X],[q50], '.',color='w',markersize=4,markeredgecolor=colorVal,markeredgewidth=0.25)

js_distribution_axis.set_ylim([0.1,1])
js_distribution_axis.plot([1.4,1.4],[0,1],'k:',linewidth=0.5)
js_distribution_axis.set_xticks([-1,0,1,2])
#js_distribution_axis.set_xticklabels(['No changes','Modification','Replacement','Between host'],rotation=45,ha='right')
js_distribution_axis.set_xticklabels(['None','Mod','Rep','Btwn'],rotation=90) #,ha='right')

##### 
#
# Extinction histogram plot
#
####

extinct_axis.set_xlabel('Extinctions')
#extinct_axis.set_ylabel('Fraction of hosts')

bins = numpy.arange(0,6)-0.5
bins[-1]=100
xs = numpy.arange(0,5)

for desired_idxs,color,offset in zip([null_idxs,modification_idxs,replacement_idxs],[null_color,modification_color,replacement_color],[-0.3,0,0.3]):

    num_zeros = (nextincts[desired_idxs]==0).sum()
    num_ones = (nextincts[desired_idxs]==1).sum()
    num_twogreater = (nextincts[desired_idxs]>1.5).sum()
    
    xs = numpy.array([0,1,2])
    hs = numpy.array([num_zeros,num_ones,num_twogreater])
    
    extinct_axis.bar(xs+offset,hs*1.0/hs.sum(),width=0.3,color=color)
    
extinct_axis.set_xticks([0,1,2])
extinct_axis.set_xticklabels(['0','1','2+'])
extinct_axis.set_yticks([])

##### 
#
# Supplemental extinction/invasion plot
#
####

extinction_fig = plt.figure(figsize=(5,2))
extinction_outer_grid = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.2) 

extinction_axis = plt.Subplot(extinction_fig, extinction_outer_grid[0])
extinction_fig.add_subplot(extinction_axis)

extinction_axis.set_xlabel('Extinctions')
extinction_axis.set_ylabel('Fraction of hosts')

invasion_axis = plt.Subplot(extinction_fig, extinction_outer_grid[1])
extinction_fig.add_subplot(invasion_axis)

invasion_axis.set_xlabel('Invasions')

for desired_idxs,color,offset in zip([null_idxs,modification_idxs,replacement_idxs],[null_color,modification_color,replacement_color],[-0.3,0,0.3]):

    for ns,axis in zip([nextincts,ninvades],[extinction_axis,invasion_axis]):
        
        xs = numpy.arange(0,5)
        hs = numpy.array([(ns[desired_idxs]==k).sum() for k in xs])
        
        hs[-1] = (ns[desired_idxs]>=xs[-1]).sum()
        
        hs = hs*1.0/hs.sum()
        
        #xs = xs[1:]
        #hs = hs[1:]
    
        axis.bar(xs+offset,hs,width=0.3,color=color)
    
extinction_axis.set_xticks([0,1,2,3,4])
extinction_axis.set_xticklabels(['0','1','2','3','4'])
extinction_axis.set_ylim([0,1])
invasion_axis.set_xticks([0,1,2,3,4])
invasion_axis.set_xticklabels(['0','1','2','3','4'])
invasion_axis.set_ylim([0,1])
invasion_axis.set_yticklabels([])

extinction_fig.savefig('supplemental_extinction_invasion.pdf',bbox_inches='tight')

##### 
#
# Plot KS values directly on survival plots
#
####
survival_fig = plt.figure(figsize=(4,1.25))
survival_outer_grid = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.25) 

modification_survival_axis = plt.Subplot(survival_fig, survival_outer_grid[0])
survival_fig.add_subplot(modification_survival_axis)

xs,data_survival,null_survival = analysis.ks_survivals(jsds[modification_idxs],jsds[null_idxs])

modification_survival_axis.plot(xs,null_survival,'-',color='0.7')
modification_survival_axis.plot(xs,data_survival,'-',color=modification_color)
modification_survival_axis.plot([observed_modification_jsstar, observed_modification_jsstar],[0,1],'k:')
modification_survival_axis.set_ylim([0,1])
modification_survival_axis.set_xlim([0,1])

modification_survival_axis.set_ylabel('Fraction hosts')
modification_survival_axis.set_xlabel('                                                 Min Jensen-Shannon Distance')
replacement_survival_axis = plt.Subplot(survival_fig, survival_outer_grid[1])
survival_fig.add_subplot(replacement_survival_axis)

xs,data_survival,null_survival = analysis.ks_survivals(jsds[replacement_idxs],jsds[null_idxs])

replacement_survival_axis.plot(xs,null_survival,'-',color='0.7')
replacement_survival_axis.plot(xs,data_survival,'-',color=replacement_color)
replacement_survival_axis.plot([observed_replacement_jsstar, observed_replacement_jsstar],[0,1],'k:')

replacement_survival_axis.set_ylim([0,1])
replacement_survival_axis.set_xlim([0,1])
replacement_survival_axis.set_yticklabels([])

survival_fig.savefig('supplemental_ks.pdf',bbox_inches='tight')

##### 
#
# Coarsened Distributions of JS distances
#
####

coarse_fig = plt.figure(figsize=(6,2))

coarse_outer_grid = gridspec.GridSpec(1,3, width_ratios=[1,1,1], wspace=0.25) 

for coarse_idx, coarse_type, coarse_jsds in zip(numpy.arange(0,3), ['Genus','Family','Phylum'],[genus_jsds,family_jsds,phylum_jsds]):
    
    coarse_axis = plt.Subplot(coarse_fig, coarse_outer_grid[coarse_idx])
    coarse_fig.add_subplot(coarse_axis)

    coarse_axis.set_title(coarse_type,fontsize=7)
    if coarse_idx==0:
        coarse_axis.set_ylabel('Jensen-Shannon distance')

    coarse_axis.set_ylim([0,1])

    Xs = [0,-1,1]
    for ys,colorVal,X,width,light_colorVal in zip([coarse_jsds[modification_idxs],coarse_jsds[null_idxs],coarse_jsds[replacement_idxs]],[modification_color,'0.7',replacement_color],Xs,[0.3,0.3,0.1],['#CED2ED','#ECECEC','#E8C8CA']):

        kernel = gaussian_kde(ys)

    
        theory_ys = numpy.linspace(ys.min(),ys.max(),100)
        theory_pdf = kernel(theory_ys)
    
        scale = width/theory_pdf.max()
    
    
        xs = uniform(-1,1,size=len(ys))*kernel(ys)*scale

        q25 = numpy.quantile(ys,0.25)
        q50 = numpy.quantile(ys,0.5)
        q75 = numpy.quantile(ys,0.75)
   
        coarse_axis.fill_betweenx(theory_ys, X-theory_pdf*scale,X+theory_pdf*scale,linewidth=0.25,facecolor=light_colorVal,edgecolor=colorVal)
    
        other_width = width+0.1
        coarse_axis.plot([X-other_width,X+other_width], [q25,q25],'-',color=colorVal,linewidth=1)
        coarse_axis.plot([X-other_width,X+other_width], [q50,q50],'-',color=colorVal,linewidth=1)
        coarse_axis.plot([X-other_width,X+other_width], [q75,q75],'-',color=colorVal,linewidth=1)
        coarse_axis.plot([X-other_width,X-other_width], [q25,q75],'-',color=colorVal,linewidth=1)
        coarse_axis.plot([X+other_width,X+other_width], [q25,q75],'-',color=colorVal,linewidth=1)
    
        if len(ys)<900:
            coarse_axis.plot(X+xs, ys,'.',color=colorVal,alpha=0.5,markersize=5,markeredgewidth=0.0)    

        coarse_axis.set_ylim([0,1])
        coarse_axis.set_xticks([-1,0,1])
        coarse_axis.set_xticklabels(['None','Mod','Rep'],rotation=90) 

coarse_fig.savefig('supplemental_coarse_jsd.pdf',bbox_inches='tight')
##### 
#
# Distributions of shannon entropies
#
####
pylab.figure(figsize=(2,2))

Xs = [0,-1,1,2]
for ys,colorVal,X,width,light_colorVal,type in zip([alpha0s[modification_idxs],alpha0s[null_idxs],alpha0s[replacement_idxs],between_alpha0s],[modification_color,'0.7',replacement_color,'0.3'],Xs,[0.3,0.3,0.1,0.3],['#CED2ED','#ECECEC','#E8C8CA','#D2D2D2'],['modification','null','replacement','between']):

    kernel = gaussian_kde(ys)

    theory_ys = numpy.linspace(ys.min(),ys.max(),100)
    theory_pdf = kernel(theory_ys)
    
    scale = width/theory_pdf.max()
    
    xs = uniform(-1,1,size=len(ys))*kernel(ys)*scale

    q25 = numpy.quantile(ys,0.25)
    q50 = numpy.quantile(ys,0.5)
    q75 = numpy.quantile(ys,0.75)

    num_targets = ((ys>=q25)*(ys<=q75)).sum()
    num_nulls = ((alpha0s[null_idxs]>=q25)*(alpha0s[null_idxs]<=q75)).sum()
    
    print type, "Median rate =", num_targets*1.0/(num_targets+num_nulls), "All rate = ", len(ys)*1.0/(len(ys)+null_idxs.sum()), num_targets, len(ys), num_nulls, null_idxs.sum()
    
    pylab.fill_betweenx(theory_ys, X-theory_pdf*scale,X+theory_pdf*scale,linewidth=0.25,facecolor=light_colorVal,edgecolor=colorVal)
    
    other_width = width+0.1
    pylab.plot([X-other_width,X+other_width],[q25,q25],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X+other_width],[q50,q50],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X+other_width],[q75,q75],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X-other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    pylab.plot([X+other_width,X+other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    
    if len(ys)<900:
        pylab.plot(X+xs,ys,'.',color=colorVal,alpha=0.5,markersize=5,markeredgewidth=0.0)    
    
pylab.plot([1.4,1.4],[0,1],'k:',linewidth=0.5)
pylab.xticks([-1,0,1,2])
pylab.gca().set_xticklabels([])
#pylab.savefig('evolution_alpha0_correlation.pdf',bbox_inches='tight')

##### 
#
# Distributions of delta_effs
#
####
pylab.figure(figsize=(2,2))

Xs = [0,-1,1,2]
for ys,colorVal,X,width,light_colorVal in zip([delta_effs[modification_idxs],delta_effs[null_idxs],delta_effs[replacement_idxs],between_delta_effs],[modification_color,'0.7',replacement_color,'0.3'],Xs,[0.3,0.3,0.1,0.3],['#CED2ED','#ECECEC','#E8C8CA','#D2D2D2']):

    kernel = gaussian_kde(ys)
   
    theory_ys = numpy.linspace(ys.min(),ys.max(),100)
    theory_pdf = kernel(theory_ys)
    
    scale = width/theory_pdf.max()
    
    xs = uniform(-1,1,size=len(ys))*kernel(ys)*scale

    q25 = numpy.quantile(ys,0.25)
    q50 = numpy.quantile(ys,0.5)
    q75 = numpy.quantile(ys,0.75)
 
    
    pylab.fill_betweenx(theory_ys, X-theory_pdf*scale,X+theory_pdf*scale,linewidth=0.25,facecolor=light_colorVal,edgecolor=colorVal)
       
    other_width = width+0.1
    pylab.plot([X-other_width,X+other_width],[q25,q25],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X+other_width],[q50,q50],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X+other_width],[q75,q75],'-',color=colorVal,linewidth=1)
    pylab.plot([X-other_width,X-other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    pylab.plot([X+other_width,X+other_width],[q25,q75],'-',color=colorVal,linewidth=1)
    
    if len(ys)<900:
        pylab.plot(X+xs,ys,'.',color=colorVal,alpha=0.5,markersize=5,markeredgewidth=0.0)    
    
pylab.plot([1.4,1.4],[0,1],'k:',linewidth=0.5)
pylab.xticks([-1,0,1,2])
pylab.gca().set_xticklabels([])
#pylab.savefig('evolution_deltaeff_correlation.pdf',bbox_inches='tight')


##### 
#
# Distributions of change in Shannon entropy
#
####

dalpha_fig = plt.figure(figsize=(5.25,2))

dalpha_outer_grid = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.35) 

dalpha_axis = plt.Subplot(dalpha_fig, dalpha_outer_grid[0])
dalpha_fig.add_subplot(dalpha_axis)
dalpha_axis.set_ylabel('$\\Delta$(Shannon Diversity)')

absdalpha_axis = plt.Subplot(dalpha_fig, dalpha_outer_grid[1])
dalpha_fig.add_subplot(absdalpha_axis)
absdalpha_axis.set_ylabel('|$\\Delta$(Shannon Diversity)|')

Xs = [0,-1,1,2]
for raw_ys,colorVal,X,width,light_colorVal,type in zip([dalphas[modification_idxs],dalphas[null_idxs],dalphas[replacement_idxs],between_dalphas],[modification_color,'0.7',replacement_color,'0.3'],Xs,[0.3,0.3,0.1,0.3],['#CED2ED','#ECECEC','#E8C8CA','#D2D2D2'],['modification','null','replacement','between']):

    for ys,axis in zip([raw_ys,numpy.fabs(raw_ys)],[dalpha_axis, absdalpha_axis]):
    
        kernel = gaussian_kde(ys)

        theory_ys = numpy.linspace(ys.min(),ys.max(),100)
        theory_pdf = kernel(theory_ys)
    
        scale = width/theory_pdf.max()
    
        xs = uniform(-1,1,size=len(ys))*kernel(ys)*scale

        q25 = numpy.quantile(ys,0.25)
        q50 = numpy.quantile(ys,0.5)
        q75 = numpy.quantile(ys,0.75)

        num_targets = ((ys>=q25)*(ys<=q75)).sum()
        num_nulls = ((alpha0s[null_idxs]>=q25)*(alpha0s[null_idxs]<=q75)).sum()
    
        axis.fill_betweenx(theory_ys, X-theory_pdf*scale,X+theory_pdf*scale,linewidth=0.25,facecolor=light_colorVal,edgecolor=colorVal)
    
        other_width = width+0.1
        axis.plot([X-other_width,X+other_width], [q25,q25],'-',color=colorVal,linewidth=1)
        axis.plot([X-other_width,X+other_width], [q50,q50],'-',color=colorVal,linewidth=1)
        axis.plot([X-other_width,X+other_width], [q75,q75],'-',color=colorVal,linewidth=1)
        axis.plot([X-other_width,X-other_width], [q25,q75],'-',color=colorVal,linewidth=1)
        axis.plot([X+other_width,X+other_width], [q25,q75],'-',color=colorVal,linewidth=1)
    
        if len(ys)<900:
            axis.plot(X+xs, ys,'.',color=colorVal,alpha=0.5,markersize=5,markeredgewidth=0.0)    
    
dalpha_axis.plot([1.4,1.4],[-3,3],'k:',linewidth=0.5)
absdalpha_axis.plot([1.4,1.4],[0,3],'k:',linewidth=0.5)

dalpha_axis.set_xticks([-1,0,1,2])
dalpha_axis.set_xticklabels(['None','Mod','Rep','Btwn'],rotation=90) 

absdalpha_axis.set_xticks([-1,0,1,2])
absdalpha_axis.set_xticklabels(['None','Mod','Rep','Btwn'],rotation=90) 

dalpha_fig.savefig('supplemental_dalpha.pdf',bbox_inches='tight')

##### 
#
# JSD as function of alpha0
#
####
pylab.figure(figsize=(2,2))

Xs = [0,-1,1,2]
for ys,xs,colorVal,X,width,light_colorVal,zorder,label in zip([jsds[modification_idxs],jsds[null_idxs],jsds[replacement_idxs]],[alpha0s[modification_idxs],alpha0s[null_idxs],alpha0s[replacement_idxs]],[modification_color,'0.7',replacement_color,'0.3'],Xs,[0.3,0.3,0.1,0.3],['#CED2ED','#ECECEC','#E8C8CA','#D2D2D2'],[1,0,2],['Mod','None','Rep']):

    slope, intercept, rvalue, pvalue, stderr = linregress(xs,ys)
    
    print label, slope, pvalue
    
    theory_xs = numpy.linspace(xs.min(),xs.max(),10)
    
    pylab.plot(xs,ys,'.',color=colorVal,alpha=0.5,zorder=zorder)
    pylab.plot(theory_xs,intercept+slope*theory_xs,'-',color=colorVal,zorder=zorder)
    
    
    
#pylab.savefig('supplemental_js_alpha0_correlation.pdf',bbox_inches='tight')

##### 
#
### Calculate relative probability of one event or other as function of alpha0
#
####
pylab.figure(figsize=(2,2))

alphas = numpy.linspace(1,3.5,50)

total_null = null_idxs.sum()*1.0
total_modification = modification_idxs.sum()*1.0
total_replacement = replacement_idxs.sum()*1.0

modification_pavg = total_modification/total_null
replacement_pavg = total_replacement/total_null

null_alpha0s = alpha0s[null_idxs]
null_greaters = (null_alpha0s[:,None]>alphas[None,:]).sum(axis=0)

modification_alpha0s = alpha0s[modification_idxs]
modification_greaters = (modification_alpha0s[:,None]>alphas[None,:]).sum(axis=0)

replacement_alpha0s = alpha0s[replacement_idxs]
replacement_greaters = (replacement_alpha0s[:,None]>alphas[None,:]).sum(axis=0)

modification_ps = (modification_greaters+modification_pavg)/(modification_greaters+null_greaters+1.0)

replacement_ps = (replacement_greaters+replacement_pavg)/(replacement_greaters+null_greaters+1.0)

pylab.plot(alphas,modification_ps,'-',color=modification_color)
pylab.plot(alphas,replacement_ps,'-',color=replacement_color)

#pylab.savefig('alpha_event_correlation.pdf',bbox_inches='tight')


##### 
#
### Fraction of JSD as function of delta_min 
#
####

power=2
fraction_idx = len(profile_jsdss[0,:,0])/2

bias_axis.plot([1,1e03],100*numpy.array([0.5,0.5]),'k:',linewidth=0.5)
for desired_idxs,color in zip([null_idxs,modification_idxs,replacement_idxs],[null_color,modification_color,replacement_color]):

    #fractions = numpy.power(profile_jsdss[desired_idxs,0,:,0],power).sum(axis=0)
    #positive_fractions = numpy.power(profile_jsdss[desired_idxs,1,:,0],power).sum(axis=0)
    #negative_fractions = numpy.power(profile_jsdss[desired_idxs,2,:,0],power).sum(axis=0)

    #positive_fractions = positive_fractions/fractions[0]
    #negative_fractions = negative_fractions/fractions[0]
    #fractions = fractions/fractions[0]
    
    fractions = numpy.power(profile_jsdss[desired_idxs,0,:,0],power)
    positive_fractions = numpy.power(profile_jsdss[desired_idxs,1,:,0],power)
    negative_fractions = numpy.power(profile_jsdss[desired_idxs,2,:,0],power)
    
    total_fractions = numpy.mean(fractions/(fractions[:,0])[:,None],axis=0)
    positive_fractions = numpy.mean((positive_fractions/(fractions[:,0])[:,None]),axis=0)
    negative_fractions = numpy.mean((negative_fractions/(fractions[:,0])[:,None]),axis=0)
    
    print positive_fractions/total_fractions
    
    xs = numpy.power(10,delta_mins)
    explained_axis.semilogx(xs,100*total_fractions,'-',color=color)
    bias_axis.semilogx(xs,100*positive_fractions/total_fractions,'-',color=color)
    
    #explained_axis.semilogx(xs[1:], positive_fractions[1:],'-',color=color)
    #explained_axis.semilogx(1.0/xs[1:], negative_fractions[1:],'-',color=color)

####
#
# Note to self: do weighting by x's for last figure
#
####
#xs = numpy.power(10,delta_mins)
#explained_axis.plot([xs[1],xs[1]],[0,1],'k:',linewidth=0.5)  
#explained_axis.plot([1.0/xs[1],1.0/xs[1]],[0,1],'k:',linewidth=0.5)    
explained_axis.set_ylim([0,110])
explained_axis.set_xlim([1,2e02])
explained_axis.set_ylabel('% JSD explained')
explained_axis.set_xlabel('Min fold change')

bias_axis.set_ylabel('% positive')
bias_axis.set_ylim([35,70])
bias_axis.set_xlim([1,2e02])
bias_axis.set_xticklabels([])
main_fig.savefig('figure_2.pdf',bbox_inches='tight')

sys.stderr.write("Done!\n")

#####
#
# Calculate bootstrapped pvalues!
#
####
sys.stderr.write("\nBootstrapping...\n")

# Finally do one first bootstrap test to look at symmetry of extinctions and invasions in hosts with no changes

null_extincts = (nextincts[null_idxs]>0)
null_invades = (ninvades[null_idxs]>0)

observed_num_invade = null_invades.sum()
observed_num_extinct = null_extincts.sum()
reverse_time = binomial(1,0.5,size=(num_bootstraps,len(null_extincts)))

bootstrapped_num_invades = (null_invades[None,:]*(1-reverse_time) + null_extincts[None,:]*reverse_time).sum(axis=1)
bootstrapped_num_extincts = (null_invades[None,:]*reverse_time + null_extincts[None,:]*(1-reverse_time)).sum(axis=1)

null_invasion_pvalue = ((numpy.fabs(bootstrapped_num_invades-bootstrapped_num_extincts)>=numpy.fabs(observed_num_invade-observed_num_extinct)).sum()+1.0)/(num_bootstraps+1.0)


sys.stderr.write("Creating bootstrapped within-host events...\t")
bootstrapped_within_host_eventss = analysis.calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=num_bootstraps,bootstrap_type='all')
sys.stderr.write("Done!\n")

sys.stderr.write("Collating bootstrap data...\n")

bootstrapped_modification_js_kss = []
bootstrapped_replacement_js_kss = []

bootstrapped_modification_nextinct_kss = []
bootstrapped_modification_pextincts = []

bootstrapped_replacement_nextinct_kss = []

bootstrapped_modification_ninvade_kss = []
bootstrapped_replacement_ninvade_kss = []



bootstrapped_modification_alpha0_kss = []
bootstrapped_modification_dalpha_kss = []
bootstrapped_modification_absdalpha_kss = []

bootstrapped_replacement_alpha0_kss = []
bootstrapped_replacement_dalpha_kss = []
bootstrapped_replacement_absdalpha_kss = []

bootstrapped_modification_profile_js_ksss = []
bootstrapped_replacement_profile_js_ksss = []

bootstrapped_modification_genus_js_kss = []
bootstrapped_replacement_genus_js_kss = []

bootstrapped_modification_family_js_kss = []
bootstrapped_replacement_family_js_kss = []

bootstrapped_modification_phylum_js_kss = []
bootstrapped_replacement_phylum_js_kss = []

for idx in xrange(0,num_bootstraps):
    if idx%100==0:
        print idx
        
    bootstrapped_host_classification_map = analysis.calculate_host_classifications(bootstrapped_within_host_eventss[idx])
    
    bootstrapped_null_idxs = []
    bootstrapped_modification_idxs = []
    bootstrapped_replacement_idxs = []
    
    for host_idx in xrange(0,len(host_records)):
    
        host_record = host_records[host_idx]
    
        if bootstrapped_host_classification_map[host_record]==0:
            # not a genetic change
            bootstrapped_null_idxs.append(True)
            bootstrapped_modification_idxs.append(False)
            bootstrapped_replacement_idxs.append(False)
        elif bootstrapped_host_classification_map[host_record]==1:
            # at least one modification event (no replacements)
            bootstrapped_null_idxs.append(False)
            bootstrapped_modification_idxs.append(True)
            bootstrapped_replacement_idxs.append(False)
                 
        elif bootstrapped_host_classification_map[host_record]==2:
            bootstrapped_null_idxs.append(False)
            bootstrapped_modification_idxs.append(False)
            bootstrapped_replacement_idxs.append(True)
    
    bootstrapped_null_idxs = numpy.array(bootstrapped_null_idxs)
    bootstrapped_modification_idxs = numpy.array(bootstrapped_modification_idxs)
    bootstrapped_replacement_idxs = numpy.array(bootstrapped_replacement_idxs)
    
    # Calculate KS statistics
    
    # Comparing JS distance distributions
    bootstrapped_modification_js_ks, bootstrapped_modification_jsstar = analysis.ks_2samp_greater(jsds[bootstrapped_modification_idxs], jsds[bootstrapped_null_idxs])
    bootstrapped_replacement_js_ks, bootstrapped_replacement_jsstar = analysis.ks_2samp_greater(jsds[bootstrapped_replacement_idxs], jsds[bootstrapped_null_idxs])

    # Comparing Profile JS distance distributions
    bootstrapped_modification_profile_js_kss = numpy.zeros_like(profile_jsdss[0,0])
    bootstrapped_replacement_profile_js_kss = numpy.zeros_like(profile_jsdss[0,0])
    for i in xrange(0,profile_jsdss[0].shape[1]):
        for j in xrange(0,profile_jsdss[0].shape[2]):
            ks, jsstar = analysis.ks_2samp_greater(profile_jsdss[bootstrapped_modification_idxs,0,i,j], profile_jsdss[bootstrapped_null_idxs, 0, i,j])
            bootstrapped_modification_profile_js_kss[i,j] = ks
        
            ks, jsstar = analysis.ks_2samp_greater(profile_jsdss[bootstrapped_replacement_idxs,0,i,j], profile_jsdss[bootstrapped_null_idxs,0, i,j])
            bootstrapped_replacement_profile_js_kss[i,j] = ks
    
    
    # Comparing entropy distributions
    bootstrapped_modification_alpha0_ks, bootstrapped_modification_alpha0star = analysis.ks_2samp_greater(alpha0s[bootstrapped_modification_idxs], alpha0s[bootstrapped_null_idxs])
    bootstrapped_replacement_alpha0_ks, bootstrapped_replacement_alpha0star = analysis.ks_2samp_greater(alpha0s[bootstrapped_replacement_idxs], alpha0s[bootstrapped_null_idxs])

    # Comparing change in entropy distributions
    bootstrapped_modification_dalpha_ks, bootstrapped_modification_dalphastar = analysis.ks_2samp_greater(-1*dalphas[bootstrapped_modification_idxs], -1*dalphas[bootstrapped_null_idxs])
    bootstrapped_replacement_dalpha_ks, bootstrapped_replacement_dalphastar = analysis.ks_2samp_greater(-1*dalphas[bootstrapped_replacement_idxs], -1*dalphas[bootstrapped_null_idxs])
    # Same for abs(delta entropy)
    bootstrapped_modification_absdalpha_ks, bootstrapped_modification_absdalphastar = analysis.ks_2samp_greater(numpy.fabs(dalphas[bootstrapped_modification_idxs]), numpy.fabs(dalphas[bootstrapped_null_idxs]))
    bootstrapped_replacement_absdalpha_ks, bootstrapped_replacement_absdalphastar = analysis.ks_2samp_greater(numpy.fabs(dalphas[bootstrapped_replacement_idxs]), numpy.fabs(dalphas[bootstrapped_null_idxs]))


    # Comparing num extinct distributions
    bootstrapped_modification_nextinct_ks,dummy = analysis.ks_2samp_greater(nextincts[bootstrapped_modification_idxs], nextincts[bootstrapped_null_idxs])
    bootstrapped_modification_pextinct = (nextincts[bootstrapped_modification_idxs]>0).sum()*1.0/bootstrapped_modification_idxs.sum()

    bootstrapped_replacement_nextinct_ks,dummy = analysis.ks_2samp_greater(nextincts[bootstrapped_replacement_idxs], nextincts[bootstrapped_null_idxs])
    bootstrapped_replacement_pextinct = (nextincts[bootstrapped_replacement_idxs]>0).sum()*1.0/bootstrapped_replacement_idxs.sum()
    
    bootstrapped_modification_js_kss.append(bootstrapped_modification_js_ks)
    bootstrapped_replacement_js_kss.append(bootstrapped_replacement_js_ks)
    bootstrapped_modification_nextinct_kss.append(bootstrapped_modification_nextinct_ks)
    bootstrapped_replacement_nextinct_kss.append(bootstrapped_replacement_nextinct_ks)
    bootstrapped_modification_pextincts.append(bootstrapped_modification_pextinct)

    # Now do it for num invade distributions
    bootstrapped_modification_ninvade_ks,dummy = analysis.ks_2samp_greater(ninvades[bootstrapped_modification_idxs], ninvades[bootstrapped_null_idxs])
    
    bootstrapped_replacement_ninvade_ks,dummy = analysis.ks_2samp_greater(ninvades[bootstrapped_replacement_idxs], ninvades[bootstrapped_null_idxs])
     
    bootstrapped_modification_ninvade_kss.append(bootstrapped_modification_ninvade_ks)
    bootstrapped_replacement_ninvade_kss.append(bootstrapped_replacement_ninvade_ks)

    
    bootstrapped_modification_alpha0_kss.append(bootstrapped_modification_alpha0_ks)
    bootstrapped_replacement_alpha0_kss.append(bootstrapped_replacement_alpha0_ks)
    
    bootstrapped_modification_dalpha_kss.append(bootstrapped_modification_dalpha_ks)
    bootstrapped_replacement_dalpha_kss.append(bootstrapped_replacement_dalpha_ks)
    
    bootstrapped_modification_absdalpha_kss.append(bootstrapped_modification_absdalpha_ks)
    bootstrapped_replacement_absdalpha_kss.append(bootstrapped_replacement_absdalpha_ks)
    
    bootstrapped_modification_profile_js_ksss.append( bootstrapped_modification_profile_js_kss)
    
    bootstrapped_replacement_profile_js_ksss.append( bootstrapped_replacement_profile_js_kss)
    
    
    # Do JS distributions for higher levels of taxonomic organization 
    # (didn't do this at first, but for completeness)
    bootstrapped_modification_js_ks, bootstrapped_modification_jsstar = analysis.ks_2samp_greater(genus_jsds[bootstrapped_modification_idxs], genus_jsds[bootstrapped_null_idxs])
    bootstrapped_replacement_js_ks, bootstrapped_replacement_jsstar = analysis.ks_2samp_greater(genus_jsds[bootstrapped_replacement_idxs], genus_jsds[bootstrapped_null_idxs])
    bootstrapped_modification_genus_js_kss.append(bootstrapped_modification_js_ks)
    bootstrapped_replacement_genus_js_kss.append(bootstrapped_replacement_js_ks)
    
    bootstrapped_modification_js_ks, bootstrapped_modification_jsstar = analysis.ks_2samp_greater(family_jsds[bootstrapped_modification_idxs], family_jsds[bootstrapped_null_idxs])
    bootstrapped_replacement_js_ks, bootstrapped_replacement_jsstar = analysis.ks_2samp_greater(family_jsds[bootstrapped_replacement_idxs], family_jsds[bootstrapped_null_idxs])
    bootstrapped_modification_family_js_kss.append(bootstrapped_modification_js_ks)
    bootstrapped_replacement_family_js_kss.append(bootstrapped_replacement_js_ks)
    
    bootstrapped_modification_js_ks, bootstrapped_modification_jsstar = analysis.ks_2samp_greater(phylum_jsds[bootstrapped_modification_idxs], phylum_jsds[bootstrapped_null_idxs])
    bootstrapped_replacement_js_ks, bootstrapped_replacement_jsstar = analysis.ks_2samp_greater(phylum_jsds[bootstrapped_replacement_idxs], phylum_jsds[bootstrapped_null_idxs])
    bootstrapped_modification_phylum_js_kss.append(bootstrapped_modification_js_ks)
    bootstrapped_replacement_phylum_js_kss.append(bootstrapped_replacement_js_ks)
    
print "Done!"

bootstrapped_modification_js_kss = numpy.array(bootstrapped_modification_js_kss)
bootstrapped_replacement_js_kss = numpy.array(bootstrapped_replacement_js_kss)
bootstrapped_modification_alpha0_kss = numpy.array(bootstrapped_modification_alpha0_kss)
bootstrapped_replacement_alpha0_kss = numpy.array(bootstrapped_replacement_alpha0_kss)
bootstrapped_modification_dalpha_kss = numpy.array(bootstrapped_modification_dalpha_kss)
bootstrapped_replacement_dalpha_kss = numpy.array(bootstrapped_replacement_dalpha_kss)
bootstrapped_modification_nextinct_kss = numpy.array( bootstrapped_modification_nextinct_kss)
bootstrapped_replacement_nextinct_kss = numpy.array(bootstrapped_replacement_nextinct_kss)
bootstrapped_modification_pextincts = numpy.array(bootstrapped_modification_pextincts)

bootstrapped_modification_profile_js_ksss = numpy.array(bootstrapped_modification_profile_js_ksss)
bootstrapped_replacement_profile_js_ksss = numpy.array(bootstrapped_replacement_profile_js_ksss)

# Coarse versions
bootstrapped_modification_genus_js_kss = numpy.array(bootstrapped_modification_genus_js_kss)
bootstrapped_replacement_genus_js_kss = numpy.array(bootstrapped_replacement_genus_js_kss)

bootstrapped_modification_family_js_kss = numpy.array(bootstrapped_modification_family_js_kss)
bootstrapped_replacement_family_js_kss = numpy.array(bootstrapped_replacement_family_js_kss)

bootstrapped_modification_phylum_js_kss = numpy.array(bootstrapped_modification_phylum_js_kss)
bootstrapped_replacement_phylumjs_kss = numpy.array(bootstrapped_replacement_phylum_js_kss)



output_file.write("Modification JSD KS P-value = %g, Observed JS*=%g, Observed fraction above = %g, null fraction above = %g \n" % (((bootstrapped_modification_js_kss>=observed_modification_js_ks).sum()+1.0)/(num_bootstraps+1.0), observed_modification_jsstar, (jsds[modification_idxs]>=observed_modification_jsstar).mean(),(jsds[null_idxs]>=observed_modification_jsstar).mean()))

output_file.write("Replacement JSD KS P-value = %g, Observed JS*=%g, Observed fraction above = %g, Null fraction above = %g \n" % (((bootstrapped_replacement_js_kss>=observed_replacement_js_ks).sum()+1.0)/(num_bootstraps+1.0), observed_replacement_jsstar, (jsds[replacement_idxs]>=observed_replacement_jsstar).mean(),(jsds[null_idxs]>=observed_replacement_jsstar).mean()))

output_file.write("Modification num_extinction KS P-value = %g\n" % (((bootstrapped_modification_nextinct_kss>=observed_modification_nextinct_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement num_extinction KS P-value = %g\n" % (((bootstrapped_replacement_nextinct_kss>=observed_replacement_nextinct_ks).sum()+1.0)/(num_bootstraps+1.0)))  

output_file.write("Modification num_invasion KS P-value = %g\n" % (((bootstrapped_modification_ninvade_kss>=observed_modification_ninvade_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement num_invasion KS P-value = %g\n" % (((bootstrapped_replacement_ninvade_kss>=observed_replacement_ninvade_ks).sum()+1.0)/(num_bootstraps+1.0))) 

output_file.write("Null invasion / replacement symmetry P-value = %g (%d, %d)\n" % (null_invasion_pvalue, observed_num_invade, observed_num_extinct)) 


output_file.write("Coarse versions of JSD comparisons:\n")

output_file.write("Genus Modification JSD KS P-value = %g\n" % (((bootstrapped_modification_genus_js_kss>=observed_modification_genus_js_ks).sum()+1.0)/(num_bootstraps+1.0)))
output_file.write("Genus Replacement JSD KS P-value = %g\n" % (((bootstrapped_replacement_genus_js_kss>=observed_replacement_genus_js_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Family Modification JSD KS P-value = %g\n" % (((bootstrapped_modification_family_js_kss>=observed_modification_family_js_ks).sum()+1.0)/(num_bootstraps+1.0)))
output_file.write("Family Replacement JSD KS P-value = %g\n" % (((bootstrapped_replacement_family_js_kss>=observed_replacement_family_js_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Phylum Modification JSD KS P-value = %g\n" % (((bootstrapped_modification_phylum_js_kss>=observed_modification_phylum_js_ks).sum()+1.0)/(num_bootstraps+1.0)))
output_file.write("Phylum Replacement JSD KS P-value = %g\n" % (((bootstrapped_replacement_phylum_js_kss>=observed_replacement_phylum_js_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Modification alpha0 KS P-value = %g\n" % (((bootstrapped_modification_alpha0_kss>=observed_modification_alpha0_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement alpha0 KS P-value = %g\n" % (((bootstrapped_replacement_alpha0_kss>=observed_replacement_alpha0_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Modification dalpha KS P-value = %g\n" % (((bootstrapped_modification_dalpha_kss>=observed_modification_dalpha_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement dalpha KS P-value = %g\n" % (((bootstrapped_replacement_dalpha_kss>=observed_replacement_dalpha_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Modification |dalpha| KS P-value = %g\n" % (((bootstrapped_modification_absdalpha_kss>=observed_modification_absdalpha_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement |dalpha| KS P-value = %g\n" % (((bootstrapped_replacement_absdalpha_kss>=observed_replacement_absdalpha_ks).sum()+1.0)/(num_bootstraps+1.0)))
     