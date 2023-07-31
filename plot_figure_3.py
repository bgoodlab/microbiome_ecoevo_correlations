import numpy
import pylab
import parse_data
import analysis
from numpy.random import choice,shuffle,normal,uniform
from scipy.stats import ks_2samp
from scipy.stats import gaussian_kde, binom, fisher_exact
import pylab as plt
import sys
from scipy.special import xlogy
from math import log,log10

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

output_file = open("figures/figure_3_output.txt","w")

import bacterial_phylogeny_utils
genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

def min_taxonomic_distance(x,ys):
	return min(bacterial_phylogeny_utils.calculate_taxonomic_distances_from_focal_species(x, ys,genus_family_map,genus_phylum_map))
	
	4 - (species_taxon_map['phylum'][x]==species_taxon_map['phylum'][y]) - (species_taxon_map['family'][x]==species_taxon_map['family'][y]) - (species_taxon_map['genus'][x]==species_taxon_map['genus'][y])-(x==y)


plotted_taxonomic_level = 'family'

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'
mpl.rcParams['xtick.major.pad']='2'
mpl.rcParams['ytick.major.pad']='2'

maxX=1
cNorm  = colors.Normalize(vmin=-maxX, vmax=maxX)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pylab.get_cmap('coolwarm') )

null_color = scalarMap.to_rgba(0)
modification_color = scalarMap.to_rgba(-1)
replacement_color = scalarMap.to_rgba(1)

# DEBUG
num_bootstraps = 10000
#num_bootstraps = 10000

#####
#
# Set up figure
#
#########


main_fig = plt.figure(figsize=(7,2))

outer_grid = gridspec.GridSpec(1,4, width_ratios=[1.1,1.35,1.0,0.8], wspace=0.6) 

left_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[1.8,1],subplot_spec=outer_grid[0],hspace=0.25)

bias_axis = plt.Subplot(main_fig, left_grid[1])
main_fig.add_subplot(bias_axis)
bias_axis.set_ylabel('% positive')

cdf_axis = plt.Subplot(main_fig, left_grid[0])
main_fig.add_subplot(cdf_axis)
cdf_axis.set_ylabel('Fraction species')
bias_axis.set_xlabel('Min fold change')

middle_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[1,1.8],subplot_spec=outer_grid[1],hspace=0.25)

explained_axis = plt.Subplot(main_fig, middle_grid[0])
main_fig.add_subplot(explained_axis)
explained_axis.set_ylabel('% explained')

js_axis = plt.Subplot(main_fig, middle_grid[1])
main_fig.add_subplot(js_axis)
js_axis.set_ylabel("Jensen-Shannon")
js_axis.set_xlabel("Fold change \n of focal species")

right_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[1,1],subplot_spec=outer_grid[2],hspace=0.2)

species_axis = plt.Subplot(main_fig, right_grid[0])
main_fig.add_subplot(species_axis)

genus_axis = plt.Subplot(main_fig, right_grid[1])
main_fig.add_subplot(genus_axis)
genus_axis.set_ylabel('                           Fraction hosts')
genus_axis.set_xlabel('Min % explained')

supplemental_taxonomic_fig = plt.figure(figsize=(2,5))
taxonomic_grid = taxonomic_outer_grid = gridspec.GridSpec(1,1)
taxonomic_grid = gridspec.GridSpecFromSubplotSpec(4,1,height_ratios=[1,1,1,1],subplot_spec=taxonomic_outer_grid[0],hspace=0.2) 

taxonomic_species_axis = plt.Subplot(supplemental_taxonomic_fig, taxonomic_grid[0])
supplemental_taxonomic_fig.add_subplot(taxonomic_species_axis)

levels = ['genus','family','phylum']
taxon_axes = {}
for i in range(0,len(levels)):
	level = levels[i]
	taxon_axis = plt.Subplot(supplemental_taxonomic_fig, taxonomic_grid[i+1])
	supplemental_taxonomic_fig.add_subplot(taxon_axis)
	taxon_axes[level] = taxon_axis
	
taxon_axes['family'].set_ylabel('                       Fraction of hosts')
taxon_axes['phylum'].set_xlabel('Min % explained by\nspecies in focal taxa')

#other_fig = plt.figure(figsize=(1.5,2))
#other_grid = gridspec.GridSpec(1,1)
extinction_grid = gridspec.GridSpecFromSubplotSpec(3,1,height_ratios=[0.3,1,1],subplot_spec=outer_grid[3],hspace=0.15) 

legend_axis = plt.Subplot(main_fig, extinction_grid[0])
main_fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])

modification_extinction_axis = plt.Subplot(main_fig, extinction_grid[1])
main_fig.add_subplot(modification_extinction_axis)

replacement_extinction_axis = plt.Subplot(main_fig, extinction_grid[2])
main_fig.add_subplot(replacement_extinction_axis)
replacement_extinction_axis.set_ylabel('                     Number of extinctions')

sample_metadata_map = parse_data.parse_sample_metadata_map() 
within_host_changes = parse_data.parse_within_host_changes()
abundance_matrix,speciess,samples = parse_data.parse_abundances()
# Array version for fast indexing
speciess_array = numpy.array(speciess)	
within_host_events = analysis.calculate_within_host_events(within_host_changes)

host_event_map = analysis.collate_within_host_events(within_host_events)

host_records = list(host_event_map.keys())

host_classification_map = analysis.calculate_host_classifications(within_host_events)
host_ecological_map = analysis.calculate_host_ecological_changes(host_classification_map,abundance_matrix,speciess,samples) 
js_item_map = analysis.calculate_js_item_map(host_classification_map,abundance_matrix,speciess,samples)	   

species_event_map = analysis.collate_species_events(within_host_events)

# Calculate species-species distance matrices for all species in big abundance matrix
#taxonomic_distance_matrix = bacterial_phylogeny_utils.calculate_taxonomic_distance_matrix(speciess,genus_family_map,genus_phylum_map)

null_fs = []
null_jss = []
modification_fs = []
modification_jss = []
replacement_fs = []
replacement_jss = []

species_species_map = {}
species_taxon_map = {'genus':{},'family':{},'phylum':{}}


######################################################
#
# Loop over all quasi-phaseable species/host combos
# and calculate fold change of focal species
#
######################################################
for record_idx in range(0,len(within_host_events)):
	
	cohort,subject,t0,t1,species,event = within_host_events[record_idx]
		
	host_record = (cohort,subject,t0,t1)

	if species not in species_species_map:
		
		species_idxs = bacterial_phylogeny_utils.get_idxs_from_same_taxon(species,speciess,'species',genus_family_map,genus_phylum_map)
		species_species_map[species] = species_idxs
		
		for level in species_taxon_map:	
			species_idxs = bacterial_phylogeny_utils.get_idxs_from_same_taxon(species,speciess,level,genus_family_map,genus_phylum_map)
			species_taxon_map[level][species] = species_idxs
		
	# get f0 and f1
	f0 = analysis.get_frequency(abundance_matrix,speciess,samples,species,t0)
	f1 = analysis.get_frequency(abundance_matrix,speciess,samples,species,t1)
	
	js = host_ecological_map[host_record][2]
	#print species,t0,t1,f0,f1,event
	
	# Purposely excludes null events in replacement/modification hosts
	# as well as modifications in replacement hosts
	if host_classification_map[host_record]==0 and event==0:
		null_fs.append((f0,f1))
		null_jss.append(js)
	elif host_classification_map[host_record]==1 and event==1:
		modification_fs.append((f0,f1))
		modification_jss.append(js)
	elif host_classification_map[host_record]==2 and event==2:
		replacement_fs.append((f0,f1))
		replacement_jss.append(js)
		
null_fs = numpy.array(null_fs)
modification_fs = numpy.array(modification_fs)
replacement_fs = numpy.array(replacement_fs)

null_jss = numpy.array(null_jss)
modification_jss = numpy.array(modification_jss)
replacement_jss = numpy.array(replacement_jss)

null_fs = numpy.clip(null_fs,1e-05,1-1e-05)
replacement_fs = numpy.clip(replacement_fs,1e-05,1-1e-05)
modification_fs = numpy.clip(modification_fs,1e-05,1-1e-05)

def calculate_logits(fs):
	#return numpy.log(fs[:,1]/(1-fs[:,1])*(1-fs[:,0])/fs[:,0])
	return fs[:,1]/fs[:,0]
	#return fs[:,1]/(1-fs[:,1])*(1-fs[:,0])/fs[:,0]
	
def calculate_logfoldchanges(fs):
	return numpy.log(fs[:,1])-numpy.log(fs[:,0])
	
def calculate_abs_logfoldchange(rs):
	return numpy.exp(numpy.abs(numpy.log(rs)))
	
def calculate_jsd_contributions(fs):
	
	favgs = (fs[:,1]+fs[:,0])/2.0
	dfs = fs[:,1]-fs[:,0]
	
	xs = dfs/2/favgs
	
	deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
	
	return deltas*favgs
	
null_logits = calculate_logits(null_fs)
modification_logits = calculate_logits(modification_fs)
replacement_logits = calculate_logits(replacement_fs)

null_contributions = calculate_jsd_contributions(null_fs)/numpy.square(null_jss)
modification_contributions = calculate_jsd_contributions(modification_fs)/numpy.square(modification_jss)
replacement_contributions = calculate_jsd_contributions(replacement_fs)/numpy.square(replacement_jss)

#########################################################################
#
# Loop over hosts to calculate contribution from focal species and focal genus/family/phylum
#
#########################################################################
species_ps = []
taxon_ps = {'genus':[],'family':[],'phylum':[]}
modification_idxs = []
replacement_idxs = []


bootstrapped_species_pss = []
bootstrapped_taxon_pss = {'genus':[],'family':[],'phylum':[]}

modification_extinction_distances = []
replacement_extinction_distances = []

bootstrapped_modification_extinction_distances = []
bootstrapped_replacement_extinction_distances = []

for host_idx in range(0,len(host_records)):
	
	host_record = host_records[host_idx]
	
	js_items = js_item_map[host_record]
	
	extinctions = host_ecological_map[host_record][-2]
	
	# Only do stuff within replacement and modification events
	if host_classification_map[host_record]>0:
		
		focal_species_idxs = numpy.zeros_like(js_items)
		taxon_idxs = {level: numpy.zeros_like(js_items) for level in taxon_ps}
		# TODO: LEFT OFF HERE!
		
		host_species_weights = []
		host_species_taxon_idxs = {level: [] for level in taxon_ps}
		host_species_species_idxs = []
		host_num_events = 0
		
		for species,event in host_event_map[host_record]:
			
			if event==host_classification_map[host_record]:
				
				focal_species_idxs += species_species_map[species]
				for level in taxon_idxs:
					taxon_idxs[level] += species_taxon_map[level][species]
				host_num_events += 1
				
			# Create list of species to draw future events from...	  
			host_species_weights.append( species_event_map[species][event]*1.0/sum(species_event_map[species]))
			host_species_species_idxs.append(species_species_map[species])
			for level in host_species_taxon_idxs:
				host_species_taxon_idxs[level].append(species_taxon_map[level][species])
				
		species_p = js_items[focal_species_idxs>0].sum()/js_items.sum()
		species_ps.append(species_p)
		
		for level in taxon_ps:
			taxon_p = js_items[taxon_idxs[level]>0].sum()/js_items.sum()
			taxon_ps[level].append(taxon_p)		   
		
		taxonomic_distances = []
		if len(extinctions)>0:
			focal_species_list = list(speciess_array[focal_species_idxs>0])
			for extinct_species in extinctions:
				taxonomic_distances.append( min_taxonomic_distance(extinct_species,focal_species_list))
			 
		host_species_idxs = numpy.arange(0,len(host_species_weights))
		host_species_weights = numpy.array(host_species_weights)
		host_species_weights = host_species_weights/host_species_weights.sum()
		host_species_species_idxs = numpy.array(host_species_species_idxs)
		for level in host_species_taxon_idxs:
			host_species_taxon_idxs[level] = numpy.array(host_species_taxon_idxs[level])
		
		if host_classification_map[host_record]==1:
			modification_idxs.append(True)
			replacement_idxs.append(False)
			modification_extinction_distances.extend(taxonomic_distances)
		else:
			modification_idxs.append(False)
			replacement_idxs.append(True)
			replacement_extinction_distances.extend(taxonomic_distances)
		
		# Now do bootstrapped version
		bootstrapped_taxon_ps = {level: [] for level in taxon_ps}
		bootstrapped_species_ps = []
		for bootstrap_idx in range(0,num_bootstraps):
			
			bootstrapped_event_idxs = choice(host_species_idxs,size=host_num_events,replace=False,p=host_species_weights)
			
			bootstrapped_species_idxs = host_species_species_idxs[bootstrapped_event_idxs,:].sum(axis=0)
			
			bootstrapped_species_p = js_items[bootstrapped_species_idxs>0].sum()/js_items.sum()
			
			bootstrapped_species_ps.append(bootstrapped_species_p)
			
			for level in bootstrapped_taxon_ps:			  
				bootstrapped_taxon_idxs = host_species_taxon_idxs[level][bootstrapped_event_idxs,:].sum(axis=0)
			
				bootstrapped_taxon_p = js_items[bootstrapped_taxon_idxs>0].sum()/js_items.sum()
			
				bootstrapped_taxon_ps[level].append(bootstrapped_taxon_p)
			
			
			if len(extinctions)>0:
				bootstrapped_taxonomic_distances = []
				focal_species_list = list(speciess_array[bootstrapped_species_idxs>0])
				for extinct_species in extinctions: 
					bootstrapped_taxonomic_distances.append( min_taxonomic_distance( extinct_species,focal_species_list))
		
				if host_classification_map[host_record]==1:
					bootstrapped_modification_extinction_distances.extend( bootstrapped_taxonomic_distances)
				else:
					bootstrapped_replacement_extinction_distances.extend( bootstrapped_taxonomic_distances)
				
		bootstrapped_species_pss.append(bootstrapped_species_ps)	   
		for level in bootstrapped_taxon_pss:
			bootstrapped_taxon_pss[level].append(bootstrapped_taxon_ps[level])		 
		
replacement_idxs = numpy.array(replacement_idxs)
modification_idxs = numpy.array(modification_idxs)
species_ps = numpy.array(species_ps)
bootstrapped_species_pss = numpy.array(bootstrapped_species_pss)
for level in taxon_ps:
	taxon_ps[level] = numpy.array(taxon_ps[level])
	bootstrapped_taxon_pss[level] = numpy.array(bootstrapped_taxon_pss[level])

##############
#
# Plot figures
#
#############


theory_rs = numpy.logspace(0,1,20)

bias_axis.semilogx(theory_rs,50*numpy.ones_like(theory_rs),'k:',linewidth=0.5)
for foldchanges,color,label in zip([null_logits, modification_logits, replacement_logits],[null_color, modification_color, replacement_color],['None','Mod','Rep']):
	
	fraction_positives = (foldchanges[:,None]>=theory_rs[None,:]).sum(axis=0)*1.0/len(foldchanges)
	
	fraction_negatives = (foldchanges[:,None]<=1.0/theory_rs[None,:]).sum(axis=0)*1.0/len(foldchanges)
	
	fraction_total = fraction_positives+fraction_negatives
	cdf_axis.semilogx(theory_rs, fraction_total,'-',color=color,label=label)
	
	bias_axis.semilogx(theory_rs[fraction_total>0], 100*fraction_positives[fraction_total>0]/fraction_total[fraction_total>0],'-',color=color)
	
	#print label, fraction_positives[fraction_total>0]/fraction_total[fraction_total>0]
cdf_axis.legend(loc='upper right',frameon=False)
cdf_axis.set_xlim([1,10])
bias_axis.set_xlim([1,10])
bias_axis.set_ylim([25,75])
bias_axis.set_xticklabels([],minor=True)
#bias_axis.set_xticklabels([],minor=False)
cdf_axis.set_xticklabels([],minor=True)
cdf_axis.set_xticklabels([],minor=False)



#####
#
# Fold change vs scatter plot
#
####
js_axis.semilogx([1,1],[0,1],'k:',linewidth=0.5)
js_axis.plot(null_logits,null_jss,'.',color='0.7',markersize=5,markeredgewidth=0.5,alpha=0.5)
js_axis.plot(modification_logits,modification_jss,'.',color=modification_color,markersize=7,markeredgewidth=0.5,alpha=0.5)
js_axis.plot(replacement_logits,replacement_jss,'.',color=replacement_color,markersize=7,markeredgewidth=0.5,alpha=0.5)

js_axis.set_ylim([0,1])
js_axis.set_yticks([0,0.25,0.5,0.75,1])
js_axis.set_yticklabels(['0','.25','.50','.75','1'])
js_axis.set_xlim([3e-02,3e01])

explained_axis.semilogx([1,1],[0,100],'k:',linewidth=0.5)
explained_axis.plot(null_logits,100*null_contributions,'.',color='0.7',markersize=5,markeredgewidth=0.5,alpha=0.5)
explained_axis.plot(modification_logits,100*modification_contributions,'.',color=modification_color,markersize=7,markeredgewidth=0.5,alpha=0.5)
explained_axis.plot(replacement_logits,100*replacement_contributions,'.',color=replacement_color,markersize=7,markeredgewidth=0.5,alpha=0.5)

explained_axis.set_ylim([0,60])
explained_axis.set_xlim([3e-02,3e01])
explained_axis.set_xticklabels([])

observed_fraction_positive = (modification_logits>1).sum()*1.0/len(modification_logits) 
observed_null_fraction_positive = (null_logits>1).sum()*1.0/len(null_logits) 
pvalue = binom.sf((modification_logits>1).sum(), len(modification_logits) , 0.5)

output_file.write("%d positive fold changes in %d modification events (%g, p-value=%g)\n" % ((modification_logits>1).sum(), len(modification_logits), observed_fraction_positive, pvalue))

observed_modification_df_ks, observed_modification_dfstar = analysis.ks_2samp_greater(modification_logits, null_logits)

observed_replacement_df_ks, observed_replacement_dfstar = analysis.ks_2samp_greater(replacement_logits, null_logits)

observed_modification_absdf_ks, observed_modification_dfstar = analysis.ks_2samp_greater(numpy.fabs(numpy.log(modification_logits)), numpy.fabs(numpy.log(null_logits)))

observed_replacement_absdf_ks, observed_replacement_dfstar = analysis.ks_2samp_greater(numpy.fabs(numpy.log(replacement_logits)), numpy.fabs(numpy.log(null_logits)))

#print observed_modification_df_ks, observed_modification_dfstar

######
#
# Plot fraction explained by focal
#
######

bootstrapped_ps = bootstrapped_species_pss.flatten()
species_axis.plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'k:',linewidth=0.5)

bootstrapped_ps = bootstrapped_species_pss[modification_idxs].flatten()
observed_mean = species_ps[modification_idxs].mean()
bootstrapped_means = []
for bootstrap_idx in range(0,num_bootstraps):
	mean = bootstrapped_species_pss[modification_idxs,bootstrap_idx].mean() 
	bootstrapped_means.append(mean)	  
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("modification focal species percent explained mean pvalue = %g\n" % pvalue)

bootstrapped_ps = bootstrapped_species_pss[replacement_idxs].flatten()
observed_mean = species_ps[replacement_idxs].mean()
bootstrapped_means = []
for bootstrap_idx in range(0,num_bootstraps):
	mean = bootstrapped_species_pss[replacement_idxs,bootstrap_idx].mean() 
	bootstrapped_means.append(mean)	  
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("replacement focal species percent explained mean pvalue = %g\n" % pvalue)


for level in bootstrapped_taxon_pss:
	bootstrapped_ps = bootstrapped_taxon_pss[level][modification_idxs].flatten()
	observed_mean = taxon_ps[level][modification_idxs].mean()
	bootstrapped_means = []
	for bootstrap_idx in range(0,num_bootstraps):
		mean = bootstrapped_taxon_pss[level][modification_idxs,bootstrap_idx].mean() 
		bootstrapped_means.append(mean)	  
	bootstrapped_means = numpy.array(bootstrapped_means)
	pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0) * 1.0/(len(bootstrapped_means)+1.0)
	effect_size = (observed_mean-bootstrapped_means.mean())
	output_file.write("modification focal %s percent explained mean pvalue = %g, effect size = %g\n" % (level, pvalue, effect_size))

	bootstrapped_ps = bootstrapped_taxon_pss[level][replacement_idxs].flatten()
	observed_mean = taxon_ps[level][replacement_idxs].mean()
	bootstrapped_means = []
	for bootstrap_idx in range(0,num_bootstraps):
		mean = bootstrapped_taxon_pss[level][replacement_idxs,bootstrap_idx].mean() 
		bootstrapped_means.append(mean)	  
	bootstrapped_means = numpy.array(bootstrapped_means)
	pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0) * 1.0/(len(bootstrapped_means)+1.0)
	effect_size = (observed_mean-bootstrapped_means.mean())
	output_file.write("replacement focal %s percent explained mean pvalue = %g, effect size=%g\n" % (level, pvalue, effect_size))

species_axis.plot(100*numpy.sort(species_ps[modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[modification_idxs]))),'-',color=modification_color)

species_axis.plot(100*numpy.sort(species_ps[replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[replacement_idxs]))),'-',color=replacement_color)

species_axis.set_xlim([0,50])

level = plotted_taxonomic_level
bootstrapped_ps = bootstrapped_taxon_pss[level].flatten()
genus_axis.plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'k:',linewidth=0.5)

genus_axis.plot(100*numpy.sort(taxon_ps[level][modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(taxon_ps[level][modification_idxs]))),'-',color=modification_color)

genus_axis.plot(100*numpy.sort(taxon_ps[level][replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(taxon_ps[level][replacement_idxs]))),'-',color=replacement_color)

#bootstrapped_ps = bootstrapped_taxon_pss[replacement_idxs].flatten()
#genus_axis.plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'-',color=replacement_color,alpha=0.5)


genus_axis.set_xlim([0,50])
genus_axis.set_xticks([0,25,50])

species_axis.set_xticklabels([])

species_axis.plot([200],[1],'k-',label='Observed')
species_axis.plot([200],[1],'k:',linewidth=0.5,label='Expected')
species_axis.legend(loc='upper right',frameon=False)

######
#
# Plot extinction part
#
######

ks = numpy.arange(1,5)

counts = numpy.histogram(modification_extinction_distances, bins=(numpy.arange(0,5)+0.5))[0]

null_counts = numpy.histogram(bootstrapped_modification_extinction_distances, bins=(numpy.arange(0,5)+0.5))[0]
null_counts = null_counts*(counts.sum()/null_counts.sum())

modification_extinction_axis.bar(ks-0.4,counts,width=0.4,align='edge',color=modification_color)
modification_extinction_axis.bar(ks,null_counts,width=0.4,color='0.7',align='edge')

counts = numpy.histogram(replacement_extinction_distances, bins=(numpy.arange(0,5)+0.5))[0]

null_counts = numpy.histogram(bootstrapped_replacement_extinction_distances, bins=(numpy.arange(0,5)+0.5))[0]
null_counts = null_counts*(counts.sum()/null_counts.sum())

replacement_extinction_axis.bar(ks-0.4,counts,width=0.4,align='edge',color=replacement_color)
replacement_extinction_axis.bar(ks,null_counts,width=0.4,color='0.7',align='edge')

modification_extinction_axis.set_xlim([0.5,4.5])
replacement_extinction_axis.set_xlim([0.5,4.5])
modification_extinction_axis.set_xticks([1,2,3,4])
modification_extinction_axis.set_xticklabels([])

replacement_extinction_axis.set_xticks([1,2,3,4])
replacement_extinction_axis.set_xticklabels(['genus','family','phylum','domain'],rotation=45,rotation_mode='anchor',fontsize=7,ha='right')
#other_fig.savefig('supplemental_extinction_distances.pdf',bbox_inches='tight')

legend_axis.plot([0],[-5],'s',color='k',label='Observed')
legend_axis.plot([0],[-5],'s',color='0.7',label='Expected')
legend_axis.legend(loc='center',frameon=False,numpoints=1,markerscale=0.5)

main_fig.savefig('figures/figure_3.pdf',bbox_inches='tight')


######
#
# Plot supplemental fraction explained for different taxonomic distances
#
######

taxonomic_species_axis.plot(100*numpy.sort(species_ps[modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[modification_idxs]))),'-',color=modification_color)

taxonomic_species_axis.plot(100*numpy.sort(species_ps[replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[replacement_idxs]))),'-',color=replacement_color)

taxonomic_species_axis.set_xlim([0,100])

for level in taxon_axes:
	
	bootstrapped_ps = bootstrapped_taxon_pss[level].flatten()
	taxon_axes[level].plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'k:',linewidth=0.5)

	taxon_axes[level].plot(100*numpy.sort(taxon_ps[level][modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(taxon_ps[level][modification_idxs]))),'-',color=modification_color)

	taxon_axes[level].plot(100*numpy.sort(taxon_ps[level][replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(taxon_ps[level][replacement_idxs]))),'-',color=replacement_color)

	taxon_axes[level].set_xlim([0,100])

	taxon_axes[level].text(95,0.85,level,ha='right')
	if level!='phylum':
		taxon_axes[level].set_xticklabels([])
	
taxon_axes['phylum'].set_xticks([0,25,50,75,100])

taxonomic_species_axis.text(95,0.85,'species',ha='right')
taxonomic_species_axis.set_xticklabels([])

taxon_axes['phylum'].plot([200],[1],'k-',label='Observed')
taxon_axes['phylum'].plot([200],[1],'k:',linewidth=0.5,label='Expected')
taxon_axes['phylum'].legend(loc='lower left',frameon=False)

supplemental_taxonomic_fig.savefig('figures/supplemental_explained.pdf',bbox_inches='tight')


######
#
# Plot fraction explained vs phylogenetic distance
#
######



# Bootstrapping!
print("Bootstrapping...")
bootstrapped_within_host_eventss = analysis.calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=num_bootstraps,bootstrap_type='genus')
print("Processing bootstraps!")
bootstrapped_modification_df_kss = []
bootstrapped_replacement_df_kss = []
bootstrapped_modification_absdf_kss = []
bootstrapped_replacement_absdf_kss = []

bootstrapped_fraction_positives = []

for idx in range(0,num_bootstraps):	

	if idx%100==0:
		print(idx)
	
	within_host_events = bootstrapped_within_host_eventss[idx]
	host_classification_map = analysis.calculate_host_classifications(within_host_events)

	null_fs = []
	modification_fs = []
	replacement_fs = []

	for record_idx in range(0,len(within_host_events)):
	
		cohort,subject,t0,t1,species,event = within_host_events[record_idx]
	
		host_record = (cohort,subject,t0,t1)
	
		# get f0 and f1
		f0 = analysis.get_frequency(abundance_matrix,speciess,samples,species,t0)
		f1 = analysis.get_frequency(abundance_matrix,speciess,samples,species,t1)
	
		if host_classification_map[host_record]==0 and event==0:
			null_fs.append((f0,f1))
		elif host_classification_map[host_record]==1 and event==1:
			modification_fs.append((f0,f1))
		elif host_classification_map[host_record]==2 and event==2:
			replacement_fs.append((f0,f1))
		
	null_fs = numpy.array(null_fs)
	modification_fs = numpy.array(modification_fs)
	replacement_fs = numpy.array(replacement_fs)

	null_fs = numpy.clip(null_fs,1e-05,1-1e-05)
	replacement_fs = numpy.clip(replacement_fs,1e-05,1-1e-05)
	modification_fs = numpy.clip(modification_fs,1e-05,1-1e-05)

	null_logits = calculate_logits(null_fs)
	modification_logits = calculate_logits(modification_fs)
	replacement_logits = calculate_logits(replacement_fs)
	
	abs_null_logits = numpy.fabs(numpy.log(null_logits))
	abs_modification_logits = numpy.fabs(numpy.log(modification_logits))
	abs_replacement_logits = numpy.fabs(numpy.log(replacement_logits))
	

	bootstrapped_fraction_positives.append( (modification_logits>1).sum()*1.0/len(modification_logits) )
	
	modification_df_ks, dummy = analysis.ks_2samp_greater(modification_logits, null_logits)

	replacement_df_ks, dummy = analysis.ks_2samp_greater(replacement_logits, null_logits)

	bootstrapped_modification_df_kss.append(modification_df_ks)
	bootstrapped_replacement_df_kss.append(replacement_df_ks)

	abs_modification_df_ks, dummy = analysis.ks_2samp_greater(abs_modification_logits, abs_null_logits)

	abs_replacement_df_ks, dummy = analysis.ks_2samp_greater(abs_replacement_logits, abs_null_logits)

	bootstrapped_modification_absdf_kss.append(abs_modification_df_ks)
	bootstrapped_replacement_absdf_kss.append(abs_replacement_df_ks)


bootstrapped_modification_df_kss = numpy.array(bootstrapped_modification_df_kss)
bootstrapped_replacement_df_kss = numpy.array(bootstrapped_replacement_df_kss)

output_file.write("Modification other fraction positive P-value = %g\n" % (((bootstrapped_fraction_positives>=observed_fraction_positive).sum()+1.0)/(num_bootstraps+1.0)))


output_file.write("Modification delta_f KS P-value = %g\n" % (((bootstrapped_modification_df_kss>=observed_modification_df_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement delta_f KS P-value = %g\n" % (((bootstrapped_replacement_df_kss>=observed_replacement_df_ks).sum()+1.0)/(num_bootstraps+1.0)))

bootstrapped_modification_absdf_kss = numpy.array(bootstrapped_modification_absdf_kss)
bootstrapped_replacement_absdf_kss = numpy.array(bootstrapped_replacement_absdf_kss)

output_file.write("Modification abs(delta_f) KS P-value = %g\n" % (((bootstrapped_modification_absdf_kss>=observed_modification_absdf_ks).sum()+1.0)/(num_bootstraps+1.0)))

output_file.write("Replacement abs(delta_f) KS P-value = %g" % (((bootstrapped_replacement_absdf_kss>=observed_replacement_absdf_ks).sum()+1.0)/(num_bootstraps+1.0)))
