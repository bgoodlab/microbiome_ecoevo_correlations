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

output_file = open("figure_3_output.txt","w")

import bacterial_phylogeny_utils
genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()
taxonomic_level = 'family'

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

num_bootstraps = 10000

#####
#
# Set up figure
#
#########


main_fig = plt.figure(figsize=(6,2))

outer_grid = gridspec.GridSpec(1,3, width_ratios=[1.2,1.45,1.05], wspace=0.6) 

left_grid = gridspec.GridSpecFromSubplotSpec(2,1,height_ratios=[1,1.8],subplot_spec=outer_grid[0],hspace=0.25)

bias_axis = plt.Subplot(main_fig, left_grid[0])
main_fig.add_subplot(bias_axis)
bias_axis.set_ylabel('% positive')

cdf_axis = plt.Subplot(main_fig, left_grid[1])
main_fig.add_subplot(cdf_axis)
cdf_axis.set_ylabel('Fraction species')
cdf_axis.set_xlabel('Min fold change')

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


sample_metadata_map = parse_data.parse_sample_metadata_map() 
within_host_changes = parse_data.parse_within_host_changes()
abundance_matrix,speciess,samples = parse_data.parse_abundances()
    
within_host_events = analysis.calculate_within_host_events(within_host_changes)

host_event_map = analysis.collate_within_host_events(within_host_events)

host_records = host_event_map.keys()

host_classification_map = analysis.calculate_host_classifications(within_host_events)
host_ecological_map = analysis.calculate_host_ecological_changes(host_classification_map,abundance_matrix,speciess,samples) 
js_item_map = analysis.calculate_js_item_map(host_classification_map,abundance_matrix,speciess,samples)    

species_event_map = analysis.collate_species_events(within_host_events)
    

null_fs = []
null_jss = []
modification_fs = []
modification_jss = []
replacement_fs = []
replacement_jss = []

species_species_map = {}
species_taxon_map = {}

######################################################
#
# Loop over all quasi-phaseable species/host combos
# and calculate fold change of focal species
#
######################################################
for record_idx in xrange(0,len(within_host_events)):
    
    cohort,subject,t0,t1,species,event = within_host_events[record_idx]
        
    host_record = (cohort,subject,t0,t1)

    if species not in species_species_map:
        
        species_idxs = bacterial_phylogeny_utils.get_idxs_from_same_taxon(species,speciess,'species',genus_family_map,genus_phylum_map)
        species_species_map[species] = species_idxs
        
        species_idxs = bacterial_phylogeny_utils.get_idxs_from_same_taxon(species,speciess,taxonomic_level,genus_family_map,genus_phylum_map)
        species_taxon_map[species] = species_idxs
        
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
taxon_ps = []
modification_idxs = []
replacement_idxs = []


bootstrapped_species_pss = []
bootstrapped_taxon_pss = []


for host_idx in xrange(0,len(host_records)):
    
    host_record = host_records[host_idx]
    
    js_items = js_item_map[host_record]
    
    
    
    if host_classification_map[host_record]>0:
        
        focal_species_idxs = numpy.zeros_like(js_items)
        taxon_idxs = numpy.zeros_like(js_items)
        
        host_species_weights = []
        host_species_taxon_idxs = []
        host_species_species_idxs = []
        host_num_events = 0
        
        for species,event in host_event_map[host_record]:
            
            if event==host_classification_map[host_record]:
                
                focal_species_idxs += species_species_map[species]
                taxon_idxs += species_taxon_map[species]
                host_num_events += 1
                
                
            host_species_weights.append( species_event_map[species][event]*1.0/sum(species_event_map[species]))
            host_species_species_idxs.append(species_species_map[species])
            host_species_taxon_idxs.append(species_taxon_map[species])
                
        species_p = js_items[focal_species_idxs>0].sum()/js_items.sum()
        taxon_p = js_items[taxon_idxs>0].sum()/js_items.sum()
        
        species_ps.append(species_p)
        taxon_ps.append(taxon_p)        
        
        host_species_idxs = numpy.arange(0,len(host_species_weights))
        host_species_weights = numpy.array(host_species_weights)
        host_species_weights = host_species_weights/host_species_weights.sum()
        host_species_taxon_idxs = numpy.array(host_species_taxon_idxs)
        host_species_species_idxs = numpy.array(host_species_species_idxs)
        
        if host_classification_map[host_record]==1:
            modification_idxs.append(True)
            replacement_idxs.append(False)
        else:
            modification_idxs.append(False)
            replacement_idxs.append(True) 
            
        # Now do bootstrapped version
        bootstrapped_taxon_ps = []
        bootstrapped_species_ps = []
        for bootstrap_idx in xrange(0,num_bootstraps):
            
            bootstrapped_event_idxs = choice(host_species_idxs,size=host_num_events,replace=False,p=host_species_weights)
            
                        
            bootstrapped_taxon_idxs = host_species_taxon_idxs[bootstrapped_event_idxs,:].sum(axis=0)
            
            bootstrapped_taxon_p = js_items[bootstrapped_taxon_idxs>0].sum()/js_items.sum()
            
            bootstrapped_taxon_ps.append(bootstrapped_taxon_p)
            
            bootstrapped_species_idxs = host_species_species_idxs[bootstrapped_event_idxs,:].sum(axis=0)
            
            bootstrapped_species_p = js_items[bootstrapped_species_idxs>0].sum()/js_items.sum()
            
            bootstrapped_species_ps.append(bootstrapped_species_p)
            
        bootstrapped_taxon_pss.append(bootstrapped_taxon_ps)       
        bootstrapped_species_pss.append(bootstrapped_species_ps)       

replacement_idxs = numpy.array(replacement_idxs)
modification_idxs = numpy.array(modification_idxs)
species_ps = numpy.array(species_ps)
taxon_ps = numpy.array(taxon_ps)
bootstrapped_taxon_pss = numpy.array(bootstrapped_taxon_pss)
bootstrapped_species_pss = numpy.array(bootstrapped_species_pss)


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
cdf_axis.set_xticklabels([],minor=True)
bias_axis.set_xlim([1,10])
bias_axis.set_ylim([25,75])
bias_axis.set_xticklabels([],minor=True)
bias_axis.set_xticklabels([],minor=False)


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
for bootstrap_idx in xrange(0,num_bootstraps):
    mean = bootstrapped_species_pss[modification_idxs,bootstrap_idx].mean() 
    bootstrapped_means.append(mean)   
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("modification focal species percent explained mean pvalue = %g\n" % pvalue)

bootstrapped_ps = bootstrapped_species_pss[replacement_idxs].flatten()
observed_mean = species_ps[replacement_idxs].mean()
bootstrapped_means = []
for bootstrap_idx in xrange(0,num_bootstraps):
    mean = bootstrapped_species_pss[replacement_idxs,bootstrap_idx].mean() 
    bootstrapped_means.append(mean)   
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("replacement focal species percent explained mean pvalue = %g\n" % pvalue)

bootstrapped_ps = bootstrapped_taxon_pss[modification_idxs].flatten()
observed_mean = taxon_ps[modification_idxs].mean()
bootstrapped_means = []
for bootstrap_idx in xrange(0,num_bootstraps):
    mean = bootstrapped_taxon_pss[modification_idxs,bootstrap_idx].mean() 
    bootstrapped_means.append(mean)   
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("modification focal %s percent explained mean pvalue = %g\n" % (taxonomic_level,pvalue))

bootstrapped_ps = bootstrapped_taxon_pss[replacement_idxs].flatten()
observed_mean = taxon_ps[replacement_idxs].mean()
bootstrapped_means = []
for bootstrap_idx in xrange(0,num_bootstraps):
    mean = bootstrapped_taxon_pss[replacement_idxs,bootstrap_idx].mean() 
    bootstrapped_means.append(mean)   
bootstrapped_means = numpy.array(bootstrapped_means)
pvalue = ((bootstrapped_means>=observed_mean).sum()+1.0)*1.0/(len(bootstrapped_means)+1.0)
output_file.write("replacement focal %s percent explained mean pvalue = %g\n" % (taxonomic_level,pvalue))

species_axis.plot(100*numpy.sort(species_ps[modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[modification_idxs]))),'-',color=modification_color)

species_axis.plot(100*numpy.sort(species_ps[replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[replacement_idxs]))),'-',color=replacement_color)

species_axis.set_xlim([0,50])

bootstrapped_ps = bootstrapped_taxon_pss.flatten()
genus_axis.plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'k:',linewidth=0.5)

genus_axis.plot(100*numpy.sort(taxon_ps[modification_idxs]), numpy.flip(numpy.linspace(0, 1, len(species_ps[modification_idxs]))),'-',color=modification_color)

genus_axis.plot(100*numpy.sort(taxon_ps[replacement_idxs]), numpy.flip(numpy.linspace(0, 1, len(taxon_ps[replacement_idxs]))),'-',color=replacement_color)

#bootstrapped_ps = bootstrapped_taxon_pss[replacement_idxs].flatten()
#genus_axis.plot(100*numpy.sort(bootstrapped_ps), numpy.flip(numpy.linspace(0, 1, len(bootstrapped_ps))),'-',color=replacement_color,alpha=0.5)


genus_axis.set_xlim([0,50])
genus_axis.set_xticks([0,25,50])

species_axis.set_xticklabels([])

species_axis.plot([200],[1],'k-',label='Observed')
species_axis.plot([200],[1],'k:',linewidth=0.5,label='Expected')
species_axis.legend(loc='upper right',frameon=False)

main_fig.savefig('figure_3.pdf',bbox_inches='tight')


# Bootstrapping!
print "Bootstrapping..."
bootstrapped_within_host_eventss = analysis.calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=num_bootstraps,bootstrap_type='genus')
print "Processing bootstraps!"
bootstrapped_modification_df_kss = []
bootstrapped_replacement_df_kss = []
bootstrapped_modification_absdf_kss = []
bootstrapped_replacement_absdf_kss = []

bootstrapped_fraction_positives = []

for idx in xrange(0,num_bootstraps):    

    if idx%100==0:
        print idx
    
    within_host_events = bootstrapped_within_host_eventss[idx]
    host_classification_map = analysis.calculate_host_classifications(within_host_events)

    null_fs = []
    modification_fs = []
    replacement_fs = []

    for record_idx in xrange(0,len(within_host_events)):
    
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
