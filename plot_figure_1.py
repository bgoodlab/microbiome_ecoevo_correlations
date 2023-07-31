import sys
import numpy
import pylab
import parse_data
import analysis
from numpy.random import choice,shuffle,normal,uniform
from scipy.stats import ks_2samp
from scipy.stats import gaussian_kde
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import config
import stats_utils

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'

output_file = open("figures/figure_1_output.txt","w")

table_file = open("supplemental_data/TableS3.txt","w")
# Write header
table_file.write(", ".join(["Regression","Sample size","Num positive", "beta_phylum", "pvalue","beta_initial_diversity","pvalue"]))
table_file.write("\n")


#####
#
# Set up figure
#
#########

main_fig = plt.figure(figsize=(7,3.3))

outer_grid = gridspec.GridSpec(1,2, width_ratios=[5,2], wspace=0.3) 

inner_grid = gridspec.GridSpecFromSubplotSpec(3,2,width_ratios=[4,1],height_ratios=[0.05,1,1.2],wspace=0.1,subplot_spec=outer_grid[0])

right_grid = gridspec.GridSpecFromSubplotSpec(4,1,height_ratios=[0.15,0.9,0.33,0.82],hspace=0.1,subplot_spec=outer_grid[1])

family_axis = plt.Subplot(main_fig, inner_grid[1,0])
main_fig.add_subplot(family_axis)
family_axis.set_ylabel('Probability of event')
family_axis.semilogy([1],[2],'.')
family_axis.set_ylim([7e-03,1])

phylum_axis = plt.Subplot(main_fig, inner_grid[1,1])
main_fig.add_subplot(phylum_axis)
phylum_axis.semilogy([1],[2],'.')
phylum_axis.set_ylim([7e-03,1])
phylum_axis.set_yticklabels([])


entropy_axis = plt.Subplot(main_fig, right_grid[1])
main_fig.add_subplot(entropy_axis)
entropy_axis.set_ylabel('Pr[modification]')
entropy_axis.set_xlabel('Community diversity (Shannon)')
entropy_axis.set_xlim([2,5])

event_axis = plt.Subplot(main_fig, right_grid[3])
main_fig.add_subplot(event_axis)
event_axis.set_ylabel('Number of hosts')
event_axis.set_xlabel('Number of events')
event_axis.set_xlim([-4,4])

richness_fig = plt.figure(figsize=(2.5,2))
richness_grid = gridspec.GridSpec(1,1) 
richness_axis = plt.Subplot(richness_fig, richness_grid[0])
richness_fig.add_subplot(richness_axis)
richness_axis.set_ylabel('Pr[modification]')
richness_axis.set_xlabel('Community diversity (richness)')

bf_fig = plt.figure(figsize=(2.5,2))
bf_grid = gridspec.GridSpec(1,1) 
bf_axis = plt.Subplot(bf_fig, richness_grid[0])
bf_fig.add_subplot(bf_axis)
bf_axis.set_ylabel('Pr[modification]')
bf_axis.set_xlabel('log10(Firmicutes:Bacteroidetes)')

replacement_fig = plt.figure(figsize=(5,2))
replacement_grid = gridspec.GridSpec(1,2,width_ratios=[1,1],wspace=0.2) 

bac_replacement_axis = plt.Subplot(replacement_fig, replacement_grid[0])
replacement_fig.add_subplot(bac_replacement_axis)
bac_replacement_axis.set_xlabel('Shannon diversity, H')
bac_replacement_axis.set_ylabel('Fraction >= H')

firm_replacement_axis = plt.Subplot(replacement_fig, replacement_grid[1])
replacement_fig.add_subplot(firm_replacement_axis)
firm_replacement_axis.set_xlabel('Shannon diversity, H')


maxX=1
cNorm  = colors.Normalize(vmin=-maxX, vmax=maxX)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pylab.get_cmap('coolwarm') )

null_color = scalarMap.to_rgba(0)
modification_color = scalarMap.to_rgba(-1)
replacement_color = scalarMap.to_rgba(1)

import bacterial_phylogeny_utils

genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()
taxonomic_level="phylum"

sys.stderr.write("Loading data...\t")
sample_metadata_map = parse_data.parse_sample_metadata_map() 

# First look at all pairs (default is consecutive)
within_host_changes = parse_data.parse_within_host_changes("within_host_changes_all_pairs.txt")
within_host_events = analysis.calculate_within_host_events(within_host_changes)
host_classification_map = analysis.calculate_host_classifications(within_host_events)
species_event_map = analysis.collate_species_events(within_host_events)

num_events = numpy.array([0,0,0])
for species in species_event_map:
	num_events += species_event_map[species]

unique_hosts = set([])
for cohort,subject,sample_1,sample_2 in host_classification_map:
	unique_hosts.add((cohort,subject))
		
output_file.write("Full dataset = %d replacements, %d modifications, %d timepoint pairs, %s species, %d hosts\n" % (num_events[2],num_events[1],num_events.sum(),len(species_event_map), len(unique_hosts)))	

within_host_changes = parse_data.parse_within_host_changes()
abundance_matrix,speciess,samples = parse_data.parse_abundances()

# Now do for consecutive pairs (dereplicated)
within_host_events = analysis.calculate_within_host_events(within_host_changes)
host_classification_map = analysis.calculate_host_classifications(within_host_events)
host_records = host_classification_map.keys()
host_tabulations_map = analysis.calculate_host_tabulations(within_host_events)
host_event_map = analysis.collate_within_host_events(within_host_events)
host_ecological_map = analysis.calculate_host_ecological_changes(host_records,abundance_matrix,speciess,samples)
species_event_map = analysis.collate_species_events(within_host_events)
sys.stderr.write("Done!\n")	 

num_events = numpy.array([0,0,0])
for species in species_event_map:
	num_events += species_event_map[species]

unique_hosts = set([])
for cohort,subject,sample_1,sample_2 in host_classification_map:
	unique_hosts.add((cohort,subject))
		
output_file.write("Dereplicated dataset = %d replacements, %d modifications, %d timepoint pairs, %s species, %d hosts\n" % (num_events[2],num_events[1],num_events.sum(),len(species_event_map), len(unique_hosts)))	


# Bacteroides firmicutes ratio
bf0s = []
# Used to calculate B:F ratio
bacteroidetes_f0_idxs = []
firmicutes_f0_idxs = []
for species in speciess:
	
	phylum = bacterial_phylogeny_utils.get_taxon_name(species,taxonomic_level='phylum',genus_family_map=genus_family_map,genus_phylum_map=genus_phylum_map)
	
	bacteroidetes_f0_idxs.append(phylum=="Bacteroidetes")
	firmicutes_f0_idxs.append(phylum=="Firmicutes")
	
bacteroidetes_f0_idxs = numpy.array(bacteroidetes_f0_idxs)
firmicutes_f0_idxs = numpy.array(firmicutes_f0_idxs)


alpha0s = []
richness0s = []
modification_idxs = []
replacement_idxs = []

event_vector = []

phyla = []
firmicutes_idxs = []
bacteroidetes_idxs = []


host_alpha0s = []
allowed_phyla = ["Bacteroidetes","Firmicutes"]

species_names = []

for host_record in host_records:
	
	f0s = abundance_matrix[:,samples.index(host_record[2])]
	
	B_fraction = f0s[bacteroidetes_f0_idxs].sum()
	F_fraction = f0s[firmicutes_f0_idxs].sum()
	
	bf0 = numpy.log10(F_fraction/B_fraction)
	
	
	alpha0 = host_ecological_map[host_record][0]
	host_alpha0s.append(alpha0)
	
	richness0 = host_ecological_map[host_record][4]
	
	for species,event in host_event_map[host_record]:
		
		phylum = bacterial_phylogeny_utils.get_taxon_name(species,taxonomic_level=taxonomic_level,genus_family_map=genus_family_map,genus_phylum_map=genus_phylum_map)
		
		#if phylum not in allowed_phyla:
		#	 pass
			#continue
			 
		if phylum=='Firmicutes':
			firmicutes_idxs.append(True)
		else:
			firmicutes_idxs.append(False)
			
		if phylum=='Bacteroidetes':
			bacteroidetes_idxs.append(True)
		else:	
			bacteroidetes_idxs.append(False)
		
		phyla.append(phylum)   
		alpha0s.append(alpha0)
		richness0s.append(richness0)
		species_names.append(species)
		bf0s.append(bf0)
		
		event_vector.append(event)
		
		if event==0:
			modification_idxs.append(False)
			replacement_idxs.append(False)
		elif event==1: # modification
			modification_idxs.append(True)
			replacement_idxs.append(False)
		elif event==2: # replacement   
			modification_idxs.append(False)
			replacement_idxs.append(True)

alpha0s = numpy.array(alpha0s)
richness0s = numpy.array(richness0s)
bf0s = numpy.array(bf0s)

modification_idxs = numpy.array(modification_idxs)
replacement_idxs = numpy.array(replacement_idxs)
null_idxs = numpy.logical_not(numpy.logical_or(modification_idxs,replacement_idxs))
event_vector = numpy.array(event_vector)

firmicutes_idxs = numpy.array(firmicutes_idxs)
bacteroidetes_idxs = numpy.array(bacteroidetes_idxs)
host_alpha0s = numpy.array(host_alpha0s)
phyla = numpy.array(phyla)

sys.stderr.write("Generating bootstrapped event vector...\n")
bootstrapped_event_vectors = analysis.calculate_bootstrapped_event_vector(event_vector,species_names)


bootstrapped_null_idxss = (bootstrapped_event_vectors==0)
bootstrapped_modification_idxss = (bootstrapped_event_vectors==1)
bootstrapped_replacement_idxss = (bootstrapped_event_vectors==2)

sys.stderr.write("Done!\n")

def logistic_regression(ys, x1s, x2s, bootstrapped_ys=[]):
	
	
	logit_mod = sm.Logit(ys, sm.add_constant(numpy.c_[x1s,x2s]))
	logit_res = logit_mod.fit()

	intercept = logit_res.params[0]
	beta_alpha0 = logit_res.params[1]
	pvalue_alpha0 = logit_res.pvalues[1]
	beta_firmicutes = logit_res.params[2]
	pvalue_firmicutes = logit_res.pvalues[2]

	if len(bootstrapped_ys)>0:
		sys.stderr.write("Calculating permutation pvalue for logistic regression...\n")
		observed_pvalue = pvalue_alpha0
		num_bootstraps = bootstrapped_ys.shape[0]
		bootstrapped_pvalues = []
		for bootstrap_idx in range(0,num_bootstraps):
			if bootstrap_idx%1000==0:
				sys.stderr.write("%d\n" % bootstrap_idx)
			logit_mod = sm.Logit(bootstrapped_ys[bootstrap_idx,:], sm.add_constant(numpy.c_[x1s,x2s]))
			logit_res = logit_mod.fit(disp=0)
			bootstrapped_pvalues.append( logit_res.pvalues[1] )
		
		bootstrapped_pvalues = numpy.array(bootstrapped_pvalues)
	
		permuted_pvalue_alpha0 = ((bootstrapped_pvalues<=observed_pvalue).sum()+1.0)/(num_bootstraps+1.0)
		
		sys.stderr.write("Done!\n")
	
	else:
		permuted_pvalue_alpha0 = 1
	
	return intercept, beta_firmicutes, pvalue_firmicutes, beta_alpha0, pvalue_alpha0, permuted_pvalue_alpha0

sys.stderr.write("Analyzing rates per taxon...\n")
num_bootstraps = config.default_num_bootstraps # 100,000 in longest one
for type in ['species','genus','family','phylum']:
	
	# make coarse-grained version of data:
	
	label_id_map = {}
	id_label_map = []
	id_phylum_map = []
	ids = []
	events = []
	for cohort,subject,t0,t1,species,event in within_host_events:
		
		family = bacterial_phylogeny_utils.get_family_name(species, genus_family_map=genus_family_map)
		phylum = bacterial_phylogeny_utils.get_phylum_name(species, genus_phylum_map=genus_phylum_map)
		
		# get label
		if type=='family':
			label=family
		elif type=='phylum':
			label=phylum
		elif type=='genus':
			label = bacterial_phylogeny_utils.get_genus_name(species)	 
		else:
			label=species
		
		#if not phylum in good_phyla:
			#print species, family, phylum
			#continue
			
		#if label=='Veillonellaceae':
		#	 print cohort,subject,t0,t1,species,event
			
		if label not in label_id_map:
			id = len(id_label_map)
			label_id_map[label] = id
			id_label_map.append(label)
			id_phylum_map.append(phylum)
		
		id = label_id_map[label]
		
		ids.append(id)
		events.append(event)
	
	observed_event_matrix = analysis.calculate_matrix_from_lists(ids,events)
	observed_L = analysis.calculate_loglikelihood(observed_event_matrix)
	
	#print observed_L
	
	bootstrapped_Ls = []
	for bootstrap_idx in range(0,num_bootstraps):
		
		shuffle(events) # shuffle events across categories
		bootstrapped_event_matrix = analysis.calculate_matrix_from_lists(ids,events)
		bootstrapped_L = analysis.calculate_loglikelihood(bootstrapped_event_matrix)
		bootstrapped_Ls.append(bootstrapped_L)
		
	bootstrapped_Ls = numpy.array(bootstrapped_Ls)
	
	pvalue = ((bootstrapped_Ls>=observed_L).sum()+1.0)/(len(bootstrapped_Ls)+1.0)
	
	output_file.write("Loglikelihood test for %s rates: %d categories, p-value = %g\n" % (type, len(id_label_map), pvalue))
	
	# Now plot (only plot at family and phylum level)
	if type in ['family','phylum']:
		
		ntots = observed_event_matrix.sum(axis=1)
		n0s = observed_event_matrix[:,0]*1.0
		idxs = numpy.arange(0,len(ntots))
		
		sorted_phyla, negative_sorted_ntots, negative_sorted_n0s, sorted_idxs = zip(*sorted(zip(id_phylum_map, -ntots, -n0s, idxs),reverse=False))
		xticks = []
		xticklabels = []
		
		sorted_ntots = -1*numpy.array(negative_sorted_ntots)
		sorted_n0s = -1*numpy.array(negative_sorted_n0s)
		sorted_idxs = numpy.array(sorted_idxs)
		pooled_ntot = ntots.sum()
		pooled_n0 = n0s.sum()
		pooled_p = (pooled_ntot-pooled_n0)*1.0/pooled_ntot
		
		if type=='phylum':
			axis = phylum_axis
		else:
			axis = family_axis
			
		
		phylum_x_map = {}
		for x in range(0,len(sorted_idxs)):
			phylum = sorted_phyla[x]
			
			if phylum not in phylum_x_map:
				phylum_x_map[phylum] = []
			
			phylum_x_map[phylum].append(x)
			
		grey = False
		for phylum in sorted(phylum_x_map):
			min_x = min(phylum_x_map[phylum])
			max_x = max(phylum_x_map[phylum])
			
			if grey:
				axis.fill_between([min_x-0.5,max_x+0.5],[1e-07,1e-07],[1,1],color='0.8')
			
			grey = not grey
			
			
		# Now do another one with separate rep & mod
		for x in range(0,len(sorted_idxs)):
			
			
			ntot = sorted_ntots[x]
			id = sorted_idxs[x]
			
			nmod = observed_event_matrix[id,1]*1.0
			pmod = nmod*1.0/ntot
			plower,pupper = analysis.calculate_poisson_confidence_interval(nmod,ntot,alpha=0.05)
			
			plower = max([1e-07,plower])
			pmod = max([1e-07,pmod])
			
			#pylab.bar([x-0.2],[pmod],width=0.4,color='b')
			axis.plot([x-0.15,x-0.15], [plower,pupper],'-',color=modification_color,linewidth=0.5)
			axis.plot([x-0.15],[pmod],'.',color=modification_color)
			
			nrep = observed_event_matrix[id,2]*1.0
			prep = nrep*1.0/ntot
			plower,pupper = analysis.calculate_poisson_confidence_interval(nrep,ntot,alpha=0.05)
			
			plower = max([1e-07,plower])
			prep = max([1e-07,prep])
			
			#pylab.bar([x+0.2],[prep],width=0.4,color='r')
			axis.plot([x+0.15,x+0.15], [plower,pupper],'-',color=replacement_color,linewidth=0.5)
			axis.semilogy([x+0.15],[prep],'.',color=replacement_color)
			
			#print sorted_phyla
			
			xticks.append(x)
			xticklabels.append('%s (n=%d)' % (id_label_map[id],ntot))
			
		axis.set_xlim([xticks[0]-0.5,xticks[-1]+0.5])
		axis.set_xticks(xticks)
		axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=6)

phylum_axis.set_yticklabels([])

family_axis.plot([0],[2],'.',color=replacement_color,label='Replacement')
family_axis.plot([0],[2],'.',color=modification_color,label='Modification')

#family_axis.legend(loc=(0,1),frameon=False,ncol=2,numpoints=1,handletextpad=0.1)
family_axis.legend(loc=(0,0.72),frameon=False,numpoints=1,handletextpad=0.1)

sys.stderr.write("Done!\n")
		
#############################
#
# Now look at evolutionary rate vs community entropy w/ logistic regression
#
##############################

sys.stderr.write("Analyzing correlation w/ entropy...\n")

# Look within Bacteroidetes and firmicutes
good_idxs = numpy.logical_or(bacteroidetes_idxs,firmicutes_idxs)

intercept, beta_firmicutes, pvalue_firmicutes, beta_alpha0, pvalue_alpha0, permuted_pvalue_alpha0 = logistic_regression(replacement_idxs[good_idxs],alpha0s[good_idxs],firmicutes_idxs[good_idxs],bootstrapped_replacement_idxss[:,good_idxs])

output_file.write("Logistic regression for replacements:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_alpha0 = %g, pvalue = %g, permutation_pvalue = %g\n" % (beta_alpha0, pvalue_alpha0, permuted_pvalue_alpha0))

table_file.write(", ".join(["logit p(replacement) ~ phylum + shannon","%d" % good_idxs.sum(),"%d" % replacement_idxs[good_idxs].sum(), "%0.3f" % beta_firmicutes, "%0.3f" % pvalue_firmicutes, "%0.3f" % beta_alpha0, "%0.3f" % pvalue_alpha0]))
table_file.write("\n")

# Now do same thing for modifications
intercept, beta_firmicutes, pvalue_firmicutes, beta_alpha0, pvalue_alpha0, permuted_pvalue_alpha0 = logistic_regression(modification_idxs[good_idxs],alpha0s[good_idxs],firmicutes_idxs[good_idxs],bootstrapped_modification_idxss[:,good_idxs])

output_file.write("Logistic regression for modifications:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_alpha0 = %g, pvalue = %g, permutation_pvalue = %g\n" % (beta_alpha0, pvalue_alpha0, permuted_pvalue_alpha0))

table_file.write(", ".join(["logit p(modification) ~ phylum + shannon","%d" % good_idxs.sum(),"%d" % modification_idxs[good_idxs].sum(), "%0.3f" % beta_firmicutes, "%0.3f" % pvalue_firmicutes, "%0.3f" % beta_alpha0, "%0.3f" % pvalue_alpha0]))
table_file.write("\n")

# Now do same thing, but exclude 10th percentile 

alpha0_star = numpy.quantile(alpha0s[good_idxs],0.9)
good_idxs = good_idxs*(alpha0s<=alpha0_star)
intercept2, beta_firmicutes2, pvalue_firmicutes2, beta_alpha02, pvalue_alpha02, permuted_pvalue_alpha02 = logistic_regression(modification_idxs[good_idxs],alpha0s[good_idxs],firmicutes_idxs[good_idxs])

output_file.write("Logistic regression for modifications, excluding top 0.1 of entropy values (H=%g):\n" % alpha0_star)
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes2, pvalue_firmicutes2))
output_file.write("beta_alpha0 = %g, pvalue = %g\n" % (beta_alpha02, pvalue_alpha02))

table_file.write(", ".join(["logit p(modification) ~ phylum + shannon (bottom 90%)","%d" % good_idxs.sum(),"%d" % modification_idxs[good_idxs].sum(), "%0.3f" % beta_firmicutes, "%0.3f" % pvalue_firmicutes, "%0.3f" % beta_alpha0, "%0.3f" % pvalue_alpha0]))
table_file.write("\n")

# Now plot results for visualization

# Supplemental Fig
all_alphas = alpha0s[bacteroidetes_idxs]
replacement_alphas = alpha0s[bacteroidetes_idxs*replacement_idxs]

all_xs, all_survivals = stats_utils.calculate_unnormalized_survival_from_vector(all_alphas, min_x=all_alphas.min(), max_x=all_alphas.max())

replacement_xs, replacement_survivals = stats_utils.calculate_unnormalized_survival_from_vector(replacement_alphas, min_x=all_alphas.min(), max_x=all_alphas.max())
 
bac_replacement_axis.step(all_xs,all_survivals/all_survivals[0],color='0.7',where='pre')
bac_replacement_axis.step(replacement_xs,replacement_survivals/replacement_survivals[0],color='#08519c',where='pre')

all_alphas = alpha0s[firmicutes_idxs]
replacement_alphas = alpha0s[firmicutes_idxs*replacement_idxs]

all_xs, all_survivals = stats_utils.calculate_unnormalized_survival_from_vector(all_alphas, min_x=all_alphas.min(), max_x=all_alphas.max())

replacement_xs, replacement_survivals = stats_utils.calculate_unnormalized_survival_from_vector(replacement_alphas, min_x=all_alphas.min(), max_x=all_alphas.max())
 
firm_replacement_axis.step(all_xs,all_survivals/all_survivals[0],color='0.7',where='pre')
firm_replacement_axis.step(replacement_xs,replacement_survivals/replacement_survivals[0],color='#6baed6',where='pre')

replacement_fig.savefig('figures/supplemental_replacement.pdf',bbox_inches='tight')

# Main Fig
for phylum_idxs,Q,color,label,firmicutes_indicator in zip([bacteroidetes_idxs,firmicutes_idxs],[10,3],['#08519c','#6baed6'],['Bacteroidetes','Firmicutes'],[0,1]):

	# For legend
	entropy_axis.plot([0],[1],'.',color=color,label=label)
	
	quantiles = numpy.array([numpy.quantile(alpha0s[phylum_idxs],(1.0/Q)*i) for i in range(0,Q+1)])
	
	print(quantiles[-2])
	numerators = ((alpha0s[:,None]>=quantiles[None,0:-1]) * (alpha0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]*modification_idxs[:,None]).sum(axis=0)
	
	denominators = ((alpha0s[:,None]>=quantiles[None,0:-1]) * (alpha0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]).sum(axis=0)
	
	rates = numerators*1.0/denominators
	sigmas = numpy.sqrt(rates*(1-rates)/denominators)
		
	xs = (quantiles[0:-1]+quantiles[1:])/2.0
	
	line, = entropy_axis.plot(xs,rates,'.',color=color)
	color = pylab.getp(line,'color')
	for q_idx in range(0,len(rates)):
		
		plower,pupper = analysis.calculate_poisson_confidence_interval(numerators[q_idx],denominators[q_idx],alpha=0.05)
			
		entropy_axis.plot([xs[q_idx],xs[q_idx]],[plower, pupper],'-',color=color)
	
	logit_mod = sm.Logit(modification_idxs[phylum_idxs],	 sm.add_constant(alpha0s[phylum_idxs]))
	logit_res = logit_mod.fit()
	pvalue_alpha0 = logit_res.pvalues[1]
	print(" ")
	theory_alphas = numpy.linspace(alpha0s.min(),alpha0s.max(),50)

	theory_logits = intercept+beta_alpha0*theory_alphas+beta_firmicutes*firmicutes_indicator
	
	theory_rates = numpy.exp(theory_logits)/(1+numpy.exp(theory_logits))
	
	entropy_axis.plot(theory_alphas,theory_rates,'-',color=color)
	
	numerator = ((alpha0s<=quantiles[-2])*phylum_idxs*modification_idxs).sum(axis=0)
	
	denominator = ((alpha0s<=quantiles[-2])*phylum_idxs).sum(axis=0)
	
	rate = numerator*1.0/denominator
	entropy_axis.plot(theory_alphas, rate*numpy.ones_like(theory_alphas),':',color=color)
	
entropy_axis.set_ylim([0,0.35])
entropy_axis.set_xlim([1.5,5])
entropy_axis.xaxis.set_tick_params(pad=2)
entropy_axis.yaxis.set_tick_params(pad=2)

entropy_axis.legend(loc=(-0.05,1),frameon=False,ncol=2,numpoints=1,handletextpad=0.1,columnspacing=0.1)

#############################
#
# Supplemental version w/ richness instead of entropy 
#
##############################

# Look within Bacteroidetes and firmicutes
good_idxs = numpy.logical_or(bacteroidetes_idxs,firmicutes_idxs)

# Make logistic regression using entropy and firmicutes idxs as parameters
logit_mod = sm.Logit(replacement_idxs[good_idxs], sm.add_constant(numpy.c_[richness0s[good_idxs],firmicutes_idxs[good_idxs]]))
logit_res = logit_mod.fit()

intercept = logit_res.params[0]
beta_richness0 = logit_res.params[1]
pvalue_richness0 = logit_res.pvalues[1]
beta_firmicutes = logit_res.params[2]
pvalue_firmicutes = logit_res.pvalues[2]


output_file.write("Logistic regression for replacements:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_richness0 = %g, pvalue = %g\n" % (beta_richness0, pvalue_richness0))

table_file.write(", ".join(["logit p(replacement) ~ phylum + richness","%d" % good_idxs.sum(),"%d" % replacement_idxs[good_idxs].sum(), "%0.3f" % beta_firmicutes, "%0.3f" % pvalue_firmicutes, "%0.3f" % beta_richness0, "%0.3f" % pvalue_richness0]))
table_file.write("\n")

# Now do same thing for modifications
logit_mod = sm.Logit(modification_idxs[good_idxs], sm.add_constant(numpy.c_[richness0s[good_idxs],firmicutes_idxs[good_idxs]]))
logit_res = logit_mod.fit()

intercept = logit_res.params[0]
beta_richness0 = logit_res.params[1]
pvalue_richness0 = logit_res.pvalues[1]
beta_firmicutes = logit_res.params[2]
pvalue_firmicutes = logit_res.pvalues[2]

output_file.write("Logistic regression for modifications:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_richness0 = %g, pvalue = %g\n" % (beta_richness0, pvalue_richness0))

table_file.write(", ".join(["logit p(modification) ~ phylum + richness","%d" % good_idxs.sum(),"%d" % modification_idxs[good_idxs].sum(), "%0.3f" % beta_firmicutes, "%0.3f" % pvalue_firmicutes, "%0.3f" % beta_richness0, "%0.3f" % pvalue_richness0]))
table_file.write("\n")

# Now plot results for visualization
for phylum_idxs,Q,color,label,firmicutes_indicator in zip([bacteroidetes_idxs,firmicutes_idxs],[10,3],['#08519c','#6baed6'],['Bacteroidetes','Firmicutes'],[0,1]):

	# For legend
	richness_axis.plot([0],[1],'.',color=color,label=label)
	
	quantiles = numpy.array([numpy.quantile(richness0s[phylum_idxs],(1.0/Q)*i) for i in range(0,Q+1)])
	
	print(quantiles[-2])
	numerators = ((richness0s[:,None]>=quantiles[None,0:-1]) * (richness0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]*modification_idxs[:,None]).sum(axis=0)
	
	denominators = ((richness0s[:,None]>=quantiles[None,0:-1]) * (richness0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]).sum(axis=0)
	
	rates = numerators*1.0/denominators
	sigmas = numpy.sqrt(rates*(1-rates)/denominators)
		
	xs = (quantiles[0:-1]+quantiles[1:])/2.0
	
	line, = richness_axis.plot(xs,rates,'.',color=color)
	color = pylab.getp(line,'color')
	for q_idx in range(0,len(rates)):
		
		plower,pupper = analysis.calculate_poisson_confidence_interval(numerators[q_idx],denominators[q_idx],alpha=0.05)
			
		richness_axis.plot([xs[q_idx],xs[q_idx]],[plower, pupper],'-',color=color)
	
	logit_mod = sm.Logit(modification_idxs[phylum_idxs],	 sm.add_constant(richness0s[phylum_idxs]))
	logit_res = logit_mod.fit()
	pvalue_alpha0 = logit_res.pvalues[1]
	
	print(label, pvalue_alpha0)

	theory_richnesses = numpy.linspace(richness0s.min(),richness0s.max(),50)

	theory_logits = intercept+beta_richness0*theory_richnesses+beta_firmicutes*firmicutes_indicator
	
	theory_rates = numpy.exp(theory_logits)/(1+numpy.exp(theory_logits))
	
	richness_axis.plot(theory_richnesses,theory_rates,'-',color=color)
	
	numerator = ((richness0s<=quantiles[-2])*phylum_idxs*modification_idxs).sum(axis=0)
	
	denominator = ((richness0s<=quantiles[-2])*phylum_idxs).sum(axis=0)
	
	rate = numerator*1.0/denominator
	richness_axis.plot(theory_richnesses, rate*numpy.ones_like(theory_richnesses),':',color=color)
	
richness_axis.set_ylim([0,0.35])
#richness_axis.set_xlim([1.5,5])
#richness_axis.xaxis.set_tick_params(pad=2)
#entropy_axis.yaxis.set_tick_params(pad=2)

#richness_axis.legend(loc=(-0.05,1),frameon=False,ncol=2,numpoints=1,handletextpad=0.1,columnspacing=0.1)

#############################
#
# Supplemental version w/ bacteroides:firmicutes fraction instead of entropy 
#
##############################

# Look within Bacteroidetes and firmicutes
good_idxs = numpy.logical_or(bacteroidetes_idxs,firmicutes_idxs)

# Make logistic regression using entropy and firmicutes idxs as parameters
logit_mod = sm.Logit(replacement_idxs[good_idxs], sm.add_constant(numpy.c_[bf0s[good_idxs],firmicutes_idxs[good_idxs]]))
logit_res = logit_mod.fit()

intercept = logit_res.params[0]
beta_bf0 = logit_res.params[1]
pvalue_bf0 = logit_res.pvalues[1]
beta_firmicutes = logit_res.params[2]
pvalue_firmicutes = logit_res.pvalues[2]


output_file.write("Logistic regression for replacements:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_bf0 = %g, pvalue = %g\n" % (beta_bf0, pvalue_bf0))

# Now do same thing for modifications
logit_mod = sm.Logit(modification_idxs[good_idxs], sm.add_constant(numpy.c_[bf0s[good_idxs],firmicutes_idxs[good_idxs]]))
logit_res = logit_mod.fit()

intercept = logit_res.params[0]
beta_bf0 = logit_res.params[1]
pvalue_bf0 = logit_res.pvalues[1]
beta_firmicutes = logit_res.params[2]
pvalue_firmicutes = logit_res.pvalues[2]

output_file.write("Logistic regression for modifications:\n")
output_file.write("beta_firmicutes = %g, pvalue = %g\n" % (beta_firmicutes, pvalue_firmicutes))
output_file.write("beta_bf0 = %g, pvalue = %g\n" % (beta_bf0, pvalue_bf0))

# Now plot results for visualization
for phylum_idxs,Q,color,label,firmicutes_indicator in zip([bacteroidetes_idxs,firmicutes_idxs],[10,3],['#08519c','#6baed6'],['Bacteroidetes','Firmicutes'],[0,1]):

	# For legend
	bf_axis.plot([0],[1],'.',color=color,label=label)
	
	quantiles = numpy.array([numpy.quantile(bf0s[phylum_idxs],(1.0/Q)*i) for i in range(0,Q+1)])
	
	print(quantiles[-2])
	numerators = ((bf0s[:,None]>=quantiles[None,0:-1]) * (bf0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]*modification_idxs[:,None]).sum(axis=0)
	
	denominators = ((bf0s[:,None]>=quantiles[None,0:-1]) * (bf0s[:,None]<=quantiles[None,1:])*phylum_idxs[:,None]).sum(axis=0)
	
	rates = numerators*1.0/denominators
	sigmas = numpy.sqrt(rates*(1-rates)/denominators)
		
	xs = (quantiles[0:-1]+quantiles[1:])/2.0
	
	line, = bf_axis.plot(xs,rates,'.',color=color)
	color = pylab.getp(line,'color')
	for q_idx in range(0,len(rates)):
		
		plower,pupper = analysis.calculate_poisson_confidence_interval(numerators[q_idx],denominators[q_idx],alpha=0.05)
			
		bf_axis.plot([xs[q_idx],xs[q_idx]],[plower, pupper],'-',color=color)
	
	logit_mod = sm.Logit(modification_idxs[phylum_idxs],	 sm.add_constant(richness0s[phylum_idxs]))
	logit_res = logit_mod.fit()
	pvalue_bf0 = logit_res.pvalues[1]
	
	print(label, pvalue_bf0)

	dq = (quantiles[-1]-quantiles[0])/10.0
	theory_bfs = numpy.linspace(quantiles[0]-dq,quantiles[-1]+dq,50)

	theory_logits = intercept+beta_bf0*theory_bfs+beta_firmicutes*firmicutes_indicator
	
	theory_rates = numpy.exp(theory_logits)/(1+numpy.exp(theory_logits))
	
	bf_axis.plot(theory_bfs,theory_rates,'-',color=color)
	
	numerator = ((bf0s<=quantiles[-2])*phylum_idxs*modification_idxs).sum(axis=0)
	
	denominator = ((bf0s<=quantiles[-2])*phylum_idxs).sum(axis=0)
	
	rate = numerator*1.0/denominator
	bf_axis.plot(theory_bfs, rate*numpy.ones_like(theory_bfs),':',color=color)
	
bf_axis.set_ylim([0,0.35])
#richness_axis.set_xlim([1.5,5])
#richness_axis.xaxis.set_tick_params(pad=2)
#entropy_axis.yaxis.set_tick_params(pad=2)

#richness_axis.legend(loc=(-0.05,1),frameon=False,ncol=2,numpoints=1,handletextpad=0.1,columnspacing=0.1)

########
#
# Plot number of events of each type
#
########
	
totals = []
modification_positives = []
replacement_positives = []
ks = numpy.arange(1,20)
for k in ks:
	
	num_hosts = 0
	num_modification_hosts = 0
	num_replacement_hosts = 0
	
	for host_record in host_tabulations_map:
		
		ntot,nmod,nrep = host_tabulations_map[host_record]
		
		if ntot>=k:
			num_hosts+=1
		
		if (nmod==k) and (nrep==0):
			num_modification_hosts += 1	  
			
		if ((nmod+nrep)==k) and (nrep>0): 
			num_replacement_hosts += 1
			
			if k>3:
					print(ntot,nmod,nrep)
					print(host_record)
					print(host_event_map[host_record])
			
	modification_positives.append(num_modification_hosts)
	totals.append(num_hosts)
	replacement_positives.append(num_replacement_hosts)
	
observed_modification_positives = numpy.array(modification_positives)
observed_totals = numpy.array(totals)
observed_replacement_positives = numpy.array(replacement_positives)
print(observed_totals)
print(observed_modification_positives)
print(observed_replacement_positives)
print("Boostrapping...")
num_bootstraps=1000
bootstrapped_within_host_eventss = analysis.calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=num_bootstraps,bootstrap_type='species')
print("Processing bootstraps!")
bootstrapped_modification_positives = []
bootstrapped_replacement_positives = []

for idx in range(0,num_bootstraps):
	#print idx
	#host_event_map = analysis.collate_within_host_events(bootstrapped_within_host_eventss[idx])
	host_tabulations_map = analysis.calculate_host_tabulations(bootstrapped_within_host_eventss[idx])

	totals = []
	modification_positives = []
	replacement_positives = []
	for k in ks:
	
		num_hosts = 0
		num_modification_hosts = 0
		num_replacement_hosts = 0
	
		for host_record in host_tabulations_map:
			
			ntot,nmod,nrep = host_tabulations_map[host_record]
		
			if ntot>=k:
				num_hosts+=1
		
			if (nmod==k) and (nrep==0):
				num_modification_hosts += 1	  
			
			if ((nmod+nrep)==k) and (nrep>0): 
				num_replacement_hosts += 1
				
				
			
		modification_positives.append(num_modification_hosts)
		totals.append(num_hosts)
		replacement_positives.append(num_replacement_hosts)
	
	bootstrapped_modification_positives.append(modification_positives)
	bootstrapped_replacement_positives.append(replacement_positives)
	
bootstrapped_modification_positives = numpy.array(bootstrapped_modification_positives)
bootstrapped_replacement_positives = numpy.array(bootstrapped_replacement_positives)

null_modification_positives = bootstrapped_modification_positives.mean(axis=0)
null_replacement_positives = bootstrapped_replacement_positives.mean(axis=0)

# Set up figure

event_axis.set_xlim([-5,5])
event_axis.set_xticks([-4,-3,-2,-1,0,1,2,3,4])
event_axis.set_xticklabels([4,3,2,1,0,1,2,3,4])

event_axis.bar(-ks,observed_modification_positives,width=0.8,color=modification_color)
event_axis.bar(ks,observed_replacement_positives,width=0.8,color=replacement_color)
event_axis.bar([0],totals[0]-observed_replacement_positives.sum()-observed_modification_positives.sum(),width=0.8,color=null_color)

null_xs = numpy.hstack([-ks[::-1],[0],ks])
null_ys = numpy.hstack([null_modification_positives[::-1],[totals[0]-null_replacement_positives.sum()-null_modification_positives.sum()],null_replacement_positives])

#print null_xs
#print null_ys
event_axis.semilogy(null_xs,null_ys,'k.:',linewidth=0.5,markeredgewidth=1,markeredgecolor='k',markerfacecolor='w')

event_axis.set_ylim([2e-01,200])

event_axis.xaxis.set_tick_params(pad=2)
event_axis.yaxis.set_tick_params(pad=2)

main_fig.savefig('figures/figure_1.pdf',bbox_inches='tight')   
richness_fig.savefig('figures/supplemental_fig_2.pdf',bbox_inches='tight')	 
bf_fig.savefig('figures/supplemental_bf_ratio.pdf',bbox_inches='tight')	 