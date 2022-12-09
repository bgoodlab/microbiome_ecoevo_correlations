import numpy
from numpy.random import shuffle
from scipy.spatial.distance import jensenshannon 
from scipy.stats import entropy,gamma
import diversity_utils
from math import log
import sys
import bacterial_phylogeny_utils
import config


replacement_threshold = config.default_replacement_threshold

NO_EVENT = 0
MODIFICATION_EVENT = 1
REPLACEMENT_EVENT = 2



def ks_2samp_greater(data,null):
    
    data = numpy.array(data)
    null = numpy.array(null)
    data_n = len(data)
    null_n = len(null)
    
    xs = numpy.array(sorted(set(data)|set(null)))
    
    data_survival = (data[:,None]>=xs[None,:]).sum(axis=0)*1.0/data_n
    null_survival = (null[:,None]>=xs[None,:]).sum(axis=0)*1.0/null_n
    
    xstar = xs[(data_survival-null_survival).argmax()]
    ks_distance = (data_survival-null_survival).max()
    
    return ks_distance, xstar
    
def ks_survivals(data,null):
    data = numpy.array(data)
    null = numpy.array(null)
    data_n = len(data)
    null_n = len(null)
    
    xs = numpy.array(sorted(set(data)|set(null)))
    
    data_survival = (data[:,None]>=xs[None,:]).sum(axis=0)*1.0/data_n
    null_survival = (null[:,None]>=xs[None,:]).sum(axis=0)*1.0/null_n
    
    return xs,data_survival,null_survival

def calculate_within_host_events(within_host_changes):
    
    within_host_events = []
    for record_idx in range(0,len(within_host_changes)):
            
        cohort,subject,t0,t1,species,Lsnp,ksnp,Lprivate,kprivate = within_host_changes[record_idx]
            
        event = 0
        if ksnp > 0: # or kgene>0:
            event = 1
            if ksnp > replacement_threshold:
                event = 2
                    
        within_host_events.append((cohort,subject,t0,t1,species,event))
    
    return within_host_events  

def get_genus(species_name):
    return species_name.split("_")[0]
    

def calculate_bootstrapped_event_vector(event_vector,species_names,num_bootstraps=config.default_num_bootstraps,bootstrap_type=config.default_taxonomic_level):
    
    genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
    genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

    species_event_vector_map = {}
    species_idx_map = {}
    for idx in range(0,len(event_vector)):
        event = event_vector[idx]
        species = species_names[idx]
        
        if bootstrap_type=='all':
            species_record = 'all'
        else:
            taxon = bacterial_phylogeny_utils.get_taxon_name(species,bootstrap_type,genus_family_map,genus_phylum_map)
            species_record = taxon
                
        if species_record not in species_event_vector_map:
            species_event_vector_map[species_record] = []
            species_idx_map[species_record] = []

        species_event_vector_map[species_record].append(event)
        species_idx_map[species_record].append(idx)
    
    bootstrapped_event_vectors = []
    # Now can permute within each taxon
    for bootstrap_idx in range(0,num_bootstraps):
        
        bootstrapped_event_vector = numpy.zeros_like(event_vector)
        
        for species_record in species_event_vector_map:
            shuffle(species_event_vector_map[species_record])
            
            for event, idx in zip(species_event_vector_map[species_record], species_idx_map[species_record]):
                
                bootstrapped_event_vector[idx] = event
            
        bootstrapped_event_vectors.append(bootstrapped_event_vector)
    
    bootstrapped_event_vectors = numpy.array(bootstrapped_event_vectors)
    return bootstrapped_event_vectors
    
def calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=config.default_num_bootstraps,bootstrap_type=config.default_taxonomic_level):

    genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
    genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

    # can be 'all', 'species','genus'
    
    species_event_map = {}
    for record_idx in range(0,len(within_host_events)):
        cohort,subject,t0,t1,species,event = within_host_events[record_idx]
        
        if bootstrap_type=='all':
            species_record = (cohort)
        else:
            taxon = bacterial_phylogeny_utils.get_taxon_name(species,bootstrap_type,genus_family_map,genus_phylum_map)
            species_record = (cohort, taxon)
                
        if species_record not in species_event_map:
            species_event_map[species_record] = [[],[]]

            
        species_event_map[species_record][0].append((cohort,subject,t0,t1,species))
        species_event_map[species_record][1].append(event)
        
    bootstrapped_within_host_eventss = []
        
    for bootstrapped_idx in range(0,num_bootstraps):
        
        bootstrapped_within_host_events = []
        # shuffle events across species
        for species_record in species_event_map:
            shuffle(species_event_map[species_record][1])
                
            for i in range(0,len(species_event_map[species_record][0])):
                cohort,subject,t0,t1,species = species_event_map[species_record][0][i]
                event = species_event_map[species_record][1][i]
                bootstrapped_within_host_events.append((cohort, subject,t0,t1,species,event))
            
        bootstrapped_within_host_eventss.append(bootstrapped_within_host_events)
    
    return bootstrapped_within_host_eventss

# Turns list of within host events in to map from host to events
def collate_within_host_events(within_host_events):
    
    host_event_map = {}
    
    for record_idx in range(0,len(within_host_events)):
        
        cohort,subject,t0,t1,species,event = within_host_events[record_idx]
        
        host_record = (cohort,subject,t0,t1)
        
        if host_record not in host_event_map:
            host_event_map[host_record] = []
        
        host_event_map[host_record].append((species,event))
        
    return host_event_map
    
# Turns list of within host events in to map from host to events
def collate_within_host_events_by_species(within_host_events):
    
    host_event_map = {}
    
    for record_idx in range(0,len(within_host_events)):
        
        cohort,subject,t0,t1,species,event = within_host_events[record_idx]
        
        host_record = (cohort,subject,t0,t1)
        
        if host_record not in host_event_map:
            host_event_map[host_record] = {}
        
        host_event_map[host_record][species] = event
        
    return host_event_map
    
def calculate_host_classifications(within_host_events):
    
    host_event_map = collate_within_host_events(within_host_events)
    host_classification_map = {}
    for host_record in host_event_map:
    
        host_event = 0
        for species,event in host_event_map[host_record]:
            if event>0 and host_event==0:
                host_event = event
            if event==2:
                host_event = event
        
        host_classification_map[host_record] = host_event
    
    return host_classification_map

def calculate_host_tabulations(within_host_events):
    
    host_event_map = collate_within_host_events(within_host_events)
    host_tabulation_map = {}
    for host_record in host_event_map:
    
        num_modifications = 0
        num_replacements = 0
        num_total = 0
        for species,event in host_event_map[host_record]:
            if event==1:
                num_modifications+=1
            elif event==2:
                num_replacements+=1
            num_total += 1
        
        host_tabulation_map[host_record] = (num_total,num_modifications,num_replacements)
    
    return host_tabulation_map    
    
    
def calculate_host_ecological_changes(host_records,abundance_matrix,speciess,samples):
    
    speciess_array = numpy.array(speciess)
    host_ecological_map = {}
    for host_record in host_records:
        cohort,subject,t0,t1 = host_record
        
        f0s = abundance_matrix[:,samples.index(t0)]
        f1s = abundance_matrix[:,samples.index(t1)]
        
        D = jensenshannon(f0s,f1s,base=config.default_entropy_base)
        alpha0 = entropy(f0s,base=config.default_entropy_base)
        alpha1 = entropy(f1s,base=config.default_entropy_base)
        extinctions = speciess_array[(f0s>config.default_extinction_f0)*(f1s<config.default_extinction_f1)]
        invasions = speciess_array[(f1s>config.default_extinction_f0)*(f0s<config.default_extinction_f1)]

        jsd, neff_delta, neff_avg, delta_eff = diversity_utils.calculate_neffs(f0s,f1s,delta_min=0)
        
        profile_jsds, delta_mins, fmins = diversity_utils.calculate_profile_jsds(f0s,f1s)
        
        #x_resolved_jsds, x_bins = diversity_utils.calculate_x_resolved_jsds(f0s,f1s)
        
        richness0 = ((1.0-numpy.exp(-f0s/config.default_richness_fmin))).sum()
        #print richness0
        
        host_ecological_map[host_record] = (alpha0,alpha1,jsd,(profile_jsds,delta_mins,fmins), richness0, neff_delta, neff_avg,delta_eff, extinctions,invasions)
    
    return host_ecological_map
    
def calculate_js_item_map(host_records, abundance_matrix,speciess,samples):
    
    speciess_array = numpy.array(speciess)
    js_item_map = {}
    for host_record in host_records:
        cohort,subject,t0,t1 = host_record
        
        f0s = abundance_matrix[:,samples.index(t0)]
        f1s = abundance_matrix[:,samples.index(t1)]
        
        js_items = diversity_utils.calculate_js_items(f0s,f1s)    
        js_item_map[host_record] = js_items
    
    return js_item_map

def get_frequency(abundance_matrix,speciess,samples,species,sample):
    #print species,sample,speciess.index(species),samples.index(sample)
    return abundance_matrix[speciess.index(species),samples.index(sample)]


def collate_species_events(within_host_events):
    
    species_event_map = {}
    
    for cohort,subject,t0,t1,species,event in within_host_events:
        if cohort!='hmp':
            continue
        
        if species not in species_event_map:
            species_event_map[species] = [0,0,0]
        
        species_event_map[species][event] += 1
        
    return species_event_map    

def prune_within_host_changes(within_host_changes):
    
    within_host_events = calculate_within_host_events(within_host_changes)
    host_event_map = collate_within_host_events(within_host_events)
    host_classification_map = calculate_host_classifications(within_host_events)
    # prunes list of within host changes so that you don't have overlapping samples
    import parse_data
    sample_metadata_map = parse_data.parse_sample_metadata_map()
    
    subject_classification_vector_map = {}
    for host_record in host_classification_map:
        cohort,subject,sample_i,sample_j = host_record
        
        subject_record = (cohort,subject)
        if subject_record not in subject_classification_vector_map:
            subject_classification_vector_map[subject_record] = [[("",""),("",""),("","")],[-1,-1,-1]]
            
        t_i = sample_metadata_map[sample_i][-1]
        t_j = sample_metadata_map[sample_j][-1]
        
        if t_i==1 and t_j==2:
            idx=0
        elif t_i==2 and t_j==3:
            idx=1
        elif t_i==1 and t_j==3:
            idx=2
        else:
            sys.stderr.write("Weird sample ordering: %s, %s. Shouldn't happen!\n" % (str(t_i),str(t_j)))
            continue
        
        subject_classification_vector_map[subject_record][0][idx] = (sample_i,sample_j)
        subject_classification_vector_map[subject_record][1][idx] = host_classification_map[host_record]
    
    # Done compiling triplet counts
    allowed_host_records = set()
    
    num_ones = 0
    num_twos = 0
    num_threes = 0
    three_zero_hosts = {}
    two_zero_hosts = {}
    
    
    num_shorts = 0
    num_longs = 0
    # First go through and look at the positives
    for subject_record in subject_classification_vector_map:
        cohort,subject = subject_record
        sample_pairs = numpy.array(subject_classification_vector_map[subject_record][0])
        classification_vector = numpy.array(subject_classification_vector_map[subject_record][1])
        
        num_masked = (classification_vector<-0.5).sum()
        
        if num_masked==3:
            sys.stderr.write("All masked.. shouldn't get here!\n")
        elif num_masked==2:
            # Just one timepoint pair, pretty easy
            good_idx = (classification_vector>-0.5)
            sample_i, sample_j = sample_pairs[good_idx][0]
            allowed_host_records.add((cohort,subject,sample_i,sample_j))
            
            if good_idx[2]:
                num_longs+=1
            else:
                num_shorts+=1
        
        elif num_masked==1:
            # one masked timepoint pair. need to make a choice 
            if classification_vector[2]==-1:
                # easy, take both
                sample_i, sample_j = sample_pairs[0]
                allowed_host_records.add((cohort,subject,sample_i,sample_j))
                sample_i, sample_j = sample_pairs[0]
                allowed_host_records.add((cohort,subject,sample_i,sample_j))
                num_shorts+=2
            else:
                if classification_vector[0]==-1:
                    good_idx = 1
                else:
                    good_idx = 0
                    
                #if classification_vector[good_idx]>0.5 or classification_vector[2]==0:
                    # Take the shorter one!
                #    sample_i, sample_j = sample_pairs[good_idx]
                #    allowed_host_records.add((cohort,subject,sample_i,sample_j))
                #    num_shorts+=1
                #else:
                    # take the longer one!
                #    sample_i, sample_j = sample_pairs[2]
                #    allowed_host_records.add((cohort,subject,sample_i,sample_j))
                #    num_longs+=1

                if classification_vector[2]>0.5 or classification_vector[good_idx]==0:
                    # take the longer one!
                    sample_i, sample_j = sample_pairs[2]
                    allowed_host_records.add((cohort,subject,sample_i,sample_j))
                    num_longs+=1
                    
                else:
                    # Take the shorter one!
                    sample_i, sample_j = sample_pairs[good_idx]
                    allowed_host_records.add((cohort,subject,sample_i,sample_j))
                    num_shorts+=1

        
        elif num_masked==0:
            # no masked timepoint pairs, now we really have to make a choice. 
            #if classification_vector[0]>0.5 or classification_vector[1]>0.5 or classification_vector[2]==0:
                # take both shorter ones!
            #    sample_i, sample_j = sample_pairs[0]
            #    allowed_host_records.add((cohort,subject,sample_i,sample_j))
            #    sample_i, sample_j = sample_pairs[1]
            #    allowed_host_records.add((cohort,subject,sample_i,sample_j))
            #    num_shorts+=2
            #else:
                # take the long one! 
            #    sample_i, sample_j = sample_pairs[2]
            #    allowed_host_records.add((cohort,subject,sample_i,sample_j))
            #    num_longs+=1
            
            # no masked timepoint pairs, now we really have to make a choice. 
            if classification_vector[2]>0.5 or (classification_vector[0]==0 and classification_vector[1]==0):
                # take the long one! 
                sample_i, sample_j = sample_pairs[2]
                allowed_host_records.add((cohort,subject,sample_i,sample_j))
                num_longs+=1
            else:
                # take both shorter ones!
                sample_i, sample_j = sample_pairs[0]
                allowed_host_records.add((cohort,subject,sample_i,sample_j))
                sample_i, sample_j = sample_pairs[1]
                allowed_host_records.add((cohort,subject,sample_i,sample_j))
                num_shorts+=2
    
    print("Finished with %d shorts sand %d longs" % (num_shorts,num_longs))
    
    new_within_host_changes = []
    for record_idx in range(0,len(within_host_changes)):
            
        cohort,subject,t0,t1,species,Lsnp,ksnp,Lgene,kgene = within_host_changes[record_idx]
        
        host_record = (cohort,subject,t0,t1)
        if host_record in allowed_host_records:
            new_within_host_changes.append(within_host_changes[record_idx])
    
    return new_within_host_changes

##########
#
# Do a logistic regression w/ entropy and phylum as input variables
#
##########
def calculate_matrix_from_lists(row_idxs,col_idxs):

    matrix = numpy.zeros((max(row_idxs)+1,max(col_idxs)+1))
    
    for row_idx, col_idx in zip(row_idxs,col_idxs):
        matrix[row_idx,col_idx]+=1
        
    return matrix
    
def calculate_loglikelihood(ns):
    
    pooled_ns = ns.sum(axis=0) # sums across species 
    pooled_ntot = pooled_ns.sum()
    pooled_ps = pooled_ns*1.0/pooled_ntot
    
    ntots = ns.sum(axis=1) # sums across event types
    ps = ns*1.0/ntots[:,None]
    
    return (ns*numpy.log(ps/pooled_ps[None,:]+(ns==0))).sum()

####
#
# Calculates "confidence intervals" on rate from Poisson distribution 
# based on n>=0 counts at L>=0 sites.
#
####
def calculate_poisson_confidence_interval(n,L,alpha=0.5): # by default use a 50% confidence interval
    
    if n<0.5:
        # No counts. Have some info on upper bound, but none on lower bound.
        plower = 0
        pupper = log(2/alpha)/L
    
    else:
        # Posterior distribution is Gamma with shape n-1 and scale 1/L 
        # Get confidence intervals from tail probabilities
        plower = gamma.ppf(alpha/2, n)/L
        pupper = gamma.ppf(1-alpha/2,n)/L
        
    return plower,pupper  
    
    
def coarse_grain_abandance_matrix(abundance_matrix,speciess,samples,taxonomic_level="family"):
    
    genus_family_map = bacterial_phylogeny_utils.get_genus_family_map()
    genus_phylum_map = bacterial_phylogeny_utils.get_genus_phylum_map()

    # can be 'all', 'species','genus' or phylum 
    
    taxon_abundance_map = {}
    
    for species_idx in range(0,len(speciess)):
        
        taxon = bacterial_phylogeny_utils.get_taxon_name(speciess[species_idx],taxonomic_level,genus_family_map,genus_phylum_map)
        if taxon not in taxon_abundance_map:
            taxon_abundance_map[taxon] = numpy.zeros(len(samples))*1.0
        
        taxon_abundance_map[taxon] += abundance_matrix[species_idx,:]
        
    taxa = []
    coarse_abundance_matrix = []
    for taxon in taxon_abundance_map:
        taxa.append(taxon)
        coarse_abundance_matrix.append(taxon_abundance_map[taxon])       
    
    coarse_abundance_matrix = numpy.array(coarse_abundance_matrix)
    taxa = numpy.array(taxa)
    
    return coarse_abundance_matrix, taxa
            
if __name__=='__main__':
    
    import parse_data
    sample_metadata_map = parse_data.parse_sample_metadata_map() 
    within_host_changes = parse_data.parse_within_host_changes(filename='within_host_changes_all_pairs.txt')
    abundance_matrix,speciess,samples = parse_data.parse_abundances()
    
    print("Loaded %d within host changes" % len(within_host_changes))
    within_host_events = calculate_within_host_events(within_host_changes)
    host_classification_map = calculate_host_classifications(within_host_events)
    
    pruned_within_host_changes = prune_within_host_changes(within_host_changes)
    
    print("Done! Left with %d within-host changes" % len(pruned_within_host_changes))
    
    #print "Boostrapping..."
    #bootstrapped_within_host_eventss = calculate_species_bootstrapped_within_host_events(within_host_events,num_bootstraps=1000)
    #print "Done!"
    
    
            
            
        