import numpy
import config
import sys

# identified to have swapped samples (Thanks Zhiru Liu)
BAD_SUBJECTS = set(['763536994','763880905'])

def parse_sample_metadata_map():
    
    sample_metadata_map = {}
    
    # First load HMP metadata
    #file = open(config.scripts_directory+"HMP_ids_order.txt","r")
    file = open(config.default_metadata_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        merged_sample_ids = [subitem.strip() for subitem in items[2].split(",")]
        accession_ids = [subitem.strip() for subitem in items[3].split(",")]
        country = items[4].strip()
        continent = items[5].strip()
        order = int(items[6].strip())
        
        sample_metadata_map[sample_id] = (subject_id, sample_id, accession_ids, country, continent, order)
        
    file.close()
    
    return sample_metadata_map
	
def parse_within_host_changes(filename=config.default_within_host_filename,allowed_cohorts=["hmp"],remove_bad_subjects=True):
	
	
	within_host_changes = []
	
	sample_metadata_map = parse_sample_metadata_map()
	file = open(filename,"r")
	file.readline() # header
	for line in file:
		items = line.split()
		cohort = items[0].strip()
		sample_t0 = items[1].strip()
		sample_t1 = items[2].strip()
		species = items[3].strip()
		Lsnp = int(items[4])
		ksnp = int(items[5])
		Lgene = int(items[6])
		kgene = int(items[7])
		
		Lprivate = int(items[8])
		kprivate = int(items[9])
		
		if cohort not in allowed_cohorts:
			continue
			
		if sample_t0 not in sample_metadata_map:
			sys.stderr.write("Sample not in metadata map. Shouldn't happen!\n")
			continue
			
		subject = sample_metadata_map[sample_t0][0]
		
		if remove_bad_subjects and (subject in BAD_SUBJECTS):
			#print "removing subject", subject
			continue
			
		within_host_changes.append((cohort, subject,sample_t0,sample_t1,species,Lsnp,ksnp,Lprivate,kprivate))
		
		
	
	file.close()
	return within_host_changes 

#def parse_abundances(filename="coverage.txt"):
def parse_abundances(filename=config.default_coverage_filename):	
	file = open(filename,"r")
	header = file.readline() 
	samples = []
	for sample in header.split()[1:]:
		if sample[-1]=='c':
			samples.append(sample[:-1])
		else:
			samples.append(sample)
	
	#samples = numpy.array(samples)
	coverage_matrix = []
	speciess = []
	for line in file:
		items = line.split()
		
		species = items[0]
		speciess.append(species)
		coverages = numpy.array([float(item) for item in items[1:]])
		
		#if species.startswith('Odoribacter_laneus_62216'):
		#	 print coverages[coverages>10]
		#	 print samples[coverages>10]
		coverage_matrix.append(coverages)
		
	coverage_matrix = numpy.array(coverage_matrix)
	speciess = numpy.array(speciess)
	abundance_matrix = coverage_matrix*1.0/coverage_matrix.sum(axis=0)
	#abundance_matrix = coverage_matrix
	
	
	# filter for things that are present anywhere
	good_species_idxs = ((abundance_matrix>1e-04).sum(axis=1)>0)
	
	abundance_matrix = abundance_matrix[good_species_idxs,:]
	speciess = list(speciess[good_species_idxs])
	return abundance_matrix,speciess,samples		


if __name__=='__main__':
	
	# Test these things
	sample_metadata_map = parse_sample_metadata_map()
	
	within_host_changes_consecutive = parse_within_host_changes()
	within_host_changes_all = parse_within_host_changes(filename=config.default_within_host_all_pairs_filename,remove_bad_subjects=True)
	
	abundance_matrix,speciess,samples = parse_abundances()
	
	metadata_samples = set()
	metadata_subjects = set()
	
	for sample_id in sample_metadata_map:
		
		subject_id, sample_id, accession_id, country, continent, order = sample_metadata_map[sample_id]
		
		metadata_subjects.add(subject_id)
		metadata_samples.add(sample_id)
		
		
	output_file = open("figures/figure_0_output.txt","w")
	output_file.write("%d total samples from %d total subjects\n" % (len(metadata_samples),len(metadata_subjects)))
	
	for type,within_host_changes in zip(['all timepoint pairs','consecutive timepoints'],[within_host_changes_all,within_host_changes_consecutive]):
		
		num_replacements=0
		num_modifications=0
		within_host_samples = set()
		for idx in range(0,len(within_host_changes)):
			cohort, subject,sample_t0,sample_t1,species,Lsnp,ksnp,Lprivate,kprivate = within_host_changes[idx]
		
			within_host_samples.add(sample_t0)
			within_host_samples.add(sample_t1)
			if ksnp > config.default_replacement_threshold:
				num_replacements+=1
			elif ksnp > 0:
				num_modifications+=1
		output_file.write("%s:\n" % type)
		output_file.write("%d total modification events\n" % num_modifications)		
		output_file.write("%d total replacement events\n" % num_replacements) 
  
	output_file.close()
	
	
	