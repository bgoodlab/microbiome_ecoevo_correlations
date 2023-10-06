import parse_data
import numpy
import config
import sys
	
def parse_within_host_changes(filename=config.default_within_host_filename,allowed_cohorts=["hmp"]):
	
	
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
		
		if subject in BAD_SUBJECTS:
			#print "removing subject", subject
			continue
			
		within_host_changes.append((cohort, subject,sample_t0,sample_t1,species,Lsnp,ksnp,Lprivate,kprivate))
		
		
	
	file.close()
	return within_host_changes 

#def parse_abundances(filename="coverage.txt"):
def parse_abundances(filename="species_coverage.txt"):
	
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
	
	# Load up the within host changes
	within_host_changes = parse_data.parse_within_host_changes()
	within_host_changes_all = parse_data.parse_within_host_changes(filename=config.default_within_host_all_pairs_filename, remove_bad_subjects=False)
	
	print("Within host samples:", len(within_host_changes),len(within_host_changes_all),'(all)')
	
	# Figure out which samples are in the set we're looking at
	within_host_samples = set()
	for idx in range(0,len(within_host_changes_all)):
		cohort, subject,sample_t0,sample_t1,species,Lsnp,ksnp,Lprivate,kprivate = within_host_changes_all[idx]
		
		within_host_samples.add(sample_t0)
		within_host_samples.add(sample_t1)
	
	print("Num within host samples", len(within_host_samples))	
	# Now restrict metadata file to ones that we're actually looking at
	full_metadata_file = open("HMP1-2_ids_order.txt","r")
	line = full_metadata_file.readline() # header
	
	sample_record_map = {}
	record_data_map = {}
	for line in full_metadata_file:
		
		items = line.split("\t")
		subject_id = items[0].strip()
		sample_id = items[1].strip()
		accession = items[2].strip()
		country = items[3].strip()
		continent = items[4].strip()
		visno = items[5].strip()
		
		record = (subject_id,country,continent,visno)
		
		sample_record_map[sample_id] = record
		
		if record not in record_data_map:
			record_data_map[record]={'sample_ids':[], 'accessions':[]}
		
		record_data_map[record]['sample_ids'].append(sample_id)
		record_data_map[record]['accessions'].append(accession)
		
	record_output_map = {}
	for sample_id in within_host_samples:
		record = sample_record_map[sample_id]
		subject_id, country, continent, visno = record	
		
		sample_ids = record_data_map[record]['sample_ids']
		accessions = record_data_map[record]['accessions']
		
		output_str = "\t".join([subject_id, sample_id, ",".join(sample_ids),",".join(accessions),country, continent, visno])
		
		record_output_map[record] = output_str
	
	output_metadata_file = open(config.default_metadata_filename,"w")
	header = "\t".join(["subject_id","sample_id","merged_sample_ids","run_accessions","country","continent","VISNO"])
	output_metadata_file.write(header)
	output_metadata_file.write("\n")
	
	for record in sorted(record_output_map):
		output_metadata_file.write(record_output_map[record])
		output_metadata_file.write("\n")
	output_metadata_file.close()
	
	
	# Now do the same thing with the species coverage file
	full_coverage_file = open("garud_good_etal_species_coverage.txt","r")
	coverage_file = open(config.default_coverage_filename,"w")
	line = full_coverage_file.readline()
	items = line.split()
	
	relabeled_items = []
	for sample_id in items:
		if sample_id[-1]=='c':
			sample_id = sample_id[:-1]
		relabeled_items.append(sample_id)
	
	abundance_samples = set(relabeled_items[1:])
	
	good_idxs = [0]
	for i in range(1,len(relabeled_items)):
		sample_id = relabeled_items[i]
		if sample_id in within_host_samples:
			good_idxs.append(i)
	
	relabeled_items = numpy.array(relabeled_items)
	good_idxs = numpy.array(good_idxs)
	
	print("Num abundance samples,", len(abundance_samples))
	print("In evolution not in abundance", within_host_samples-abundance_samples)
			
	print("Number items", len(relabeled_items))
	print("Number good idxs", len(good_idxs))
	
	# Write header
	coverage_file.write("\t".join(relabeled_items[good_idxs]))
	coverage_file.write("\n")
	# Write everything else
	for line in full_coverage_file:
		items = numpy.array(line.split())
		coverage_file.write("\t".join(items[good_idxs]))
		coverage_file.write("\n")
	
	coverage_file.close()