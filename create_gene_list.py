import config
import parse_data
import gzip

# First load the list of all within-host changes we are looking at
within_host_changes_all = parse_data.parse_within_host_changes(filename=config.default_within_host_all_pairs_filename, remove_bad_subjects=False)

species_sample_pair_map = {}

for idx in range(0,len(within_host_changes_all)):
	cohort, subject,sample_t0,sample_t1,species,Lsnp,ksnp,Lprivate,kprivate = within_host_changes_all[idx]
	
	if ksnp==0:
		continue
		
	if species not in species_sample_pair_map:
		species_sample_pair_map[species] = []
	
	species_sample_pair_map[species].append( (species, sample_t0, sample_t1) )

output_file = gzip.open("within_host_changes_all_genes.txt.gz","wt")

# Write header
output_file.write("\t".join(["Species","Sample_t0","Sample_t1","SNVs [gene_name;contig;position;var_type;A0;D0;A1;D1, ...]"]))
output_file.write("\n")

# Now go load the stuff from the Garud & Good et al paper
for species in species_sample_pair_map:

	species_sample_pair_map[species] = set(species_sample_pair_map[species])
	
	file = gzip.open("garud_good_etal_within_host_changes/%s.txt.gz" % species, 'rt')
	line = file.readline() # header
	
	for line in file:
		items = line.split(",")
		species = items[0].strip()
		sample_t0 = items[1].strip()
		sample_t1 = items[2].strip()
		type = items[3].strip()
		L = float(items[4])
		Perr = float(items[5])
		mutation_list = items[6:]
		
		if type!='snps':
			continue
			
		record = (species,sample_t0,sample_t1)
		
		#print(record)
		
		if record in species_sample_pair_map[species]:
			# We want to print one of these
			
			mutation_str = ", ".join(mutation_list)
			output_str = "\t".join([species,sample_t0,sample_t1,mutation_str])
			
			output_file.write(output_str)
			#output_file.write("\n")
			
	file.close()

output_file.close()