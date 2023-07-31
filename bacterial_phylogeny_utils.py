import numpy
import sys
from Bio import Phylo

def get_pretty_name(species_name):
	return " ".join(species_name.split("_")[:-1])
	
def get_genus_name(species_name):
	return species_name.split("_")[0]
	
def get_genus_family_map():

	genus_family_map = {}
	
	file = open("midas_genome_taxonomy.txt","r")
	file.readline() # header
	for line in file:
		items =	 line.split("\t")
		phylum = items[4].strip()
		family = items[7].strip()
		genus = items[8].strip()
		species = items[9].strip()
		
		if genus=='' and species=='':
			continue
			
		if genus=='':
			genus = species.split()[0]
		
		if family=='':
			family = 'Other %s' % phylum
			
		#if family=='Ruminococcaceae':
		#	print(species, genus, family)
		if family=='Oscillospiraceae':
			print(species, genus, family)
		
		#if family=='Eubacteriaceae':
		#	print(species, genus, family)
		#if genus=='Eubacterium':
		#	print(species, genus, family)
		
		if genus in genus_family_map:
			if genus_family_map[genus] != family:
				print("Genus family conflict!", genus, family, genus_family_map[genus])
				
		genus_family_map[genus] = family
		
	return genus_family_map
	
def get_family_name(species_name, genus_family_map={}):
	
	if len(genus_family_map)==0:
		genus_family_map = get_genus_family_map()
		
	genus_name = get_genus_name(species_name)
	if genus_name in genus_family_map:
		return genus_family_map[genus_name]
	else:
		return genus_name
		

def get_genus_phylum_map():

	genus_phylum_map = {}
	
	file = open("midas_genome_taxonomy.txt","r")
	file.readline() # header
	for line in file:
		items =	 line.split("\t")
		
		phylum = items[4].strip()
		genus = items[8].strip()
		species = items[9].strip()
		
		if genus=='' and species=='':
			continue
			
		if genus=='':
			genus = species.split()[0]
		
		genus_phylum_map[genus] = phylum
	
	genus_phylum_map['Bacteroidales'] = 'Bacteroidetes' 
	genus_phylum_map['Burkholderiales'] = 'Proteobacteria'
	   
	return genus_phylum_map
	
def get_phylum_name(species_name, genus_phylum_map={}):
	
	if len(genus_phylum_map)==0:
		genus_phylum_map = get_genus_phylum_map()
		
	genus_name = get_genus_name(species_name)
	if genus_name in genus_phylum_map:
		return genus_phylum_map[genus_name]
	else:
		if genus_name=='Guyana':
			return 'TODO'
		elif genus_name== 'Xylanibacter':
			return 'TODO'
		elif genus_name== 'Jeddahella':
			return 'TODO'
		#sys.stderr.write("Couldn't find phylum name for %s\n" % species_name)
		return 'Unknown Bacteria' #genus_name
		
def get_taxon_name(species_name,taxonomic_level='species',genus_family_map={},genus_phylum_map={}):
	
	if taxonomic_level=='species':
		return species_name
	elif taxonomic_level=='genus':
		return get_genus_name(species_name)
	elif taxonomic_level=='family':
		return get_family_name(species_name,genus_family_map)
	elif taxonomic_level=='phylum':
		return get_phylum_name(species_name,genus_phylum_map)
	else:
		print("Wrong taxonomic level!", taxonomic_level)
		return species_name
		

def get_idxs_from_same_taxon(focal_species,speciess,taxonomic_level='species',genus_family_map={},genus_phylum_map={}):

	taxon_names = []
	for species_name in speciess:
		
		taxon = get_taxon_name(species_name,taxonomic_level,genus_family_map,genus_phylum_map)
		taxon_names.append(taxon)
	
	taxon_names = numpy.array(taxon_names)
			
	return (taxon_names==taxon_names[speciess.index(focal_species)])
	
def calculate_taxonomic_distances_from_focal_species(focal_species, speciess,genus_family_map={},genus_phylum_map={}):

	focal_genus = get_taxon_name(focal_species, 'genus',genus_family_map,genus_phylum_map)
	focal_family = get_taxon_name(focal_species, 'family',genus_family_map,genus_phylum_map)
	focal_phylum = get_taxon_name(focal_species, 'phylum',genus_family_map,genus_phylum_map)
	
	distances = []
	# 0 = same species
	# 1 = same genus
	# 2 = same family
	# 3 = same phyla
	# 4 = same kingdom
	for species_name in speciess:
		
		genus = get_taxon_name(species_name,'genus',genus_family_map,genus_phylum_map)
		family = get_taxon_name(species_name,'family',genus_family_map,genus_phylum_map)
		phylum = get_taxon_name(species_name,'phylum',genus_family_map,genus_phylum_map)
	 
		if species_name==focal_species:
			distance=0
		elif genus==focal_genus:
			distance=1
		elif family==focal_family:
			distance=2
		elif phylum==focal_phylum:
			distance=3
		else:
			distance=4
			
		distances.append(distance) 
	
	return numpy.array(distances)

def calculate_taxonomic_distance_matrix(speciess, genus_family_map={},genus_phylum_map={}):
	
	distance_matrix = []
	for focal_species in speciess:
		distances = calculate_taxonomic_distances_from_focal_species(focal_species, speciess,genus_family_map,genus_phylum_map)
		distance_matrix.append(distances)
		
	return numpy.array(distance_matrix)

def calculate_phylogenetic_distance_matrix(speciess):
	
		
	species_tree = Phylo.read("midas_species_tree.newick",'newick')

	terminals = species_tree.get_terminals() 

	species_idxs = []
	taxon_ids = []
	for species_name in speciess:
		taxon_id = species_name.split("_")[-1]
		taxon_ids.append(taxon_id)
		
		#if taxon_id not in terminals:
		#	print("ERROR: TAXON NOT IN TERMINALS", species_name, taxon_id)
			
		#idxs = terminals.index(taxon_id)
		#species_idxs.append(idx)

	distance_matrix = numpy.zeros((len(speciess),len(speciess)))*1.0
	
	for i in range(0,len(speciess)):
		
		if i%10==0:
			print(i)
		
		#x = terminals[species_idxs[i]]
		x = taxon_ids[i]
		
		for j in range(i+1,len(speciess)):
			#y = terminals[species_idxs[j]]
			y = taxon_ids[j]
			
			d = species_tree.distance(x, y)
			distance_matrix[i,j] = d
			distance_matrix[j,i] = d
	
	return distance_matrix

if __name__=='__main__':
	genus_family_map = get_genus_family_map()
	genus_phylum_map = get_genus_phylum_map()
	
	speciess = ['Bacteroides_uniformis_57318', 'Bacteroides_vulgatus_57955', 'Alistipes_putredinis_61533', 'Eubacterium_rectale_56927','Eubacterium_eligens_61678']
	
	taxonomic_distance_matrix = calculate_taxonomic_distance_matrix(speciess, genus_family_map, genus_phylum_map)
	print(taxonomic_distance_matrix)
	
	
	phylogenetic_distance_matrix = calculate_phylogenetic_distance_matrix(speciess)
	print(phylogenetic_distance_matrix)
	
	
	
	
	