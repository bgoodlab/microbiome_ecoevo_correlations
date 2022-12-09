import numpy
import sys

def get_pretty_name(species_name):
    return " ".join(species_name.split("_")[:-1])
    
def get_genus_name(species_name):
    return species_name.split("_")[0]
    
def get_genus_family_map():

    genus_family_map = {}
    
    file = open("midas_genome_taxonomy.txt","r")
    file.readline() # header
    for line in file:
        items =  line.split("\t")
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
        items =  line.split("\t")
        
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
        sys.stderr.write("Couldn't find phylum name for %s\n" % species_name)
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
    