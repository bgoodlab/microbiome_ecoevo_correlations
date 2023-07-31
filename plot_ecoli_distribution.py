import pylab
import numpy
import parse_data
import stats_utils

abundance_matrix,speciess,samples = parse_data.parse_abundances()

for i in range(0,len(speciess)):
	
	if speciess[i].startswith('Escherichia_coli'):
		
		
		abundances = abundance_matrix[i,:]
		
		print(speciess[i],numpy.median(abundances[abundances>1e-05]))
		xs,Ss = stats_utils.calculate_unnormalized_survival_from_vector(abundances)
		
		pylab.semilogx(xs,Ss/Ss[0])
		
pylab.savefig('ecoli_distribution.pdf',bbox_inches='tight')