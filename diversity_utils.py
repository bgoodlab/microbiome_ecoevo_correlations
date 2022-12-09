import numpy
import pylab
from scipy.special import xlogy 
from math import log,log10

def calculate_jsd_precursors(f0s,f1s):
    
    favgs = (f0s+f1s)/2
    dfs = (f1s-f0s)
    xs = dfs/2/favgs
    
    return favgs,xs
    
def calculate_jsd_stats(f0s,f1s,neff_percentage=0.5):
    
    favgs,xs = calculate_jsd_precursors(f0s,f1s)
    
    deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    
    jsd = (deltas*favgs).sum()
    
    sorted_idxs = numpy.flip(numpy.argsort(deltas))
    
    neff = 1+ numpy.nonzero(numpy.cumsum(deltas(sorted_idxs))>(neff_percentage*total_jsd))[0]
    
    return total_jsd**0.5, neff
    
def calculate_neffs(f0s,f1s,neff_percentage=0.9,delta_min=0):
    
    favgs = (f0s+f1s)/2
    dfs = (f1s-f0s)
    xs = dfs/2/(favgs+(favgs==0))
    deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    deltas = deltas*(deltas>=delta_min)
    
    js_items = deltas*favgs
    jsd = (js_items).sum()
    
    delta_eff = (js_items*deltas).sum()/(js_items.sum())
    sorted_idxs = numpy.flip(numpy.argsort(js_items))
    
    neff_delta = 1+numpy.nonzero(numpy.cumsum(js_items[sorted_idxs])>(neff_percentage*jsd))[0][0]
    
    sorted_idxs = numpy.flip(numpy.argsort(favgs))
    neff_avg = 1+numpy.nonzero(numpy.cumsum(favgs[sorted_idxs])>(neff_percentage))[0][0]
    
    return jsd**0.5, neff_delta, neff_avg, delta_eff

def calculate_js_items(f0s,f1s):
    
    favgs = (f0s+f1s)/2
    dfs = (f1s-f0s)
    xs = dfs/2/(favgs+(favgs==0))
    deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    js_items = favgs*deltas
    
    return js_items
        
def calculate_profile_jsds(f0s,f1s):
    
    #foldchange_mins = numpy.hstack([[1],numpy.logspace(log10(2),3,20)])
    foldchange_mins = numpy.logspace(0,2.5,20)
    xs = (foldchange_mins-1)/(foldchange_mins+1)
    delta_mins = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    
    #delta_mins = numpy.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    favg_mins = numpy.array([0,1e-03])
    favgs = (f0s+f1s)/2
    dfs = (f1s-f0s)
    xs = dfs/2/(favgs+(favgs==0))
    deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    
    profile_jsds_all = numpy.sqrt(((deltas*favgs)[:,None,None]*(deltas[:,None,None]>=delta_mins[None,:,None])*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0)/(favgs[:,None,None]*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0))
     
    profile_jsds_positive = numpy.sqrt(((deltas*favgs)[:,None,None]*(deltas[:,None,None]>=delta_mins[None,:,None])*(xs[:,None,None]>0)*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0)/(favgs[:,None,None]*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0))
    
    profile_jsds_negative = numpy.sqrt(((deltas*favgs)[:,None,None]*(deltas[:,None,None]>=delta_mins[None,:,None])*(xs[:,None,None]<0)*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0)/(favgs[:,None,None]*(favgs[:,None,None]>=favg_mins[None,None,:])).sum(axis=0))
    
    profile_jsds = numpy.array([profile_jsds_all, profile_jsds_positive, profile_jsds_negative]) 
    
    return profile_jsds, numpy.log10(foldchange_mins), favg_mins
    
def calculate_x_resolved_jsds(f0s,f1s):
    
    foldchange_bins = numpy.logspace(-2,2,30)
    x_bins = (foldchange_bins-1)/(foldchange_bins+1)
    x_bins[0]=-1
    x_bins[-1]=1+1e-06
    
    favgs = (f0s+f1s)/2
    dfs = (f1s-f0s)
    xs = dfs/2/(favgs+(favgs==0))
    deltas = (xlogy(1-xs,1-xs)+xlogy(1+xs,1+xs))/(2*log(2))
    
    x_resolved_jsds = ((deltas*favgs)[:,None]*(xs[:,None]>=(x_bins[:-1][None,:]))*(xs[:,None]<(x_bins[1:][None,:]))).sum(axis=0)
    
    return x_resolved_jsds, numpy.log10(foldchange_bins)
    