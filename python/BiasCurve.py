# -*- coding: utf-8 -*-
"""
Computes circular bias in responses towards for secondary orientation
here called 'distractors'

@author: Sammi, Jasper (Dec 2020)
"""
import numpy as np
import scipy as sp
import pandas as pd
import pycircstat as circstat

def circ_bini(x, nbin, pbin):
    
    quantbeg = np.mod(np.linspace(0-pbin/2, 1-1/nbin - pbin/2, nbin),1)
    quantend = np.mod(quantbeg+pbin,1)
    xbinbeg = np.quantile(x, quantbeg)
    xbinend = np.quantile(x, quantend)
    
    ibin = np.full(shape = (nbin, np.size(x)), fill_value = False)
    
    for i in range(nbin):
        if quantbeg[i] < quantend[i]: #no wrapping needed
            ibin[i,:] = np.logical_and(np.greater_equal(x,xbinbeg[i]), np.less_equal(x, xbinend[i]))
        else: #wrapping
            ibin[i,:] = np.logical_or(np.greater_equal(x, xbinbeg[i]), np.less_equal(x, xbinend[i]))
    
    return ibin

def calculate_bias(signed_error, target, other, subid, confidence = None):
    targori = np.radians(target)
    otherori = np.radians(other)
    
    angdiff = circstat.cdiff(otherori, targori)
    
    #demean error to remove any angular biases in overall responding (agnostic to any other data)
    signed_error = np.subtract(signed_error, signed_error.mean())
    
    nbin = 64 #64 bins
    pbin = 0.25 #1/4 fo the data per bin
    
    bins = circ_bini(angdiff, nbin, pbin)
    
    biases    = np.full(shape = nbin, fill_value = np.nan)
    bincentre = np.full(shape = nbin, fill_value = np.nan)
    precs     = np.full(shape = nbin, fill_value = np.nan)
    binid     = np.full(shape = nbin, fill_value = 0)
    
    for i in range(nbin): #loop over bins
        biases[i]     = np.mean(signed_error[bins[i,:]]) #sp.stats.circmean(signed_error[bins[i,:]], high = np.pi, low = -np.pi)
        bincentre[i]  = sp.stats.circmean(angdiff[bins[i,:]], high = np.pi, low = -np.pi)
        precs[i]      = sp.stats.circstd(signed_error[bins[i,:]], high = np.pi, low = -np.pi)
        binid[i]      = i
    
    df = pd.DataFrame()
    df['subid'] = np.full(shape = nbin, fill_value = subid)
    df['bincentre'] = bincentre
    df['bias'] = biases
    df['prec'] = precs
    df['binid'] = binid
    
    return df

def calculate_bias_covar(signed_error, target, other, subid, nbin = 64, pbin = 0.25, covariate = None, covariate_label = None, median_split = None):
    targori = np.radians(target)
    otherori = np.radians(other)
    
    angdiff = circstat.cdiff(otherori, targori)
    
    #demean error to remove any angular biases in overall responding (agnostic to any other data)
    signed_error = np.subtract(signed_error, signed_error.mean())
    
    # #params for creating data bins
    # nbin = 64 #64 bins
    # pbin = 0.25 #1/4 fo the data per bin
    
    if median_split == True and covariate is not None: #check there's a covar of interest and that we need to median split data for it
        covar_median = np.median(covariate)
        belowmed_trls = np.less_equal(covariate, covar_median)
        abovemed_trls = np.greater(covariate, covar_median)
        
        belowmed_vars = dict()
        abovemed_vars = dict()
        
        belowmed_vars['angdiff']        = angdiff[belowmed_trls]
        belowmed_vars['targori']        = targori[belowmed_trls]
        belowmed_vars['otherori']       = otherori[belowmed_trls]
        belowmed_vars['signed_error']   = signed_error[belowmed_trls]
        belowmed_vars[covariate_label]  = covariate[belowmed_trls]
        
        abovemed_vars['angdiff']        = angdiff[abovemed_trls]
        abovemed_vars['targori']        = targori[abovemed_trls]
        abovemed_vars['otherori']       = otherori[abovemed_trls]
        abovemed_vars['signed_error']   = signed_error[abovemed_trls]
        abovemed_vars[covariate_label]  = covariate[abovemed_trls]
        
        belowmed_bins = circ_bini(belowmed_vars['angdiff'], nbin, pbin)
        abovemed_bins = circ_bini(abovemed_vars['angdiff'], nbin, pbin)

        below_vars = dict()
        above_vars = dict()
        for key in ['biases', 'bincentre', 'precs', 'binid']:
            below_vars[key] = np.full(shape = nbin, fill_value = np.nan)
            above_vars[key] = np.full(shape = nbin, fill_value = np.nan)
        
        for i in range(nbin): #loop over bins
            below_vars['biases'][i]    = np.mean(belowmed_vars['signed_error'][belowmed_bins[i,:]])
            below_vars['bincentre'][i] = sp.stats.circmean( belowmed_vars['angdiff'][belowmed_bins[i,:]], high = np.pi, low = -np.pi )
            below_vars['precs'][i]     = sp.stats.circstd(  belowmed_vars['signed_error'][belowmed_bins[i,:]], high = np.pi, low = -np.pi )
            below_vars['binid'][i]     = i
    
            above_vars['biases'][i]    = np.mean(abovemed_vars['signed_error'][abovemed_bins[i,:]])
            above_vars['bincentre'][i] = sp.stats.circmean( abovemed_vars['angdiff'][abovemed_bins[i,:]], high = np.pi, low = -np.pi )
            above_vars['precs'][i]     = sp.stats.circstd(  abovemed_vars['signed_error'][abovemed_bins[i,:]], high = np.pi, low = -np.pi )
            above_vars['binid'][i]     = i
        
        above_df = pd.DataFrame(above_vars); above_df['median_split'] = 'above'; above_df.binid = above_df.binid.astype(int)
        below_df = pd.DataFrame(below_vars); below_df['median_split'] = 'below'; below_df.binid = below_df.binid.astype(int)
        
        df = pd.concat([above_df, below_df])
    else:
        bins = circ_bini(angdiff, nbin, pbin)
        
        biases    = np.full(shape = nbin, fill_value = np.nan)
        bincentre = np.full(shape = nbin, fill_value = np.nan)
        precs     = np.full(shape = nbin, fill_value = np.nan)
        binid     = np.full(shape = nbin, fill_value = 0)
        
        for i in range(nbin): #loop over bins
            biases[i]     = np.mean(signed_error[bins[i,:]]) #sp.stats.circmean(signed_error[bins[i,:]], high = np.pi, low = -np.pi)
            bincentre[i]  = sp.stats.circmean(angdiff[bins[i,:]], high = np.pi, low = -np.pi)
            precs[i]      = sp.stats.circstd(signed_error[bins[i,:]], high = np.pi, low = -np.pi)
            binid[i]      = i
        
        df = pd.DataFrame()
        df['subid'] = np.full(shape = nbin, fill_value = subid)
        df['bincentre'] = bincentre
        df['bias'] = biases
        df['prec'] = precs
        df['binid'] = binid
    
    df['subid'] = subid
    return df
