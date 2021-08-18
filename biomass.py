import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy import integrate

# These functions are needed to express the sum over R and nR analytically.
def R(n,e,l,s):
    '''Unnormalized struture function for METE.
    n,e are microscopic variables.
    l are lambdas
    s are state variables, call S, N, or E'''
    return np.exp(-l[0]*n-l[1]*n*e)

def Rsum(e,l,s):
    '''Unnormalized struture function for METE, summed over n.
    e is a microscopic variable.
    l are lambdas
    s are state variables, call S, N, or E'''
    # Define exponent for lambdas 1-2
    l12 = l[0]+l[1]*e
    num = np.exp(-l12) - np.exp(-l12*(s['N']+1))
    denom = 1 - np.exp(-l12)
    return num/denom

def nRsum(e,l,s):
    '''Unnormalized struture function for METE multiplied by n, then summed over n. 
    e is a microscopic variable.
    l are lambdas
    s are state variables, call S, N, or E'''
    # Define exponent for lambdas 1-2, with n=1 since we've done sum over n already.
    l12 = l[0]+l[1]*e
    # Split up to make it easier
    num = np.exp(-l12)-(s['N']+1)*np.exp(-l12*(s['N']+1))+s['N']*np.exp(-l12*(s['N']+2))
    denom = (1-np.exp(-l12))**2
    return num/denom

def zfun(l,s):
    '''Calculate partition function for DynaMETE. 
    This function uses the analytic sum from Rsum and single integral over e.
    The integral is done with quad over log e.
    l are lambdas
    s are state variables, call S, N, or E'''
    # Return only value not error
    return integrate.quad(lambda loge: np.exp(loge)*Rsum(np.exp(loge),l,s),0,np.log(s['E']))[0]

def constraints(l,s):
    '''This function takes in lambdas and state variables and tests the constraints numerically.
    The return is an array with the percent deviation from the exact constraint.'''
    # First get the normalization
    z = zfun(l,s)
    # Now we need the sum over n and ne
    nR = integrate.quad(lambda loge: np.exp(loge)*nRsum(np.exp(loge),l,s),0,np.log(s['E']))[0]
    neR = integrate.quad(lambda loge: np.exp(loge*2)*nRsum(np.exp(loge),l,s),0,np.log(s['E']))[0]
    # Normalize
    nR /= z
    neR /= z
    # Now get the exact values
    ns = s['N']/s['S']
    es = s['E']/s['S']
    # Now return percent difference
    return np.array([(nR-ns)/ns,(neR-es)/es])


# METE functions
def beta_constraint(b,s):
    '''This is the beta constraint in METE with give state variables. Use this as a function call to get beta.
    Inputs s as state variables, call S, N, or E
    Also inputs beta
    outputs beta constraint to minimize'''
    # Generate n array up to N
    narr = np.arange(s['N'])+1
    # This is the better constraint condition to use
    return np.exp(-narr*b).sum()/(np.exp(-narr*b)/narr).sum() - s['N']/s['S']
    # Old constraint
    #return b*np.log(1/(1-np.exp(-b)))-s['S']/s['N'] #Old

def mete_lambdas(s,b0=0.0001,thresh=0.05):
    '''This returns the METE lambdas for a given set of state variables.
    Inputs s as state variables, call S, N, or E
    Optional input of an initial beta, if we know it's going to be somewhere other than small positive.
    outputs array of lambdas
    Also checks to see if conditions are actually satisfied, and outputs a warning if the percent difference
    of the result is greater than 5% (by default)
    This can be modified to exclude this test by setting thresh to False, or changing it to a difference percent'''
    beta = fsolve(beta_constraint,b0,args=s)[0]
    l2 = s['S']/(s['E']-s['N']) 
    ls = np.array([beta-l2,l2]) # Only 2 compared to 5 from dynamete
    if thresh:
        con = constraints(ls,s)
        if np.any(con > thresh):
            print("Constraints are not satisfied within {:.0f}%".format(thresh*100))
    return ls


# Define biomass equations
def biomass(s,power=4/3):
    '''Take in state variables and numerically calculate the biomass.
    This is the equivalent of B = S Sum_n n e^(power) R 
    This function uses th exact sum over n and then does the integral over e^(power)
    The power is by default 4/3, but can be changed if we want to test other scaling relationships
    The lambdas are calculated using the function mete_lambdas above, which are approximate, 
    but should be good for most realistic state variables. The state variables should be given as a pandas Series
    s are state variables, call S, N, or E
    power is by default 4/3, but can be changed to test other scaling relationships
    '''
    # get METE lambdas
    l = mete_lambdas(s)
    # Get normalization
    z = zfun(l,s)
    # Do the required integral in log space
    b = s['S']*integrate.quad(lambda loge: np.exp((power+1)*loge)*nRsum(np.exp(loge),l,s),0,np.log(s['E']))[0]
    # Normalize and return
    return b/z

def biomass_approx(s,order=0):
    '''The approximation for the biomass equation with power=4/3
    Here s is a panda series of state variables, call S, N, or E
    order indicates whether or not to include a first order correction. The default is to not, thus order=0.
    Change order=1 if you want to include the correction.'''
    # Get mete lambdas
    l = mete_lambdas(s)
    b = 4.17*s['E']**(4/3)/(s['S']**(1/3)*np.log(1/l.sum()))
    if order == 1:
        b *= (1-1.16*l.sum()**(1/3))
    return b