#!/usr/bin/env python
# coding: utf-8

# This file will make a plot of E/N vs N/S, with colour indicating how good the analytic approximation is for the biomass expression. I will do this for both the first and second order approximations.
# 
# The equation being approximated is
# $$ B = S \int \sum n \varepsilon R$$
# 
# and the approximation is
# $$ B = 4.17 \frac{E^{4/3}}{S^{1/3} \log(1/\beta)} \left( 1-1.16 \beta^{1/3} \right) $$
# 
# Note that at this point, this includes the approximation for the normalization Z, so there are multiple approximations going on here.
# 
# Additionally, note that because we want to only have two variables, we are redefining N/S, E/S as just N and E, and B becomes B/S. In this case, the S dependence actually just drops out. So we can pick an S and plot accordingly. 
# 
# We here choose $S=50$. Initially, I picked $S=1$, however this doesn't work very well because the sum for the first constraint goes to N, but it is equal to N/S. So we want $S>1$ to make these not the same. So that this error in summation doesn't have a large effect, we need $S$ $\approx$ 30 at the smallest, but closer to $S=50$. Data points with $S$ less than this are marked with an asterix.

# In[4]:


import pandas as pd
import numpy as np
import biomass as bm

import matplotlib.pyplot as plt
# In[5]:


# Read in the data to overlay
data = pd.read_csv('data_statevariables.csv')


# In[6]:


# Define ns and en along with the number in each direction.
# Could also make this higher and then save the grid, but this is OK for now.
num_ns = 40
num_en = 40
# Log spacing since we will use a log log plot
ns = np.logspace(np.log10(4),np.log10(800),num=num_ns)
en = np.logspace(np.log10(4),np.log10(3500),num=num_en)

# Set s also
s0=50 


# In[8]:


# Numerical
biomass_num = np.zeros([num_en,num_ns]) 
# First order
biomass_1 = np.zeros([num_en,num_ns]) 
# Second order
biomass_2 = np.zeros([num_en,num_ns])

# This is slow but should be ok just for a plot
for i,e in enumerate(en):
    for j,n in enumerate(ns):
        s = pd.Series([s0,n*s0,e*n*s0],index=['S','N','E'])
        print(n,e) #just to keep track
        biomass_num[i,j] = bm.biomass(s)
        biomass_1[i,j] = bm.biomass_approx(s,order=0)
        biomass_2[i,j] = bm.biomass_approx(s,order=1)


# In[9]:


# Get percent differences. Note that in doing this, the S dependence drops out. Since really we care about B/S
biomass_1_percent = (biomass_1-biomass_num)/biomass_num
biomass_2_percent = (biomass_2-biomass_num)/biomass_num
# ignore the nans as those will be for the values where we didn't calculate the biomass

# Get max deviation in both cases to fix colorplot
max_dev_1 = np.nanmax(np.abs(biomass_1_percent))
max_dev_2 = np.nanmax(np.abs(biomass_2_percent))


# In[10]:


# Get percent values to plot
z_bm1 = -np.log10(np.abs(biomass_1_percent))
z_bm2 = -np.log10(np.abs(biomass_2_percent))


# In[11]:


# Plotting
fig,axs = plt.subplots(1,2,sharey=True,figsize=(8,3))

# Set contour range to go from the smallest difference to the biggest, in either data set
# Need one consistent range so colorscheme is accurate
contour_floor =  np.min([np.floor(np.nanmin(z_bm1)),np.floor(np.nanmin(z_bm2))])
contour_ceil =  np.max([np.ceil(np.nanmax(z_bm1)),np.ceil(np.nanmax(z_bm2))])
contour = np.arange(contour_floor,contour_ceil+1)

# Plot the contours
im1 = axs[0].contourf(np.log10(ns),np.log10(en),z_bm1,contour,cmap='winter')
im2 = axs[1].contourf(np.log10(ns),np.log10(en),z_bm2,contour,cmap='winter')

# Overlay the data
for ax in axs:
    # First separate it out based on a threshold for S. Basically with S<50 put an asterix that
    # this may not be the most accurate
    # Define the threshold
    st = 50
    # Get the indices where this is true
    ind_st = data['S'] > st
    # Plot these ones normally
    ax.scatter(np.log10(data[ind_st]['N']/data[ind_st]['S']),np.log10(data[ind_st]['E']/data[ind_st]['N']),
               c='0.9',alpha=0.8)
    # Plot the ones where S<50 with an asterix
    ax.scatter(np.log10(data[~ind_st]['N']/data[~ind_st]['S']),np.log10(data[~ind_st]['E']/data[~ind_st]['N']),
               marker='v',c='0.9',alpha=0.8)

# Adjust and add colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
fig.colorbar(im2, cax=cbar_ax)

# Now labels and titles
axs[0].set_ylabel(r'$\log_{10}(E/N)$')
for ax in axs:
    ax.set_xlabel(r'$\log_{10}(N/S)$')
axs[0].set_title('Zeroth order')
axs[1].set_title('First order')

# Save this as just logcontours.pdf to differentiate it from the other files. 
# It should be the same as some of the other figures though.
fig.savefig('Figures/FigS1.pdf',bbox_inches='tight')


# In[ ]:




