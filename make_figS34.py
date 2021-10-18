#!/usr/bin/env python
# coding: utf-8

# This file will make Figure S3 and S4 for the appendix comparing alternative scaling laws using state variable data from data_statevariables.txt and the biomass prediction code biomass.py

# In[6]:


import numpy as np
import pandas as pd
import biomass as bm
from scipy.stats import linregress
import matplotlib.pyplot as plt

# Import the data
data = pd.read_csv('data_statevariables.csv')
# Now add a column for predicted numerical biomass data
data['pBnum_2_3'] = np.zeros(len(data))
# Iterate through each row and append the biomass information
for index, row in data.iterrows():
    # Get 2/3
    dtemp = {'S': row['S'], 'N': row['N'], 'E': row['E_2_3']}
    data.loc[index,'pBnum_2_3'] = bm.biomass(dtemp,power=3/2)



# Choose color scheme
cm = 'winter'#'viridis'
# Get max and min species richness for colour scheme
smin = np.min(data['S'])
smax = np.max(data['S'])
# Set up normalization
norm = plt.Normalize(np.log(smin),np.log(smax))
# Get list of site types
stype = data['Type'].unique()
# Make marker list. Has to be same length as stype
mlist = ['s','^','D','o','X']



# Figure S2
fig,ax = plt.subplots(figsize=(4,4))

# Plot data
xdata = np.log(data['pBnum_2_3'])
ydata = np.log(data['B_2_3'])

# Loop through each site type to put a different marker
for m,s in zip(mlist,stype):
    inds = data['Type']==s
    ax.scatter(xdata[inds],ydata[inds],marker=m,c=np.log(data['S'][inds]),cmap=cm,norm=norm,edgecolor='0.3')

# Set range
ymin = np.floor(np.min(ydata))
ymax = np.ceil(np.max(ydata))
xmin = np.floor(np.min(xdata))
xmax = np.ceil(np.max(xdata))
# Set range min as min of those
rmin = np.min([ymin,xmin])
rmax = np.max([ymax,xmax])
ax.set_ylim(rmin,rmax)
ax.set_xlim(rmin,rmax)

# Labels
ax.set_xlabel('ln(predicted biomass)')
ax.set_ylabel('ln(observed biomass)')

# Add in R^2 value from regression
lin = linregress(xdata,ydata)
xlin = np.linspace(xmin,xmax)
ax.annotate(r'$R^2 = {:.3f}$'.format(lin[2]**2),(0.73,0.17),xycoords='figure fraction')

# Legend
# Plot a bunch of empty points. Not sure if this is the best way, but it's how I'm doing it!
leg = {}
for m,s in zip(mlist,stype):
    leg[s], = ax.plot([],[],c='0.3',marker=m,linestyle="None")

# Plot 1:1 line at the back and add to legend codes
xrange = np.linspace(rmin,rmax)
leg['1:1'], = ax.plot(xrange,xrange,lw=2,c='0.3',zorder=0)#,label='1:1 line') # Can add this back in
lcodes = np.insert(stype,0,'1:1')

ax.legend([leg[s] for s in lcodes],lcodes,prop={"size":7.3})#,frameon=False)#borderpad=0.0)

# Save
fig.savefig('Figures/figS3.pdf',bbox_inches='tight')


# In[16]:


# Figure S3
fig,ax = plt.subplots(figsize=(4,4))

# Plot data
xdata = data['E_2_3']/data['pBnum_2_3']**(2/3)
ydata = data['E_2_3']/data['B_2_3']**(2/3)
# Loop through each site type to put a different marker
for m,s in zip(mlist,stype):
    inds = data['Type']==s
    ax.scatter(xdata[inds],ydata[inds],marker=m,c=np.log(data['S'][inds]),cmap=cm,norm=norm,edgecolor='0.3')

# Set range
ymin = np.floor(np.min(ydata))
ymax = np.ceil(np.max(ydata))
xmin = np.floor(np.min(xdata))
xmax = np.ceil(np.max(xdata))
# Set range min as min of those
rmin = np.min([ymin,xmin])
rmax = np.max([ymax,xmax])
ax.set_ylim(rmin,rmax)
ax.set_xlim(rmin,rmax)

# Labels
ax.set_xlabel(r'Predicted ratio $E:B^{2/3}$')
ax.set_ylabel(r'Observed ratio $E:B^{2/3}$')

# Add in R^2 value from regression
lin = linregress(xdata,ydata)
xlin = np.linspace(xmin,xmax)
ax.annotate(r'$R^2 = {:.3f}$'.format(lin[2]**2),(0.18,0.82),xycoords='figure fraction')

# Legend
# Plot a bunch of empty points. Not sure if this is the best way, but it's how I'm doing it!
leg = {}
for m,s in zip(mlist,stype):
    leg[s], = ax.plot([],[],c='0.3',marker=m,linestyle="None")
    
# Plot 1:1 line at the back and add to legend codes
xrange = np.linspace(rmin,rmax)
leg['1:1'], = ax.plot(xrange,xrange,lw=2,c='0.3',zorder=0)#,label='1:1 line') # Can add this back in
lcodes = np.insert(stype,0,'1:1')

ax.legend([leg[s] for s in lcodes],lcodes,prop={"size":7.3})#,frameon=False#borderpad=0.0)

# Save
fig.savefig('Figures/figS4.pdf',bbox_inches='tight')


# In[ ]:




