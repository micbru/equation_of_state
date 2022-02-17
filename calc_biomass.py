#!/usr/bin/env python
# coding: utf-8

# This notebook is intended to show how to get biomasses from state variable data using the python code biomass.py

# In[1]:


import numpy as np
import pandas as pd
import biomass as bm

# Now how to use if we have a csv with state variable information
# Import the data
data = pd.read_csv('data_statevariables.csv')


# In[31]:
# Remove some of the extra columns
d2 = data.drop(columns=['pB0','pB1','pBnum','Lon','Lat','Elev_m','lambda_1','lambda_2','beta',
                   'Notes1','Notes2','Notes3','Notes4'])
cols = d2.columns.to_list()
cols2 = cols[:-1]
cols2.insert(11,cols[-1])

# In[35]:
d2[cols2].to_csv('data_statevariables2.csv',index=False)
# In[26]:


# Now add a column for numerical biomass data, and the approximations.
# The "p" indicates predicted, num is numerical, and 0 and 1 indicate the 0th and first order approximations.
# Note: Iterating through like this is actually bad practice for pandas, but this dataset is small, so it's fine
data['pBnum'] = np.zeros(len(data))
data['pBnum_2_3'] = np.zeros(len(data))
data['pB0'] = np.zeros(len(data))
data['pB1'] = np.zeros(len(data))
# Iterate through each row and append the biomass information
for index, row in data.iterrows():
    row23 = {'S': row['S'], 'N': row['N'], 'E': row['E_2_3']}
    data.loc[index,'pBnum_2_3'] = bm.biomass(row23,power=3/2)
    data.loc[index,'pBnum'] = bm.biomass(row)
    data.loc[index,'pB0'] = bm.biomass_approx(row)
    data.loc[index,'pB1'] = bm.biomass_approx(row,order=1)
    lambdas = bm.mete_lambdas(row)
    data.loc[index,'lambda_1'] = lambdas[0]
    data.loc[index,'lambda_2'] = lambdas[1]
    data.loc[index,'beta'] = lambdas.sum()
    # How are these lambdas doing? 
    #This should be very close to zero if the approximations in the underlying functions hold
    print(bm.constraints(lambdas,row))


# So the approximations are all holding well, though the Farewell Spit dataset has E not that much greater than N, and so the lambda calculations are not as precise (still within one percent). 
# Additionally, for this dataset, the sum of lambda1 and lambda2 (beta) is less than zero, and so our approximation for the biomass equation does not work.

# In[27]:


# To save
data.to_csv('data_biomass.csv',index=False)




