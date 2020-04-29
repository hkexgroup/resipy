#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 19:15:50 2019

@author: jkl
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from api.R2 import R2
k = R2()
k.createSurvey('./api/test/syscalFile.csv')
k.lmefit()

s = k.surveys[0]
dfg = s.df[s.df['irecip'] > 0]
recipMean = np.abs(dfg['recipMean'].values)
recipError = np.abs(dfg['recipError'].values)
array = dfg[['a','b','m','n']].values.astype(int)
data = np.vstack([recipMean, recipError]).T
data = np.hstack((data, array))
df = pd.DataFrame(data, columns=['avgR','obsErr','c1','c2','p1','p2'])
df.to_csv('/home/jkl/Downloads/df.csv')
groups = df[['c1','c2','p1','p2']].astype(int)
#print(list(set(groups)))
#x = [str(a[0]) + str(a[1]) + str(a[2]) + str(a[3]) for a in df[['c1','c2','p1','p2']].values.astype(int)]
#groups = pd.Series(x, name='g') # this works but not same number of groups ...
#groups = ['c1','c2','p1','p2']
#groups = [df['c1'].values, df['c2'].values, df['p1'].values, df['p2'].values]
md = smf.mixedlm('obsErr~avgR', df, groups=groups)
mdf = md.fit()
print('++++++++++++++++=', mdf.summary())


"""
    A mixed model with fixed effects for the columns of ``exog`` and
    independent random coefficients for the columns of ``exog_re``:

    >>> free = MixedLMParams.from_components(fe_params=np.ones(exog.shape[1]), \
                     cov_re=np.eye(exog_re.shape[1]))
    >>> model = sm.MixedLM(endog, exog, groups, exog_re=exog_re)
    >>> result = model.fit(free=free)

    A different way to specify independent random coefficients for the
    columns of ``exog_re``.  In this example ``groups`` must be a
    Pandas Series with compatible indexing with ``exog_re``, and
    ``exog_re`` has two columns.

    >>> g = pd.groupby(groups, by=groups).groups
    >>> vc = {}
    >>> vc['1'] = {k : exog_re.loc[g[k], 0] for k in g}
    >>> vc['2'] = {k : exog_re.loc[g[k], 1] for k in g}
    >>> model = sm.MixedLM(endog, exog, groups, vcomp=vc)
    >>> result = model.fit()
"""


# SOLUTION: in 0.9.0, commenting the following line in the constructor of the
# linear mixed effect did it (around line 936)
#        else:
#            groups = np.asarray(groups)
