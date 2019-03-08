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
md = smf.mixedlm('obsErr~avgR', df, groups=df[['c1','c2','p1','p2']])
mdf = md.fit()
print('++++++++++++++++=', mdf.summary())

