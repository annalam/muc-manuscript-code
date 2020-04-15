#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:41:08 2019

@author: amurtha
"""

import pandas as pd
import numpy as np
import lifelines
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
# from lifelines.utils import datetimes_to_durations
# from lifelines import NelsonAalenFitter
# from lifelines.utils import survival_table_from_events
from lifelines import AalenAdditiveFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
from lifelines.plotting import add_at_risk_counts
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['grid.linewidth'] = 0.5

filename = 'OS_ctDNA_pre1L'

# =============================================================================
# Import data
# =============================================================================

df = pd.read_excel("/groups/wyattgrp/users/amurtha/gillianFigures/data/mUC_tables.xlsx", sheet_name = 'Samples')[['Patient','Sample name','Collection date','Sample type','Cancer %', 'Sample information']]

df_clin = pd.read_excel("/groups/wyattgrp/users/amurtha/gillianFigures/data/mUC_clinical_data_master.xlsx", sheet_name = 'mUC_clean')[['CTDNA_ID','FUPDeath','OS_MetDx', 'FU_Date','DateDx_M']].rename(columns = {'CTDNA_ID':'Patient'})


# =============================================================================
# Merge dataframes, keep only highest ctDNA samples
# =============================================================================

df = df[df['Sample type'] == 'cfDNA']
df = df[(df['Sample information'].str.contains('Untreated'))|(df['Sample information'].str.contains('Pre 1L'))]
df = df.sort_values(['Cancer %','Collection date'], ascending = [False, True])
df = df.drop_duplicates(['Patient'])
df = df.merge(df_clin, how = 'left',on = 'Patient')

# =============================================================================
# Calculate number of months between collection and last follow up
# =============================================================================

df['OS'] = ((df['FU_Date'] - df['Collection date']) / np.timedelta64(1, 'M'))
df['OS'] = df['OS'].astype(float)

# =============================================================================
# Create Categorical data from ctDNA faction
# =============================================================================

df_low = df[df['Cancer %'] < 5]
df_abundant = df[df['Cancer %'] >= 5]

# =============================================================================
# Prepare plots
# =============================================================================

fig1, ax = plt.subplots(1, figsize=(3.75,2.5))

# =============================================================================
# Plot ctDNA low on kmf1
# =============================================================================

kmf1 = KaplanMeierFitter()
color = '#0571b0'
defective_patient_number = str(len(df_low))
defective_label = str(str('Low ctDNA ')+r"(n="+defective_patient_number+")")
T = df_low['OS'].round(3)
C = df_low['FUPDeath'].astype(np.int32)
kmf1.fit(T, event_observed = C, label = defective_label)
kmf1.plot(ax=ax,show_censors = True, ci_show = False, c = color, lw = 1, censor_styles={"ms":8, "clip_on":False}, clip_on=False)
ctDNA_low_median=kmf1.median_


# =============================================================================
# Plot ctDNA abundant on kmf2
# =============================================================================

kmf2 = KaplanMeierFitter()
color = '#ca0020'
defective_patient_number = str(len(df_abundant))
defective_label = str(str('High ctDNA ')+r"(n="+defective_patient_number+")")
T = df_abundant['OS'].round(3)
C = df_abundant['FUPDeath'].astype(np.int32)
kmf2.fit(T, event_observed = C, label = defective_label)
kmf2.plot(ax=ax,show_censors = True, ci_show = False, c = color, lw = 1, censor_styles={"ms":8, "clip_on":False}, clip_on=False)
ctDNA_abundant_median=kmf2.median_


# =============================================================================
# Run Cox progression
# =============================================================================

df_test = df[['OS','FUPDeath','Cancer %']].copy()
df_test.loc[df_test['Cancer %'] < 5 , 'Cancer %'] = 0
df_test.loc[df_test['Cancer %'] >=5, 'Cancer %'] = 1
df_test = df_test.dropna()
cph = CoxPHFitter()
cph.fit(df_test,'OS',event_col = 'FUPDeath', show_progress = True, step_size = 0.001)
coxPHsummary = cph.summary
coxPHsummary['Variable'] = coxPHsummary.index


# =============================================================================
# Adjust plots aethetics
# =============================================================================

ax.plot([0,0],[0,0],color = 'w',alpha=0,label = 'p='+str(round(coxPHsummary.p[0],3))
+"\nHR "+str(round(coxPHsummary['exp(coef)'][0],1))
+" (95% CI "+str(round(coxPHsummary['exp(coef) lower 95%'][0],1))
+" to "+str(round(coxPHsummary['exp(coef) upper 95%'][0],1))+")")
legend = ax.legend(fontsize = 8, loc='best')
#plt.setp(legend.get_title(),fontsize=8)

ax.set_xlim(0,36)
ax.set_xticks(np.arange(0,39,3))

ax.set_ylim(0,1)

ax.set_ylabel('Survival fraction', fontsize=8, labelpad=2)
ax.set_xlabel('OS from cfDNA collection (mo.)', fontsize=8,labelpad=2)
ax.tick_params(axis='y',which="both", length=0, pad=2, left=True, reset=False, labelleft=True)
ax.tick_params(axis='x',which="both", length=0, pad=2, left=True, reset=False, labelleft=True)
ax.grid(b=False)
ax.grid(axis='y', alpha=0.7, linewidth=0.5, linestyle='dotted', zorder=-1, clip_on=False)
ax.grid(axis='x', alpha=0.7, linewidth=0.5, linestyle='dotted', zorder=-1, clip_on=False)


plt.tick_params(axis='both', which='major', labelsize=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

#add_at_risk_counts(kmf1,kmf2,ax=ax, size=6)

plt.tight_layout()














