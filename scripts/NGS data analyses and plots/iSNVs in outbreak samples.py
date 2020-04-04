#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 19:27:59 2020

@author: viro102
"""


## Plots for outbreak samples 

# GI.2

import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
import numpy as np
import os


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/outbreak_samples/GI.2/21_japan/")
rep_a = pd.read_csv("21_japan_GI.2_final.tsv", sep="\t")
rep_a["S_NS"] = np.where(rep_a["REF_AA"] == rep_a["ALT_AA"], "S", "NS")
rep_b=rep_a[rep_a['TOTAL_DP'] >= 400]

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/outbreak_samples/GI.2/21_japan/depth")
japan1_depth = pd.read_csv("21_japan_GI.2.bad_are_masked.sorted.depth", sep = "\t", names = ["Ref", "Pos", "depth"])


sns.lmplot( x="POS", y="ALT_FREQ", data=rep_b, fit_reg=False)



fig, ax1 = plt.subplots(figsize=(15,5)) 
color="tab:green"
#isnv plot creation
ax1.set_ylabel('iSNVs', color="blue", fontsize=20)
ax1.tick_params(axis='y', labelcolor="blue")
ax1=sns.regplot(x="POS", y="ALT_FREQ", data=rep_b, ax=ax1, palette='summer')
#specify we want to share same x axis
ax2 = ax1.twinx()
color = 'tab:red'
# lineplot creation for the depth
sns.lineplot(x='Pos', y='depth', data=japan1_depth, ax=ax2)
ax2.set_ylabel('Depth', color="red")
ax2.tick_params(axis='y', labelcolor="red")
plt.show()



















fig, ax1=plt.subplots()

ax1.set_xlabel('Position')
ax1.set_ylabel('iSNVs', color="blue")
ax1.tick_params(axis='y', labelcolor="blue")
rep_a.plot(x="POS", y="ALT_FREQ", label="1_peru", ax=ax1, 
           kind="scatter", color="blue", alpha=0.6)
plt.xlabel('Position in sapovirus genome', fontsize=15)
plt.ylabel('iSNV frequency')
ax2=ax1.twinx()

ax2.set_ylabel('Depth', color="red")
ax2.tick_params(axis='y', labelcolor="red")
peru1_depth["peru1_depth"].plot(logy=True, ax = ax2, label ="peru1_depth", 
              figsize = (15,5), alpha=0.6, color="red")

plt.title('iSNVs and coverage distribution SP032X-GI.1-peru1', fontsize=20)
plt.grid(axis="both")
plt.savefig('Coverage and iSNVs.in.SP032X-GI.1-peru1.png', dpi=300)
























# GI.2
