#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:30:26 2020
iSNV per sample peru
@author: viro102
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os as os 


# PERU

# SP032X
os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP032X/2_peru/")
rep_b = pd.read_csv("2_peru_final.tsv", sep="\t")
rep_b= rep_b[rep_b['TOTAL_DP'] >= 400]
g = sns.lmplot( x="POS", y="ALT_FREQ", data=rep_b, fit_reg=False, 
               hue='S_or_NS',aspect=20/8, legend=False, 
               markers=["o", "x"], scatter_kws={"s": 90})
# legend title
new_title = 'Amino acid change'
g._legend.set_title(new_title)
# replace legend labels
new_labels = ['Synonymous', 'Non-synonymous']
for t, l in zip(g._legend.texts, new_labels): t.set_text(l)

# resize figure box to -> put the legend out of the figure
box = g.ax.get_position() # get position of figure
g.ax.set_position([box.x0, box.y0, box.width * 0.85, box.height]) # resize position

# Put a legend to the right side
g.ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
sns.plt.show(g)
plt.xlabel('Genome position', fontsize=15)
plt.ylabel('iSNV frequency', fontsize=15)
plt.title('iSNVs in sample 2_peru', fontsize=20)
plt.grid(axis="both")
plt.show()






# SP032X
os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP032X/2_peru/")
rep_b = pd.read_csv("2_peru_final.tsv", sep="\t")
rep_b= rep_b[rep_b['TOTAL_DP'] >= 400]
sns.lmplot( x="POS", y="ALT_FREQ", data=rep_b, hue='S_or_NS', fit_reg=False)
plt.xlabel('Genome position', fontsize=15)
plt.ylabel('iSNV frequency', fontsize=15)
plt.title('iSNVs in sample 2_peru', fontsize=20)
plt.show()







os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/outbreak_samples/GI.1/10_japan/")
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
rep_a = pd.read_csv("10_japan_GI.1_final.tsv", sep="\t")
rep_a["S_NS"] = np.where(rep_a["REF_AA"] == rep_a["ALT_AA"], "S", "NS")
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
rep_b=rep_a[rep_a['TOTAL_DP'] >= 400]
sns.lmplot( x="POS", y="ALT_FREQ", data=rep_b, fit_reg=False, hue='S_NS',aspect=14/5, 
           legend=True, legend_out=True)
plt.xlabel('Genome position', fontsize=15)
plt.ylabel('iSNV frequency', fontsize=15)
plt.title('iSNVs in sample 10_miyagi', fontsize=20)
plt.grid(axis="both")
plt.show()

