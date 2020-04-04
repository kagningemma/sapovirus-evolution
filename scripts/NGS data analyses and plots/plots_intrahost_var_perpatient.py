# -*- coding: utf-8 -*-
"""
Spyder Editor
This is script for plotting of iSNV of each individual on each plot to see if there
are accumulation of SNPs
"""

import matplotlib.pyplot as plt
import pandas as pd 
import os

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP032X/1_peru/cons2")
rep_a = pd.read_csv("1_peru_final.tsv", sep="\t")
rep_a=rep_a[rep_a['TOTAL_DP'] >= 400]
rep_a

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP032X/2_peru/")
rep_b = pd.read_csv("2_peru_final.tsv", sep="\t")
rep_b= rep_b[rep_b['TOTAL_DP'] >= 400]
rep_b


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP032X/5_peru/")
rep_c= pd.read_csv("5_peru_final.tsv", sep="\t")
rep_c= rep_c[rep_c['TOTAL_DP'] >= 400]
rep_c


f, ax = plt.subplots(figsize=(15,5))
rep_a.plot(x="POS", y="ALT_FREQ", label="peru_1-Day1-n3", ax=ax, s=50,
           kind="scatter", color="green", alpha=0.5)
rep_b.plot(x="POS", y="ALT_FREQ", label="peru_2-Day4-n5", ax=ax, s=50,
           kind="scatter", color="blue", alpha=0.5)
rep_c.plot(x="POS", y="ALT_FREQ", label="peru_5-Day27-n9", ax=ax, s=50,
           kind="scatter", color="red", alpha=0.5)

plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.title('sapovirus GI.1 iSNVs in child-SP031X at different time points', fontsize=25)

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))

plt.savefig('GI.1_iSNV_plots_SP031X.png', dpi=300)

# next person

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP051X/9_peru/cons2")
rep_d = pd.read_csv("9_peru_final.tsv", sep="\t")
rep_d=rep_d[rep_d['TOTAL_DP'] >= 400]
rep_d

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP051X/11_peru/")
rep_e = pd.read_csv("11_peru_final.tsv", sep="\t")
rep_e= rep_e[rep_e['TOTAL_DP'] >= 400]
rep_e

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP051X/12_peru/")
rep_f = pd.read_csv("12_peru_final.tsv", sep="\t")
rep_f= rep_f[rep_f['TOTAL_DP'] >= 400]
rep_f



f, ax = plt.subplots(figsize=(15,5))
rep_d.plot(x="POS", y="ALT_FREQ", label="peru_9-Day1-n11", ax=ax, s=50,
           kind="scatter", color="green", alpha=0.5)
rep_e.plot(x="POS", y="ALT_FREQ", label="peru_11-Day8-n8", ax=ax, s=50,
           kind="scatter", color="blue", alpha=0.5)
rep_f.plot(x="POS", y="ALT_FREQ", label="peru_12-Day16-n25", ax=ax, s=50,
           kind="scatter", color="red", alpha=0.5)

plt.title('sapovirus GI.1 SNPs in child-SP051X at different time points', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
ax.legend(loc='center left', bbox_to_anchor=(1,1))
plt.savefig('GI.1_iSNV_plots_SP051X.png', dpi=300)

# next person


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/19_peru/cons2")
rep_g = pd.read_csv("19_peru_final.tsv", sep="\t")
rep_g=rep_g[rep_g['TOTAL_DP'] >= 400]
rep_g

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/22_peru/")
rep_h = pd.read_csv("22_peru_final.tsv", sep="\t")
rep_h= rep_h[rep_h['TOTAL_DP'] >= 400]
rep_h

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/25_peru/")
rep_i = pd.read_csv("25_peru_final.tsv", sep="\t")
rep_i= rep_i[rep_i['TOTAL_DP'] >= 400]
rep_i

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/26_peru/")
rep_j = pd.read_csv("26_peru_final.tsv", sep="\t")
rep_j= rep_j[rep_j['TOTAL_DP'] >= 400]
rep_j

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/28_peru/")
rep_k = pd.read_csv("28_peru_final.tsv", sep="\t")
rep_k= rep_k[rep_k['TOTAL_DP'] >= 400]
rep_k

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP088X/30_peru/")
rep_l = pd.read_csv("30_peru_final.tsv", sep="\t")
rep_l= rep_l[rep_l['TOTAL_DP'] >= 400]
rep_l



f, ax = plt.subplots(figsize=(15,5))
rep_g.plot(x="POS", y="ALT_FREQ", label="peru_19-Day1-n5", ax=ax, s=50,
           kind="scatter", color="green", alpha=0.5)
rep_h.plot(x="POS", y="ALT_FREQ", label="peru_22-Day9-n6", ax=ax, s=50,
           kind="scatter", color="blue", alpha=0.5)
rep_i.plot(x="POS", y="ALT_FREQ", label="peru_25-Day15-n8", ax=ax, s=50,
           kind="scatter", color="black", alpha=0.5)
rep_j.plot(x="POS", y="ALT_FREQ", label="peru_26-Day21-n10", ax=ax,s=50, 
           kind="scatter", color="#e34a33", alpha=0.5)
# Put a legend to the right of the current axis

plt.title('sapovirus GI.1 SNPs in child-SP088X at different time points', fontsize=25)
plt.xlabel('Position', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)


ax.legend(loc='center left', bbox_to_anchor=(1,1))
plt.savefig('GI.1_iSNV_plots_SP088X.png', dpi=300)

plt.show()


# GI.2

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP0121X/31_peru/cons2")
rep_m = pd.read_csv("31_peru_final.tsv", sep="\t")
rep_m= rep_k[rep_k['TOTAL_DP'] >= 400]
rep_m

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP0121X/33_peru/")
rep_n = pd.read_csv("33_peru_final.tsv", sep="\t")
rep_n= rep_l[rep_l['TOTAL_DP'] >= 400]
rep_n


f, ax = plt.subplots(figsize=(15,5))
rep_m.plot(x="POS", y="ALT_FREQ", label="peru_31-Day1-n5", ax=ax, s=50,
           kind="scatter", color="green", alpha=0.5)
rep_n.plot(x="POS", y="ALT_FREQ", label="peru_33-Day9-n6", ax=ax, s=50,
           kind="scatter", color="red", alpha=0.5)

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/iSNVs_plots_per_person/iSNV plots")
plt.title('sapovirus GI.2_iSNV at different time points in child-SP121X', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(fontsize=20)


plt.savefig('GI.2_iSNV in SP121X.png', dpi=300)



# NEXT person


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP223X/14_peru/cons2")
rep_o = pd.read_csv("14_peru_final.tsv", sep="\t")
rep_o=rep_o[rep_o['TOTAL_DP'] >= 400]
rep_o

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP223X/15_peru/")
rep_p = pd.read_csv("15_peru_final.tsv", sep="\t")
rep_p= rep_p[rep_p['TOTAL_DP'] >= 400]
rep_p

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP223X/16_peru/")
rep_q = pd.read_csv("16_peru_final.tsv", sep="\t")
rep_q= rep_q[rep_q['TOTAL_DP'] >= 400]
rep_q

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP223X/17_peru/")
rep_r = pd.read_csv("17_peru_final.tsv", sep="\t")
rep_r= rep_r[rep_r['TOTAL_DP'] >= 400]
rep_r

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP223X/18_peru/")
rep_s = pd.read_csv("18_peru_final.tsv", sep="\t")
rep_s= rep_s[rep_s['TOTAL_DP'] >= 400]
rep_s


f, ax = plt.subplots(figsize=(15,5))
rep_o.plot(x="POS", y="ALT_FREQ", label="peru_14-Day1-n8", ax=ax, s=50, 
           kind="scatter", color="green", alpha=0.5)
rep_p.plot(x="POS", y="ALT_FREQ", label="peru_15-Day5-n12", ax=ax, s=50,
           kind="scatter", color="blue", alpha=0.5)
rep_q.plot(x="POS", y="ALT_FREQ", label="peru_16-Day12-n13", ax=ax, s=50,
           kind="scatter", color="grey", alpha=0.5)
rep_r.plot(x="POS", y="ALT_FREQ", label="peru_17-Day19-n18", ax=ax, s=50,
           kind="scatter", color="magenta", alpha=0.5)
rep_s.plot(x="POS", y="ALT_FREQ", label="peru_18-Day26-n30", ax=ax, s=50,
           kind="scatter", color="red", alpha=0.5)

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))


plt.title('sapovirus GI.2 iSNPs at different time points in SP223X', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(axis="both", fontsize=20)

plt.savefig('GI.2_iSNV_plots_SP223X.png', dpi=300)






# NEXT person


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP265X/20_peru/cons2")
rep_t = pd.read_csv("20_peru_final.tsv", sep="\t")
rep_t=rep_t[rep_t['TOTAL_DP'] >= 400]
rep_t

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP265X/21_peru/")
rep_u = pd.read_csv("21_peru_final.tsv", sep="\t")
rep_u= rep_u[rep_u['TOTAL_DP'] >= 400]
rep_u

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP265X/23_peru/")
rep_v = pd.read_csv("23_peru_final.tsv", sep="\t")
rep_v= rep_v[rep_v['TOTAL_DP'] >= 400]
rep_v

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SP265X/27_peru/")
rep_w = pd.read_csv("27_peru_final.tsv", sep="\t")
rep_w= rep_w[rep_w['TOTAL_DP'] >= 400]
rep_w



f, ax = plt.subplots(figsize=(15,5))
rep_t.plot(x="POS", y="ALT_FREQ", label="peru_20-Day1-n10", ax=ax, s=50,
           kind="scatter", color="green", alpha=0.5)
rep_u.plot(x="POS", y="ALT_FREQ", label="peru_21-Day5-n12", ax=ax, s=50,
           kind="scatter", color="blue", alpha=0.5)
rep_v.plot(x="POS", y="ALT_FREQ", label="peru_23-Day12-n25", ax=ax, s=50,
           kind="scatter", color="grey", alpha=0.5)
rep_w.plot(x="POS", y="ALT_FREQ", label="peru_27-Day19-n31", ax=ax, s=50,
           kind="scatter", color="magenta", alpha=0.5)


# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))

ax.legend(loc='center left', bbox_to_anchor=(1,1))
plt.title('sapovirus GI.2 iSNPs at different time points in SP265X', fontsize=25)
plt.xlabel('genome position', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(axis="both", fontsize=20)


plt.savefig('GI.2_iSNV_plots_SP265X.png', dpi=300)




# Make iSNV Plots for all GI.1 samples on same plots 



f, ax = plt.subplots(figsize=(15,5))
rep_a.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_b.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_c.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_d.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5),
rep_e.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_f.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_g.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_h.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_i.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_j.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)

rep_m.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_n.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
            kind="scatter", alpha=0.5)
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))
plt.title('Sapovirus GI.1 iSNPs in all samples from Peru', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(axis="both", fontsize=20)

plt.savefig('Intra-host.SNPs.in.all GI.1 samples from peru.png', dpi=300)

# Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))







# Make iSNV Plots for all GI.2 samples on same plots 

f, ax = plt.subplots(figsize=(15,5))
rep_o.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_p.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_q.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_r.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_s.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_t.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_u.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_v.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_w.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)


# Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1,1))

plt.title('Sapovirus GI.2 iSNPs in all samples from Peru', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(axis="both", fontsize=20)
plt.savefig('Intra-host.SNPs.in.all GI.2 samples from peru.png', dpi=300)








# Make iSNV Plots for all samples on same plots 



f, ax = plt.subplots(figsize=(15,5))
rep_a.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_b.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_c.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_d.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)

rep_e.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_f.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_g.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_h.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_i.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_j.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_m.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_n.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_o.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_p.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_q.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_r.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_s.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_t.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_u.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_v.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)
rep_w.plot(x="POS", y="ALT_FREQ", legend=None, ax=ax, s=50,
           kind="scatter", alpha=0.5)


# Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.legend(loc='center left', bbox_to_anchor=(1,1))

plt.title('Sapovirus GI.1 and GI.2 iSNPs in all samples from Peru', fontsize=25)
plt.xlabel('Genome positions', fontsize=20)
plt.ylabel('iSNV frequency', fontsize=20)
plt.grid(axis="both", fontsize=20)
plt.savefig('Intra-host.SNPs.in.all.samples.png', dpi=300)





































#Philipino samples from 4th batch



# NEXT person

import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SN091X_batch4/1997/")
rep_x = pd.read_csv("1997_phil_batch3_final.tsv", sep="\t")
rep_x=rep_x[rep_x['TOTAL_DP'] >= 400]


os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SN091X_batch4/2010/")
rep_y = pd.read_csv("2010_phil_batch3_final.tsv", sep="\t")
rep_y= rep_y[rep_y['TOTAL_DP'] >= 400]



f, ax = plt.subplots(figsize=(15,5))
rep_x.plot(x="POS", y="ALT_FREQ", label="1997_phil_batch4-day1-198", ax=ax, 
           kind="scatter", color="green", alpha=0.8)
rep_y.plot(x="POS", y="ALT_FREQ", label="2010_phil_batch4-Day10-n169", ax=ax,  
           kind="scatter", alpha=0.8)

# Set grid to use minor tick locations. 
ax.grid(which = 'minor')
plt.title('Intra-host SNPs in GI.2.SN091X')
plt.xlabel('Position')
plt.ylabel('iSNV frequency')
plt.grid()

plt.savefig('GI.2_iSNV_plots_SN091X.png', dpi=300)





# NEXT person

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SN032X_batch4/285/")
rep_z1 = pd.read_csv("285_phil_batch3_final.tsv", sep="\t")
rep_z1=rep_x[rep_x['TOTAL_DP'] >= 400]
rep_z1

os.chdir("/media/viro102/HD-ADU3/kte-data/sapo-ivar/test/SN032X_batch4/372/")
rep_z2 = pd.read_csv("372_phil_batch3_final.tsv", sep="\t")
rep_z2= rep_z2[rep_z2['TOTAL_DP'] >= 400]
rep_z2



f, ax = plt.subplots(figsize=(15,5))
rep_z1.plot(x="POS", y="ALT_FREQ", label="285_phil_batch4-day6-154", ax=ax, 
           kind="scatter", color="green", alpha=0.8)
rep_z2.plot(x="POS", y="ALT_FREQ", label="2010_phil_batch4-Day60-n118", ax=ax,  
           kind="scatter", alpha=0.8)

# Set grid to use minor tick locations. 
ax.grid(which = 'minor')
plt.title('Intra-host SNPs in GI.2.SN032X')
plt.xlabel('Position')
plt.ylabel('iSNV frequency')
plt.grid()

plt.savefig('GI.1_iSNV_plots_SN032X.png', dpi=300)






### OVERLAY OF iSNV and Coverage plot



