import pandas as pd
import numpy as np
import glob
import os
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines.statistics import logrank_test
from oebio.pipeline.coexpress import corr
import scipy
import re
import shutil
from scipy import stats
from collections import Counter

#加载t-test检验方法
stat_funct = stats.ttest_ind
#stat_funck = stats.kruskal



def ck(dirc):
    if not os.path.exists(dirc):
        os.mkdir(dirc)

def value(q):
    if q >=0.001:
        return round(q,3)
    else:
        return '%.2e' %q


def barplot_annotate_brackets(num1, num2, text,  y1,y2, yerr=None, dh=.001, barh=.01, fs=None, maxasterix=None):

    lx, ly = num1, y1
    rx, ry = num2, y2
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh = 0.01*y1

    y = max(ly, ry) + dh
    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, 1.01*y1)

    plt.plot(barx, bary, c='black',linewidth=1)
    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs
    plt.text(*mid, text, **kwargs,size=14)



def vp3(tmp1,tmp2,tmp3,name,es,typ,ylab,s):
    if typ=='T-test':
        try:
            p1 = stat_funct(tmp1, tmp2).pvalue
        except:
            p1 = 'NAN'
        try:
            p2 = stat_funct(tmp2, tmp3).pvalue
        except:
            p2 = 'NAN'
        try:
            p3 = stat_funct(tmp1, tmp3).pvalue
        except:
            p3 = 'NAN'

    y1 = sns.boxplot(tmp1).dataLim.p1[0]
    plt.close()
    y2 = sns.boxplot(tmp2).dataLim.p1[0]
    plt.close()
    y3 = sns.boxplot(tmp3).dataLim.p1[0]
    plt.close()
    if max(y1,y2,y3) == y1:
        if max(y2,y3) == y2:
            yy1 = y1*1.01
            yy2 = y2*1.01
            yy3 = y1*1.06
        else:
            yy1 = y1*1.01
            yy2 = y3*1.01
            yy3 = y1*1.06
    elif max(y1,y2,y3) == y3:
        if max(y1,y2) == y2:
            yy1 = y2*1.01
            yy2 = y3*1.01
            yy3 = y3*1.06
        else:
            yy1 = y1*1.01
            yy2 = y3*1.01
            yy3 = y3*1.06
    else:
        yy1 = y2*1.01
        yy2 = y2*1.01
        yy3 = y2*1.06


    ax = sns.boxplot(x='Group',y=name,data=es,palette=['#E07F6F','#79A323','#1AB5B7'],showfliers=False)
    ax.set_xlabel('')
    ax.set_ylabel(ylab,size=15)
    ax.set_xticklabels(['Metabolic','Basal','Mesenchymal'],size=15)
    ax.set_title('MAPK signaling pathway',size=18)
    m1,m2 = ax.get_ylim()
    m1,m2 = ax.set_ylim(m1,m2*1.1)
    plt.close()
    protein_Metabolic = tmp1.to_list()
    protein_Basal = tmp2.to_list()
    protein_Mesenchymal = tmp3.to_list()
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    box = plt.boxplot([protein_Metabolic, protein_Basal, protein_Mesenchymal], patch_artist=True, widths=0.6, showcaps=True,
                      medianprops={"color": "black", 'lw': 4}, boxprops={'lw': 2}, capprops={'lw': 2},
                      positions=[0, 1, 2], whiskerprops={'lw': 2, 'ls': '--'}, showfliers=False, labels=["Metabolic","Basal","Mesenchymal"])

    sns.stripplot(data=[protein_Metabolic, protein_Basal, protein_Mesenchymal], orient="v",
                  palette=sns.color_palette(['#E07F6F', '#79A323', '#1AB5B7'], n_colors=3), jitter=0.08, size=6, ax=ax)

    # box =sns.violinplot(data =[protein_Metabolic, protein_Basal, protein_Mesenchymal],palette=['#E07F6F', '#79A323', '#1AB5B7'],inner=None,zorder=1)
    plt.title('MAPK signaling pathway',size=20)
    plt.xticks([0,1,2],["Metabolic","Basal","Mesenchymal"])
    ax.set_ylabel(ylab,size=20)
    plt.xticks(size=20)
    plt.yticks(size=12)

    for box_, median_,color in zip(box['boxes'], box["medians"],['#E07F6F','#79A323','#1AB5B7']):
        # box_.set_color(color)
        box_.set_facecolor(color)
        box_.set_lw(2)
        box_.set_alpha(0.85)
    #y = max(ax.get_ylim())
    #ax.text(min(ax.get_xlim())*0.8,y*.92,can,size=15)
    barplot_annotate_brackets(0.05,0.95,'p = %s'%value(p1),yy1,yy1)
    #barplot_annotate_brackets(1.05,1.95,'p = %s'%value(p2),yy2,yy2)
    barplot_annotate_brackets(0.05,1.95,'p = %s'%value(p3),yy3,yy3)
    plt.savefig('%s_%s_boxplot.png' %(name,s),dpi=300,bbox_inches='tight')
    plt.savefig('%s_%s_boxplot.pdf' %(name,s),bbox_inches='tight')
    plt.close()




###part6
mm = pd.read_csv('new_order.txt',sep='\t')
mm = mm.replace('cluster1-Kinase','Metabolic')
mm = mm.replace('cluster2-Immune','Basal')
mm = mm.replace('cluster3-ECM','Mesenchymal')
dic = mm.set_index('comb')['Group'].to_dict()

es = pd.read_excel('KEGG_all_samples_ES_order_data.xlsx',index_col=0).loc[['KEGG_MAPK_SIGNALING_PATHWAY']]
es.columns = [x.split('_')[-1] for x in es.columns]
es = es.T.reset_index()
es['Group'] = es['index'].apply(lambda x:dic[x])
es['Group'] = es['Group'].astype('category')
es['Group'].cat.set_categories(['Metabolic','Basal','Mesenchymal'],inplace=True)
es = es.sort_values(by='Group')
name='KEGG_MAPK_SIGNALING_PATHWAY'

tmp1 = es[es['Group']=='Metabolic'][name]
tmp2 = es[es['Group']=='Basal'][name]
tmp3 = es[es['Group']=='Mesenchymal'][name]

vp3(tmp1,tmp2,tmp3,name,es,'T-test','ES values','ES')

nes = pd.read_excel('KEGG_all_samples_NES_order_data.xlsx',index_col=0).loc[['KEGG_MAPK_SIGNALING_PATHWAY']]
nes.columns = [x.split('_')[-1] for x in nes.columns]
nes = nes.T.reset_index()
nes['Group'] = nes['index'].apply(lambda x:dic[x])
nes['Group'] = nes['Group'].astype('category')
nes['Group'].cat.set_categories(['Metabolic','Basal','Mesenchymal'],inplace=True)
nes = nes.sort_values(by='Group')
name='KEGG_MAPK_SIGNALING_PATHWAY'

tmp1 = nes[nes['Group']=='Metabolic'][name]
tmp2 = nes[nes['Group']=='Basal'][name]
tmp3 = nes[nes['Group']=='Mesenchymal'][name]

vp3(tmp1,tmp2,tmp3,name,nes,'T-test','NES values','NES')






