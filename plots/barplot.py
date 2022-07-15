import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from decimal import Decimal
from scipy.stats import ttest_ind
from scipy import stats
from itertools import combinations
import glob

stat_funct = stats.ttest_ind

def split_type(df,col,new_col):
    df[col] = df[col].map(lambda x:';'.join(set(x.split(';'))))
    df = (df.drop(col, axis=1)
            .join(df[col].str.split(';', expand=True)
            .stack().reset_index(level=1, drop=True)
            .rename(new_col)))
    return df

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
  



def location(x):
    if x == 'Metabolic':
        x = 0
    elif x == 'Basal':
        x=1
    else:
        x=2
    return x

def vp3(tmp1,tmp2,tmp3,name,es,y1,y2,y3,ym):
    if 1:
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

# =============================================================================
#     y1 = sns.barplot(tmp1,saturation=0.8,estimator=np.mean, capsize=.05,
#                      errwidth=0.6).dataLim.p1[0]
#     plt.close()
#     y2 = sns.barplot(tmp2,saturation=0.8,estimator=np.mean, capsize=.05,
#                      errwidth=0.6).dataLim.p1[0]
#     plt.close()
#     y3 = sns.barplot(tmp3,saturation=0.8,estimator=np.mean, capsize=.05,
#                      errwidth=0.6).dataLim.p1[0]
#     plt.close()
# =============================================================================
# =============================================================================
#     if max(y1,y2,y3) == y1:
#         if max(y2,y3) == y2:
#             yy1 = y1*1.01
#             yy2 = y2*1.01
#             yy3 = y1*1.06
#         else:
#             yy1 = y1*1.01
#             yy2 = y3*1.01
#             yy3 = y1*1.06
#     elif max(y1,y2,y3) == y3:
#         if max(y1,y2) == y2:
#             yy1 = y2*1.01
#             yy2 = y3*1.01
#             yy3 = y3*1.06
#         else:
#             yy1 = y1*1.01
#             yy2 = y3*1.01
#             yy3 = y3*1.06
#     else:
#         yy1 = y2*1.01
#         yy2 = y2*1.01
#         yy3 = y2*1.06
# =============================================================================

    sns.set_theme(style="white")
    bar_mean = es.groupby('Group').mean()['Beta']
    
    y_min = min(bar_mean)
   

    ax = sns.barplot(x='Group',y='Beta',data=es,palette=['#E07F6F','#79A323','#1AB5B7'],
                     saturation=0.8,estimator=np.mean, capsize=.05, errwidth=0.6)
    ax.set_xlabel('')
    ax.set_ylabel('Mean Beta',size=15)
    ax.set_xticklabels(['Metabolic','Basal','Mesenchymal'],size=15)
    plt.yticks(size=12)
    plt.ylim([round(y_min*0.95,2),ym])

    if p1<0.05:
        barplot_annotate_brackets(0.05,0.95,'p = %s'%value(p1),yy1,yy1)
    if p2<0.05:
        barplot_annotate_brackets(1.05,1.95,'p = %s'%value(p2),yy2,yy2)
    if p3<0.05:
        barplot_annotate_brackets(0.05,1.95,'p = %s'%value(p3),yy3,yy3)
   
    plt.savefig('barplot-%s.png' %(name),dpi=300,bbox_inches='tight')
    plt.savefig('barplot-%s.pdf' %(name),bbox_inches='tight')
    plt.close()

plt.rcParams['font.sans-serif'] = ['Arial'] 
files = glob.glob("*.csv")
for file,yy1,yy2,yy3,ym in zip(files,[0.645,0.179,0.648,0.682,0.433,0.287],
                            [0.636,0.181,0.662,0.7,0.451,0.294],[0.645,0.184,0.665,0.71,0.454,0.299],
                            [0.65,0.1865,0.68,0.716,0.457,0.302]):
    loc_data = pd.read_csv(file,index_col=0)
    name = file.split('.csv')[0]
    tmp1 = loc_data[loc_data['Group']=='Metabolic']['Beta']
    tmp2 = loc_data[loc_data['Group']=='Basal']['Beta']
    tmp3 = loc_data[loc_data['Group']=='Mesenchymal']['Beta']
    
    es = pd.concat([loc_data[loc_data['Group']=='Metabolic'],
                    loc_data[loc_data['Group']=='Basal'],
                    loc_data[loc_data['Group']=='Mesenchymal']])
    vp3(tmp1,tmp2,tmp3,name,es,yy1,yy2,yy3,ym)
    


