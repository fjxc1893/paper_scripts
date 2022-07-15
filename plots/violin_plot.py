import pandas as pd
import os
from pandas.api.types import CategoricalDtype
from scipy.stats import ttest_ind
from PIL import Image
import matplotlib.patches as patches

dic = {}
def types(x):
    if x == 'cluster1-Kinase':
        return 'Metabolic'
    elif x == 'cluster2-Immune':
        return 'Basal'
    elif x == 'cluster3-ECM':
        return 'Mesenchymal'
    return 'issue'

def types1(x):
    if x == '1':
        return 'Metabolic'
    elif x == '2':
        return 'Basal'
    elif x == '3':
        return 'Mesenchymal'
    return 'issue'


def groups(x):
    return df[df['Sample'] == x]['Subtype'].to_string().split(' ')[-1]

def get_t_test(pro,prot):
    pro = pro.to_frame()
    #print(prot.columns)
    pro['group'] = prot.columns
    pro = pro.iloc[6:-3, :]
    pro['group'] = pro['group'].apply(groups)
    #print(pro)
    category = CategoricalDtype(['Metabolic', 'Basal', 'Mesenchymal'], ordered=True)
    pro['group'] = pro['group'].astype(category)
    pro = pro.sort_values(by=['group'])
    #print(pro)
    pro.columns = ['num', 'group']
    pro['num'] = pro['num'].astype('float64')
    Metabolic = pro[pro['group'] == 'Metabolic']['num']
    Basal = pro[pro['group'] == 'Basal']['num']
    Mesenchymal = pro[pro['group'] == 'Mesenchymal']['num']
    t1 = ttest_ind(Metabolic, Mesenchymal)
    t2 = ttest_ind(Basal, Mesenchymal)
    return t1[1],t2[1]

#plotting('LST')
os.chdir("F:/个性化/HT2020-19711/20211212_3plots")

data = pd.read_csv("cluster.txt", sep = '\t')
prot = pd.read_csv("post_pro_q.txt", sep = '\t')

df = pd.DataFrame(data)
df = df.sort_values('Proteomic Subtype')
df['Subtype'] = df['Proteomic Subtype'].apply(types)
df['Telomeric AI'] = df['Telomeric_AI']
#test = prot.iloc[1:10,:]
#prot['t_test1'],prot['t_test2'] = zip(*prot.apply(get_t_test,axis = 1, args= (prot,)))

genes = ['EGFR','PDGFRB','FLT1','AXL','PDGFC','GAS6','HDAC1','HDAC2']
#genes = ['PDGFC']
import plotting
for gene in genes:
    #FLT1 (VEGFR1)
    pro = prot[prot['Gene Name'] == gene]

    t_test1 = pro.iloc[0,-4]
    t_test2 = pro.iloc[0,-3]
    #print(t_test1 + ' ' + t_test2)
    pro = pro.unstack().to_frame()
    pro['group'] = prot.columns

    pro = pro.iloc[6:-7, :]
    print(pro)
    pro['group'] = pro['group'].apply(groups)
    category = CategoricalDtype(['Metabolic','Basal',  'Mesenchymal'], ordered=True)
    pro['group'] = pro['group'].astype(category)
    pro = pro.sort_values(by = ['group'])
    print(pro)
    pro[0] = pro[0].astype('float64')
    print(pro[0])
    if gene == 'AXL' :
        plotting.violin_plotting(pro, gene, top_p = 2.2,  line_height_p=2.8, t1 = t_test1, t2 = t_test2)
    elif gene == 'PDGFC':
        plotting.violin_plotting(pro, gene, top_p=1.2, t1=t_test1, t2=t_test2)
    elif gene == 'FLT1':
        plotting.violin_plotting(pro,'FLT1 (VEGFR1)', text_fs = 30, t1 = t_test1, t2 = t_test2)
    elif gene == 'VEGFA':
        plotting.violin_plotting(pro, gene, top_p = 2, t1 = t_test1, t2 = t_test2)
    elif gene == 'HDAC2':
        plotting.violin_plotting(pro, gene, top_p = 2.2,  line_height_p=2.8, t1 = t_test1, t2 = t_test2)
    elif gene == 'HDAC1':
        plotting.violin_plotting(pro, gene, top_p=1.7, line_height_p=2.8, t1 = t_test1, t2 = t_test2)
    else:
        plotting.violin_plotting(pro,gene, top_p = 2.2,  t1 = t_test1, t2 = t_test2)