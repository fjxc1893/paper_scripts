import pandas as pd
import numpy as np
from scipy import stats
import os
import sys
import re
import seaborn as sns
import matplotlib.pyplot as plt
import glob

#加载t-test检验方法
stat_funct = stats.ttest_ind

#加载kruskal检验
#stat_funck = stats.kruskal

#定义差异计算及可视化函数
def deg(df,case,con,cs,cn):
    #命名比较组
    di = cs+'-vs-'+cn
    print(di)
    #获取两个组的数据
    data = df[case+con]
    #获取两个组的样本以及组别信息
    group =  pd.DataFrame({'Sample':case+con,'Group':[cs]*len(case)+[cn]*len(con)})
    #获取样本与组别的字典
    dic = group.set_index('Sample')['Group'].to_dict()
    #对数据做过滤，过滤掉全是0值的基因
    data = data[data.sum(1)>0]
    #获取其中一组的数据
    d1 = data[case]
    #获取另外一组的数据
    d2 = data[con]
    #建一个空字典
    dic_p = {}
    #合并两组的数据
    d = pd.merge(d1,d2,left_index=True,right_index=True)

    dt = d.copy()
    #计算两组数据的T-test检验p值
    pt = stat_funct(d1, d2, axis=1).pvalue
    dt['p_value'] = pt
    #合并数据和组别信息方便查看样本来源
    dt = pd.concat([dt,group.set_index('Sample').T])
    
    #检查是否存在比较组的文件夹，没有则创建
    if not os.path.exists(di):
        os.mkdir(di)
    os.chdir(di)
    
    #对两种差异计算结果做循环处理
    for d,typ in zip([dt],['T-test']):
        d.index.name = 'Gene'
        
        #输出没有的结果
        d.to_csv('%s.xls' %(typ),sep='\t')
        #对数据做p值筛选
        diff = d[(d['p_value']<=0.05)|(d['p_value'].astype(str)=='nan')]
        pdic = diff['p_value'].to_dict()
        #保存p值小于0.05的结果
        if diff.shape[0]>0:
            diff.to_csv('%s_pvalue-0.05.xls' %(typ),sep='\t')

    os.chdir('../')

#读入数据

mm = pd.read_csv('order-RNAseq.txt',sep='\t')

dic = mm.set_index('comb')['Group'].to_dict()
lst = list(pd.read_csv('list.txt',header=None,sep='\t',index_col=0).index)
stat = pd.read_csv('diff.txt',sep='\t')



os.chdir('xCell')
df0 = pd.read_csv('xCell_101RNA-seq.txt',sep='\t',index_col=0)
df0.columns = [x.replace('.','-') for x in df0.columns]
df0 = df0[df0.index.isin(lst)]
#循环处理样本组别信息
for indes,r in stat.iterrows():
    #获取第一列的样本名称
    cassmp = r['Case'].split('(')[0].split(',')
    #获取第一列的组别信息
    casg = re.findall(r'.*\((.*)\).*',r['Case'])[0]
    #获取第二列的样本名称
    conmp = r['Control'].split('(')[0].split(',')
    #获取第二列的组别名称
    cong = re.findall(r'.*\((.*)\).*',r['Control'])[0]
    #使用自定义的函数计算差异并绘制箱线图    
    deg(df0,cassmp,conmp,casg,cong)

al = glob.glob('*/T-test_pvalue-0.05.xls')
tmp = pd.DataFrame()
for f in al:
    d = pd.read_csv(f,sep='\t')[['Gene','p_value']].dropna()
    tmp = pd.concat([tmp,d])
diff0 = tmp.sort_values(by='p_value').drop_duplicates(subset='Gene',keep='first')
pdic = diff0.set_index('Gene')['p_value'].to_dict()
diff = df0.loc[diff0['Gene']]

def barplot_annotate_brackets(num1, num2, text,  y1,y2, yerr=None, dh=.001, barh=.01, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    lx, ly = num1, y1
    rx, ry = num2, y2


    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y-dh)

    plt.plot(barx, bary, c='black',linewidth=0.5)

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)

def barplot_annotate_brackets1(num1, num2, text,  y1,y2,y3, yerr=None, dh=.001, barh=.01, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    lx, ly = num1, y1
    rx, ry = num1+num2, y2


    #ax_y0, ax_y1 = plt.gca().get_ylim()
    #dh *= (ax_y1 - ax_y0)
    #barh *= (ax_y1 - ax_y0)

    x = num1+num2*2

    barx = [lx, rx, rx, lx]
    bary = [y1, y1, y2, y2]
    mid = (x, y3)

    plt.plot(barx, bary, c='black',linewidth=0.5)

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)


focs0 = df0.iloc[:-1,:].T[['ImmuneScore', 'Macrophages M1', 'Macrophages M2']].T

focs = df0.iloc[:-1,:].T[['iDC', 'cDC', 'aDC']].T

tmp = pd.DataFrame()
for c in focs0.columns:
    d = focs0[[c]]
    d.columns = ['percent']
    d['Group'] = dic[c]
    tmp = pd.concat([tmp,d])
    
tmp.reset_index(inplace=True)



    
plt.rcParams['font.sans-serif'] = ['Arial']
plt.figure(figsize=(6,5.1))
ax = sns.barplot(x='index',y='percent',hue='Group',data=tmp,palette=['#E07F6F','#79A323','#1AB5B7'], estimator=np.mean,ci=15, capsize=.05,errwidth=0.8)
#ax.set_xticklabels(labels=['ImmuneScore', 'Macrophages M1', 'Macrophages M2'],rotation=30)
ax.set_xlabel('')
ax.set_ylabel('Immune Cell Proportion')
ax.legend(frameon=False)
ax.set_ylim(0,0.04)
barplot_annotate_brackets(0.02,0.25,'*',0.0034756,0.034756)
barplot_annotate_brackets(1,1.24,'*',0.019762,0.019762)
barplot_annotate_brackets(2,2.24,'*',0.024987,0.024987)


plt.savefig('barplot1-ver5.1.pdf',bbox_inches='tight')
plt.savefig('barplot1-ver.png',dpi=300,bbox_inches='tight')
plt.close()


tmp = pd.DataFrame()
for c in focs.columns:
    d = focs[[c]]
    d.columns = ['percent']
    d['Group'] = dic[c]
    tmp = pd.concat([tmp,d])
    
tmp.reset_index(inplace=True)

plt.rcParams['font.sans-serif'] = ['Arial']
plt.figure(figsize=(6,5.2))
ax = sns.barplot(x='index',y='percent',hue='Group',data=tmp,palette=['#E07F6F','#79A323','#1AB5B7'], estimator=np.mean,ci=15, capsize=.05,errwidth=0.8)
#ax.set_xticklabels(labels=['iDC', 'cDC', 'aDC'],rotation=30)
ax.set_xlabel('')
ax.set_ylabel('Immune Cell Proportion')
ax.legend(frameon=False)
ax.set_ylim(0,0.15)
#label_diff(0.02,0.25,0.0775,0.037,0.035,0.1,'*')
barplot_annotate_brackets(-0.28,0.24,'*',0.047,0.047)
barplot_annotate_brackets(0.72,1.24,'*',0.055544,0.055544)
barplot_annotate_brackets(1,1.24,'*',0.049075,0.049075)
barplot_annotate_brackets(1.72,2.24,'*',0.141,0.141)
barplot_annotate_brackets(2,2.24,'**',0.134291,0.134291)

plt.savefig('barplot2-ver.pdf',bbox_inches='tight')
plt.savefig('barplot2-ver.png',dpi=300,bbox_inches='tight')
plt.close()
