# -*- coding: utf-8 -*-
# Time : 2021/5/12 9:52
# Author : 张志康

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy.stats import norm
import seaborn as sns
from scipy.interpolate import interp1d

sns.set(color_codes=True)  # 导入seaborn包设定颜色
# sns.set_style(style='white')

plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

fig = plt.figure(figsize=(12, 6))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

df = pd.read_csv("gene_cor_between_mrna_protein.xls", sep='\t')
print(df)

q_lt = df[df["adjust_pvalue"] <= 0.05]["correlation"]
print(df[(df["adjust_pvalue"] <= 0.05) & (df["correlation"] > 0)]["correlation"])
q_gt = df[df["adjust_pvalue"] > 0.05]["correlation"]
# print(q_gt)
print(2775/6454)


n, bins, patches = ax.hist([q_lt / 94, q_gt / 94], bins=100, stacked=True, color=["#FFA402", "#AFCDDA"])
# # print(n[0])
sns.distplot(df["correlation"] / 94, bins=94, hist=False, rug=False,
             kde=True, color="black")  # kde=False关闭核密度分布,rug表示在x轴上每个观测上生成的小细条（边际毛毯）
# sns.distplot(df["correlation"], bins=100, hist=True, rug=False, kde=True)  # kde=False关闭核密度分布,rug表示在x轴上每个观测上生成的小细条（边际毛毯）

# 中位数
ax.axvline(x=0.18771694816540477 / 100, c="black", ls="--")

ax.set_ylabel("Counts")
ax.set_xlabel("Spearman's r")

plt.xticks([-0.005, -0.0025, 0, 0.0025, 0.005, 0.0075, 0.01], [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])

#################曲线拟合##################
# n, bins, patches = plt.hist(df["correlation"],bins=100)
#
# x = []
# for i in range(len(bins)-1):
#     x.append(np.mean([bins[i],bins[i+1]]))
#
# print(x)
# print(n)
#
# xnew = np.linspace(min(x),max(x),300)
# func = interp1d(x,n,kind='linear')
# ynew = func(xnew)
# plt.plot(xnew,ynew)
##########################################

# poly = np.polyfit(x, n, deg=2)
# y_value = np.polyval(poly, x)
# plt.plot(bins[0:-1], y_value)

# sns.distplot(df["correlation"][df["correlation"]>-0.0094938677191534],kde=False, rug=False,color="red",bins=66)#kde=False关闭核密度分布,rug表示在x轴上每个观测上生成的小细条（边际毛毯）
# sns.distplot(df["correlation"][df["correlation"]<=-0.0094938677191534],kde=False, rug=False,color="blue",bins=34)#kde=False关闭核密度分布,rug表示在x轴上每个观测上生成的小细条（边际毛毯）

# sns.distplot(x, bins=20, kde=False, rug=True)#设置了20个矩形条

#自定义图例得方法
rect1 = plt.Rectangle((0.0075, 170), 0.0005, 5, color="#FFA402")
rect2 = plt.Rectangle((0.0075, 160), 0.0005, 5, color="#AFCDDA")

legend = ax.legend((rect1,rect2),("FDR < 0.05","FDR > 0.05"))
# ax.add_patch(rect1)
# ax.add_patch(rect2)
# ax.text(0.0082, 170, "FDR<0.05", family="Arial", fontsize=10)
# ax.text(0.0082, 160, "FDR>0.05", family="Arial", fontsize=10)

plt.show()
# plt.savefig('cor_distribution.png', dpi=300)
# plt.savefig('cor_distribution.pdf')
