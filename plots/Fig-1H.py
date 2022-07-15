# -*- coding: utf-8 -*-
# Time : 2021/5/13 13:27
# Author : 张志康

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker

########################处理数据###########################
cor = pd.read_csv("gene_cor_between_mrna_protein.xls", sep='\t')
gene_cor = dict(zip(cor["gene"], cor["correlation"]))
term_cor = {}

df = pd.read_excel("作图的通路.xlsx")


def func01(row):
    ids = row["GeneSymbols"].split(";")
    cors = [gene_cor[id] for id in ids]
    term_cor[row["Term_description"]] = cors


df.apply(func01, axis=1)

term_desription = ["", ""]
term_desription.extend(df["Term_description"].to_list())
term_desription = [term_desription[i] for i in range(len(term_desription) - 1, -1, -1)]

###########################绘图###########################
sns.set(color_codes=True)  # 导入seaborn包设定颜色
#
# # sns.set()
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
#
fig = plt.figure(figsize=(8, 16))
ax = fig.add_axes([0.5, 0.04, 0.4, 0.94])
plt.yticks(np.linspace(0, 1.1, 37), term_desription, family="Arial")
plt.xticks([0, 1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6, 1], [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
ax.set_xlabel("Spearman's r", family="Arial")

axs = []
for i in range(34, -1, -1):
    ax_ = fig.add_axes([0.5, 0.04 + 0.94 / 36 * i, 0.4, 0.05])
    ax_.patch.set_facecolor('None')
    ax_.spines['top'].set_visible(False)
    ax_.spines['right'].set_visible(False)
    ax_.spines['left'].set_visible(False)
    ax_.axes.get_yaxis().set_visible(False)
    ax_.set_xlim([-0.5, 1])
    ax_.set_xticks([])
    axs.append(ax_)
#
#生成渐变色
colors = ['#FFA402', '#FCA508', '#FAA60F', '#F7A715', '#F5A91C', '#F2AA22', '#F0AB29', '#EEAC2F', '#EBAE36', '#E9AF3C',
          '#E6B043', '#E4B149', '#E2B350', '#DFB456', '#DDB55D', '#DAB663', '#D8B86A', '#D6B970', '#D3BA77', '#D1BB7D',
          '#CEBD84', '#CCBE8A', '#CABF91', '#C7C097', '#C5C29E', '#C2C3A4', '#C0C4AB', '#BEC5B1', '#BBC7B8', '#B9C8BE',
          '#B6C9C5', '#B4CACB', '#B2CCD2', '#AFCDD8', '#ADCEDF']
# #
[sns.kdeplot(term_cor[term_desription[34 - i]], shade=True, ax=axs[i], color=colors[i], edgecolor="black", alpha=1) for i in range(35)]

# cax, _= mpl.colorbar.make_axes(ax1, shrink=0.1,location='top')
color_bar_ax = fig.add_axes([0.18, 0.948, 0.3, 0.1])
color_bar_ax.axis("off")
colors = list(reversed(colors))
cmp = mpl.colors.ListedColormap(colors)
norm = mpl.colors.BoundaryNorm(np.linspace(-0.02,0.38,35), 35)
cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmp, norm=norm), ax=color_bar_ax, orientation='horizontal',pad=0.1)
tick_locator = ticker.MaxNLocator(nbins=4)  # colorbar上的刻度值个数
cb.locator = tick_locator
cb.update_ticks()

color_bar_ax.text(0.4, 0.0000001, "Medium", family="Arial")
#
# plt.show()
plt.savefig("cor_distribution_pathway.png",dpi=300)
plt.savefig("cor_distribution_pathway.pdf")
