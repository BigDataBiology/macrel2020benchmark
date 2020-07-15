import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
matplotlib.style.use('default')

data = pd.read_excel('Homology_test.data_source.Figure3.xlsx')
data.rename(columns={'NAMP/AMP': 'Unbalance'}, inplace=True)
data['Model'] = data.Model.map({
    'AMP Scanner Re-trained': 'AMP Scanner',
    'Homology Classification': 'BLAST (homology)',
    'Macrel Re-trained': 'Macrel',
    'iAMP-2L Re-trained': 'iAMP-2L',
    }.get)

matplotlib.rcParams['legend.fontsize'] = 10.0
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rc('ytick', labelsize=9)
matplotlib.rc('xtick', labelsize=9)

figsize = 4.5,6
fig,axes = plt.subplots(2, 2, figsize=figsize, sharex=True)

sns.lineplot(x='Unbalance', y='Accuracy', hue='Model', ax=axes[0][0], data=data)
axes[0][0].set_ylabel('Accuracy')
axes[0][0].set_title('a)', loc='left')

sns.lineplot(x='Unbalance', y='Precision', hue='Model', ax=axes[0][1], data=data)
axes[0][1].set_ylabel('Precision')
axes[0][1].set_title('b)', loc='left')
axes[0][1].set_ylim([0.9,1.01])

sns.lineplot(x='Unbalance', y='Sensitivity', hue='Model', ax=axes[1][0], data=data)
axes[1][0].set_ylabel('Sensitivity')
axes[1][0].set_title('c)', loc='left')
axes[1][0].set_xlabel("NAMP/AMP")

sns.lineplot(x='Unbalance', y='Specificity', hue='Model', ax=axes[1][1], data=data)
axes[1][1].set_ylabel('Specificity')
axes[1][1].set_title('d)', loc='left')
axes[1][1].set_xlabel("NAMP/AMP")
axes[1][1].set_ylim([0.9,1.01])
sns.despine(fig, trim=True)
fig.savefig('Fig3.svg')

