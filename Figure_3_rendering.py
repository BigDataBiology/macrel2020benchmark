from six import StringIO
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

data = pd.read_excel('Homology_test.data_source.Figure3.xlsx')
data.rename(columns={'NAMP/AMP': 'Unbalance'}, inplace=True)

figsize = 9,15
fig,axes = plt.subplots(2,2,figsize=figsize)

sns.lineplot(x='Unbalance', y='Accuracy', hue='Model', ax=axes[0][0], data=data)
axes[0][0].set_ylabel('Accuracy',fontsize=20)
axes[0][0].tick_params('both',labelsize = 20)
axes[0][0].set_title('a)',fontsize=20,loc='left')

sns.lineplot(x='Unbalance', y='Precision', hue='Model', ax=axes[0][1], data=data)
axes[0][1].set_ylabel('Precision',fontsize=20)
axes[0][1].tick_params('both',labelsize = 20)
axes[0][1].set_title('b)',fontsize=20,loc='left')
axes[0][1].set_ylim([0.9,1.01])

sns.lineplot(x='Unbalance', y='Sensitivity', hue='Model', ax=axes[1][0], data=data)
axes[1][0].set_ylabel('Sensitivity',fontsize=20)
axes[1][0].tick_params('both',labelsize = 20)
axes[1][0].set_title('c)',fontsize=20,loc='left')
axes[1][0].set_xlabel("NAMP/AMP",fontsize=20)
#sns.despine(fig, trim=True)

sns.lineplot(x='Unbalance', y='Specificity', hue='Model', ax=axes[1][1], data=data)
axes[1][1].set_ylabel('Specificity',fontsize=20)
axes[1][1].tick_params('both',labelsize = 20)
axes[1][1].set_title('d)',fontsize=20,loc='left')
axes[1][1].set_xlabel("NAMP/AMP",fontsize=20)
axes[1][1].set_ylim([0.9,1.01])
fig.tight_layout()
fig.savefig('Fig3.svg')

