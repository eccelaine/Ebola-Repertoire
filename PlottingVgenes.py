#input files are with headers: heavy v,percentofrep,sample

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

df_heavy=pd.read_csv('Vgenes-Combined=HC.csv')
df_light=pd.read_csv('Vgenes-Combined=LC.csv')
df_heavy.head()

fig_dims=(100,25)
fig ,ax=plt.subplots(figsize = fig_dims)


ax = sns.stripplot(x='heavy v', y='percentofrep', hue='sample', data=df_heavy, linewidth=8, size=50, edgecolor="black", palette=['#b18dfc','#c4c4c4'])
ax.set_yticklabels(ax.get_yticks(),size=50)
ax.set_xticklabels(ax.get_xticks(),size=50)

plt.show()
fig.savefig("heavyVs.png")

fig_dims=(100,25)
fig ,ax=plt.subplots(figsize = fig_dims)


ax = sns.stripplot(x='lightv', y='percentofrep', hue='sample', data=df_light, linewidth=8, size=50, edgecolor="black", palette=['#dbc2f2','#e6e6e6'])
ax.set_yticklabels(ax.get_yticks(),size=50)
ax.set_xticklabels(ax.get_xticks(),size=50)

plt.show()
fig.savefig("lightVs.png")
