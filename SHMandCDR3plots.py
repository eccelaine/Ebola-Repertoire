#input file headers are: cluster_id,count,mongo_id,heavy_cdr3_length,heavy_percent_id,heavy_mutations,light_cdr3_length,light_percent_id,light_mutations


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

df=pd.read_csv('Clusters-SHM-CDR3.csv')
df.head()

# Make the PairGrid
g = sns.PairGrid(df,x_vars=["count","heavy_cdr3_length","light_cdr3_length"], y_vars=["cluster_id"],
                 height=300, aspect=.07)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, size=70, orient="h",
      palette="ch:s=1,r=-.1,h=1_r", linewidth=9, edgecolor="black")

# Use the same x axis limits on all columns and add better labels
# g.set(xlim=(0, 30), xlabel="Crashes", ylabel="")

# Use semantically meaningful titles for the columns
titles = ["Number in Cluster", "HC CDR3 Length",
          "LC CDR3 Length"]

for ax, title in zip(g.axes.flat, titles):

    # Set a different title for each axes
    ax.set(title=title)

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)

g.savefig("CDR3s.png")
