import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# read input data
df = pd.read_csv("gstacks/summary_gstacks.txt", sep="\t")

# create plot
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sns.barplot(data=df, x="length", y="loci", color="lightgray", ax=axes[0])
sns.barplot(data=df, x="length", y="reads", color="lightgray", ax=axes[1])
sns.barplot(data=df, x="length", y="coverage", color="lightgray", ax=axes[2])
fig.suptitle("gstacks summary")
fig.subplots_adjust(wspace=0.4)
plt.savefig('gstacks/summary_gstacks.png')
