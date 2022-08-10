import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_path = sys.argv[1]
output_path = input_path.replace(".txt",".png")

# read input data
df = pd.read_csv(input_path, sep="\t")

# create plot
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
sns.barplot(data=df, x="length", y="loci", color="lightgray", ax=axes[0,0])
sns.barplot(data=df, x="length", y="sites", color="lightgray", ax=axes[0,1])
sns.barplot(data=df, x="length", y="filtered", color="lightgray", ax=axes[1,0])
sns.barplot(data=df, x="length", y="variant", color="lightgray", ax=axes[1,1])
fig.suptitle("populations summary")
plt.tight_layout()
plt.savefig(output_path)
