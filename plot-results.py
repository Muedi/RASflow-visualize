# plot rasflows results.
#%%
import seaborn as sns 
import numpy as np
import pandas as pd
import os

# %%
tpm_dict = {}
for root, dirs, files in os.walk("output/neutrophiles/trans/tpmFile/"):
    for name in files:
        if name.endswith((".tsv")):
            tpm = pd.read_csv("output/neutrophiles/trans/tpmFile/"+name, sep="\t", header=None, index_col=0)
            tpm_dict[name.split(".")[0]] = tpm
# %%
tpm_df = pd.concat(tpm_dict).T.stack()
tpm_df.index = tpm_df.index.get_level_values(1)
# %%
# both gene and trans level degs are only filtered by a padj cutoff. no lfc
# dea transcript level
deg_trans = {}
for root, dirs, files in os.walk("output/neutrophiles/trans/dea/DEA/transcript-level/"):
    for name in files:
        if name.startswith(("deg")):
            tpm = pd.read_csv("output/neutrophiles/trans/dea/DEA/transcript-level/" +
                              name, sep="\t", header=0, index_col=0)
            deg_trans[name.split(".")[0]] = tpm

#%%
# dea gene level
deg_gene = {}
for root, dirs, files in os.walk("output/neutrophiles/trans/dea/DEA/gene-level/"):
    for name in files:
        if name.startswith(("deg")):
            tpm = pd.read_csv("output/neutrophiles/trans/dea/DEA/gene-level/" +
                              name, sep="\t", header=0, index_col=0)
            deg_gene[name.split(".")[0]] = tpm

# %%
# pca of tpm
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

data = tpm_df.copy().T
data["Group"] = [
    "KO-PU1_cellline" if "PU1_Neu_v2" in x else "WT" if "WT" in x else "KO-PU1+KO-Ets2" if "creEts2" in x else "KO-PU1_prim" for x in data.index]
pca = PCA(n_components=3)
pca_result = pca.fit_transform(data.drop("Group", axis=1).values)

data['pca-one'] = pca_result[:, 0]
data['pca-two'] = pca_result[:, 1]
data['pca-three'] = pca_result[:, 2]
# prints Variation per Principal component
print('Explained variation per principal component: {}'.format(
    pca.explained_variance_ratio_))

# "D-Plot of the first and second PC"
color_number = len(set(data["Group"]))
sns.color_palette("dark")
pca_plot = plt.figure(figsize=(5, 5))
pca_plot = sns.scatterplot(
    x="pca-one",
    y="pca-two",
    hue="Group",
    palette=sns.color_palette("dark", color_number),
    data=data,
    legend="full",
    alpha=0.8
)
plt.title("PCA of TPM (transcript level) by salmon")
fig = pca_plot.get_figure()
fig.savefig("PCA of TPM transcript level by salmon.png")
#%%
# gene level abundance pca
abundance = {}
for root, dirs, files in os.walk("output/neutrophiles/trans/dea/countGroup/"):
    for name in files:
        if name.endswith(("gene_norm.tsv")):
            tpm = pd.read_csv("output/neutrophiles/trans/dea/countGroup/" +
                              name, sep="\t", header=0, index_col=0)
            abundance[name.split("_")[0]] = tpm
# %%
cols = list(abundance.keys())
df = abundance[cols[0]]
for col in cols[1:]:
    print(col)
    df = df.join(abundance[col])

df = df.fillna(0)

data = df.copy().T
data =  np.log2(data + 1)
data["Group"] = [
    "KO-PU1_cellline" if "PU1_Neu_v2" in x else "WT" if "WT" in x else "KO-PU1+KO-Ets2" if "creEts2" in x else "KO-PU1_prim" for x in data.index]
pca = PCA(n_components=3)
pca_result = pca.fit_transform(data.drop("Group", axis=1).values)

data['pca-one'] = pca_result[:, 0]
data['pca-two'] = pca_result[:, 1]
data['pca-three'] = pca_result[:, 2]
# prints Variation per Principal component
print('Explained variation per principal component: {}'.format(
    pca.explained_variance_ratio_))

# "D-Plot of the first and second PC"
color_number = len(set(data["Group"]))
sns.color_palette("dark")
pca_plot = plt.figure(figsize=(5, 5))
pca_plot = sns.scatterplot(
    x="pca-one",
    y="pca-two",
    hue="Group",
    palette=sns.color_palette("dark", color_number),
    data=data,
    legend="full",
    alpha=0.8
)
plt.title("PCA of gene-level log-norm counts by DESeq2")
fig = pca_plot.get_figure()
fig.savefig("PCA of log-normalized abundance gene level by DESeq2.png")


# %%
# hierarchical clustering of degs

