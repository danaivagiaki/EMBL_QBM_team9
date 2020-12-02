import pandas as pd
import numpy as np

colorectal = pd.read_csv("./PREDOC_COURSE/QBM/Colorectal_cancer/feat_all.tsv", sep="\t")
colorectal_meta = pd.read_csv("./PREDOC_COURSE/QBM/Colorectal_cancer/meta_all.tsv", sep="\t")

# Remove species that are absent from all samples (i.e. # 0 == ncols)
colorectal.drop(colorectal.loc[(colorectal == 0).astype(int).sum(axis=1) == colorectal.shape[1]].index).reset_index(drop=False).rename(columns={"index":"Species"}).to_csv("./PREDOC_COURSE/QBM/Colorectal_cancer/feat_all_0removed.tsv", sep="\t", header=True, index=False)

# Aggregate the counts on the genus level
taxon = pd.read_csv("./PREDOC_COURSE/QBM/motus_2.5_taxonomy.tsv", sep="\t")
colorectal = colorectal.assign(mOTU_ID = pd.Series(colorectal.index.str.extract(r'((ref|meta)_[a-zA-Z]+_v25_[0-9]+)')[0].values, index=colorectal.index))
colorectal = taxon.xs(["mOTU_ID", "genus"], axis=1).merge(colorectal, how="right", on="mOTU_ID").set_index(colorectal.index)

aggregated = pd.DataFrame(index=colorectal.genus.unique(), columns=colorectal.columns[2:])
for i in aggregated.index:
    aggregated.loc[i] = colorectal.loc[colorectal.genus==i].xs(colorectal.columns[2:], axis=1).sum(axis=0)

# Remove genera that are absent from all samples (i.e. number of 0 in a row == ncols)
aggregated.drop(aggregated.loc[(aggregated == 0).astype(int).sum(axis=1) == aggregated.shape[1]].index, inplace=True)

# Write to file
aggregated.reset_index(drop=False).rename(columns={"index":"Genus"}).to_csv("./PREDOC_COURSE/QBM/Colorectal_cancer/feat_genus.tsv", sep="\t", header=True, index=False)

# Get meta-data info; e.g. colorectal_meta.Diabetes.value_counts(dropna=False)
list(map(lambda x: colorectal_meta[x].value_counts(dropna=False), colorectal_meta.columns))

# Check whether for samples without metadata
colorectal.columns[~colorectal.columns.isin(colorectal_meta.Sample_ID)]

