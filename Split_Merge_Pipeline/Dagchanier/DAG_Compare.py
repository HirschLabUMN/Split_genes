import pandas as pd


# brohammer split merge calls
def gene_types(file):
    """Takes in a gene key filename and split our pandas df with gene
    type listed as well as other types for split and merged"""
    gene_key = pd.read_csv(file, sep="\t", index_col=False, header=None, comment='#')
    # add empty column for gene call
    gene_key.columns = ['Q_chrom',
                        'Q_gene',
                        'Q1_order',
                        'Q2_order',
                        'S_chrom',
                        'S_gene',
                        'S1_order',
                        'S2_order',
                        'pval',
                        'Syn',
                        ]
    # find duplicated genes one and use for checking later to save time
    Q_dup = gene_key["Q_gene"]
    Many_2_One = Q_dup[Q_dup.duplicated(keep=False)]
    S_dup = gene_key["S_gene"]
    One_2_Many = S_dup[S_dup.duplicated(keep=False)]

    counter = 0
    gene_type = []
    for row in gene_key.itertuples():
        if len(Many_2_One[Many_2_One.isin([row.Q_gene])]) > 0:
            gene_type.append("Many to one")
        elif len(One_2_Many[One_2_Many.isin([row.S_gene])]) > 0:
            gene_type.append("One_2_Many")
        else:
            gene_type.append("One_2_One")
            counter += 1
    gene_key['gene_type'] = gene_type
    print(counter)
    return(gene_key)


def gk_one_2_One(file1, file2):
    """Takes in two gene key filename and will designate one to one matches
    that are in reciprocal agreement betweeen the files"""
    gk1 = gene_types(file1)
    gk2 = gene_types(file2)

    one_to_one = 0
    for row in gk1.itertuples():
        if row.gene_type == "One_2_One":
            Q_gene = row.Q_gene
            S_gene = row.S_gene
            idx = gk2[gk2.Q_gene.isin([S_gene])].index.tolist()
            if idx:
                # make sure that they are also one to one
                if gk2.at[idx[0], "gene_type"] == "One_2_One":
                    gk2_S_gene = gk2.at[idx[0], "S_gene"]
                    if gk2_S_gene == Q_gene:
                        one_to_one += 1
    print(one_to_one)
    return(one_to_one)


# exmpale run commands
JMM_GK121_B73_Mo17 = gk_one_2_One("B73_PH207.filtered.aligncoords",
                                  "PH207_B73.filtered.aligncoords")


