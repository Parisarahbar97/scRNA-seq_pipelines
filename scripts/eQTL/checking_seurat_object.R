library(Seurat) 

seurat_obj <- readRDS("/data/files/updated_seurat_with_individual_id.rds")

# Check metadata
head(seurat_obj@meta.data)

# A data.frame: 6 × 20
# orig.ident	nCount_originalexp	nFeature_originalexp	label	n_counts	n_genes	celltypist.Human_AdultAged_Hippocampus.conf	batch	total_counts_mt	Sample_ID	total_counts	n_genes_by_counts	pct_counts_mt	sample	scvi.global.0.5_leiden	scvi.unknown.0.5_leiden	scvi.global.1.0_leiden	scvi.unknown.1.0_leiden	cell_type	Individual_ID
# <fct>	<dbl>	<int>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<fct>	<chr>
# S8A_S14_mapped_AAACCCAAGAGGGTAA-1	S8A	3079	1861	unknown	3079.0	1861	0.9999988261461078	S8A_S14_mapped	13.0	S8A_S14_mapped	3079.0	1861	0.42221498	S8A_S14_mapped	2	2	6	6	Astrocytes	N1540/16
# S8A_S14_mapped_AAACCCAAGCAAATGT-1	S8A	3949	2504	unknown	3949.0	2504	0.9999997351565508	S8A_S14_mapped	52.0	S8A_S14_mapped	3949.0	2504	1.316789	S8A_S14_mapped	2	2	7	7	Astrocytes	N1540/16
# S8A_S14_mapped_AAACCCAAGCGTTGTT-1	S8A	1498	1205	unknown	1498.0	1205	0.999321529376342	S8A_S14_mapped	10.0	S8A_S14_mapped	1498.0	1205	0.66755676	S8A_S14_mapped	3	3	4	4	Microglia	N1540/16
# S8A_S14_mapped_AAACCCAAGGATCATA-1	S8A	3147	2096	unknown	3147.0	2096	0.9999977238413239	S8A_S14_mapped	50.0	S8A_S14_mapped	3147.0	2096	1.5888147	S8A_S14_mapped	5	5	10	10	OPCs	N1540/16
# S8A_S14_mapped_AAACCCAAGGGTCAAC-1	S8A	4471	2537	unknown	4471.0	2537	0.8965772004027295	S8A_S14_mapped	17.0	S8A_S14_mapped	4471.0	2537	0.38022813	S8A_S14_mapped	0	0	0	0	mOli	N1540/16
# S8A_S14_mapped_AAACCCAAGGTAGACC-1	S8A	1339	1089	unknown	1339.0	1089	0.4070970920699443	S8A_S14_mapped	13.0	S8A_S14_mapped	1339.0	1089	0.9708738	S8A_S14_mapped	1	1	2	2	mGC	N1540/16

# Check the count matrix
seurat_obj[["originalexp"]]@counts[1:5, 1:5]

# 5 x 5 sparse Matrix of class "dgCMatrix"
#          S8A_S14_mapped_AAACCCAAGAGGGTAA-1 S8A_S14_mapped_AAACCCAAGCAAATGT-1
# A1BG                                     .                                 .
# A1BG-AS1                                 .                                 .
# A1CF                                     .                                 .
# A2M                                      .                                 .
# A2M-AS1                                  .                                 .
#          S8A_S14_mapped_AAACCCAAGCGTTGTT-1 S8A_S14_mapped_AAACCCAAGGATCATA-1
# A1BG                                     .                                 .
# A1BG-AS1                                 .                                 .
# A1CF                                     .                                 .
# A2M                                      .                                 .
# A2M-AS1                                  .                                 .
#          S8A_S14_mapped_AAACCCAAGGGTCAAC-1
# A1BG                                     .
# A1BG-AS1                                 .
# A1CF                                     .
# A2M                                      .
# A2M-AS1                                  .

# Verify Individual_ID distribution
table(seurat_obj$Individual_ID)

# doublet  N1001/08  N1024/14  N1066/18  N1084/22  N1174/19  N1216/06  N1220/20 
# 233095        35      9426       349     24526     18080      3030       813 
# N1229/20  N1264/16  N1270/22  N1286/22  N1311/20  N1351/20  N1373/16  N1379/15 
# 3364     33887      1886     32584     14082       853       181     29893 
# N1402/21  N1462/18  N1505/17  N1521/17  N1540/16  N1551/14  N1783/14  N2012/12 
# 1423      3203     33784     27532     21162     19675     26006     34319 
# N216/10   N236/02   N248/17   N366/21    N47/13   N478/22   N485/14   N508/10 
# 27642     20864      1203      2703       228       742     21460       530 
# N535/19   N590/16    N71/10   N737/17   N783/19 N792/2017   N812/20  N91/2018 
# 1714     27130     35599         1        40        87       655       278 
# N946/16   N969/17  N97/2017   N970/19   N989/16   N993/15 
# 22021     28242      1674      5849      3689       467 

# Verify cell-type distribution
table(seurat_obj$cell_type)


            #             Astrocytes             CA1_neurons           CA2-4_neurons 
            #                  86741                    5916                    8798 
            #             CA_neurons           Cajal-Retzius           ChoroidPlexus 
            #                  24179                    4999                     718 
            #            Endothelial               Ependymal            GABA_neurons 
            #                  18999                    5673                   19318 
            #              Microglia                    OPCs Subcul_EntorCtx_neurons 
            #                 207363                   58179                    4005 
            #                   imGC                     mGC                    mOli 
            #                   3250                   43909                  283959 