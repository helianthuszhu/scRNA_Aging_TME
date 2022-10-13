import scanpy as sc
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
#
sc.settings.set_figure_params(dpi=100, frameon=True, figsize=(3, 4), facecolor='white')

#
color_CLASS = [
                 '#0000FF', # 0
                 '#FF0000', # 1   
                 '#FF8C00', # 2
                 '#FFD700', # 3  
                 '#808000', # 4
                 '#9ACD32', # 5   
                 '#44F544', # 6
                 '#20B2AA', # 7  
                 '#2F4F4F', # 8
                 '#00BFFF', # 9   
                 '#147ADF', # 10
                 '#191970', # 11 
                 '#6200A4', # 12
                 '#8B008B', # 13  
                 '#FF00FF', # 14
                 '#FF1493', # 15 
                 '#8B4513', # 16
                 '#D2691E', # 17 
                 '#B0C4DE', # 18
                 '#696969', # 19
                 '#FF586E', # 20
                 '#FAEBD7', # 21
                 '#E2C792', # 22
                 '#FFF7A4', # 23
                 '#83715A', # 24
                 '#BC8F8F', # 25
                 '#4C718E', # 26
                 '#519395', # 27
                 '#7FFFD4', # 28  
                 '#0F830F', # 29
                 '#D4FFF0', # 30
                 '#B72222', # 31
                 '#FDBABA', # 32
                 '#1B1010', # 33
                 '#CD57FF', # 34
                 '#F0FFFF', # 35
        
]

#
#
# SPECIFY OUTPUT STEMS FOR FIGURES/PATHWAY ANALYSIS
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Fibroblasts_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Fibroblasts_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

# CREATE FIGURE DIRECTORY IF IT DOES NOT EXIST     
d = os.path.dirname(FIG_output_stem)
if not os.path.exists(d):
        os.makedirs(d)

d = os.path.dirname(DATA_output_stem)
if not os.path.exists(d):
        os.makedirs(d)
#
sc.settings.figdir=FIG_output_stem
#
#
adataall= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data/adata.harmony.overclustered.filtered.CelltypeDefined.withrawcounts.h5ad')
#
adata= adataall[adataall.obs['Defined_Cell_Type'].isin(['Fibroblasts'])]
#adata_raw=adata.copy()
#
#prepeocess the data
####
adata.layers["counts"] = adata.X.copy()
####
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata
adata.layers["logcounts"] = adata.X.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#adata= adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['nCount_RNA', 'percent_mt_merge','CC_Difference'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#batch effect
sce.pp.harmony_integrate(adata, 'Sample')
#
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.tl.louvain(adata,resolution=0.3)
#
##cluster
sc.pl.umap(adata, color=['louvain','Defined_Cell_Type','Sample','Cell_subtype_unimodel'],
           palette=color_CLASS,hspace=0.5,ncols=1,save='_louvain_Defined_Cell_Type_Sample.pdf')
#
#define age group
cut_labels_3 = ['Young', 'Intermediated', 'Aged']
cut_bins = [20, 49, 60,100]
adata.obs['intage'] = adata.obs['age'].astype(int)
adata.obs['age_group'] = pd.cut(adata.obs['intage'], bins=cut_bins, labels=cut_labels_3)
#
#qc and age
sc.pl.umap(adata, color=['intage', 'nCount_RNA', 'nFeature_RNA', 'percent_mt'],color_map= 'plasma', hspace=0.5,ncols=1,save='_intage.qc.pdf')

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['louvain'],
           hspace=0.5,ncols=1,
           legend_loc='on data',legend_fontsize='small',
           palette=color_CLASS
           ,save='louvain_clusters.pdf'
)
#
#density plot
sc.tl.embedding_density(adata, basis='umap', groupby='age_group')

sc.pl.embedding_density(adata, basis='umap', key='umap_density_age_group', group='all',hspace=0.5,ncols=1 ,save='_age_density.pdf')
#
#
marker_genes = ['RGS5','ACTA2','MEF2C','ITGA7','RERGL','ACTG2','SORBS2','COL13A1','COL14A1','ITGA8','NPNT','GSN','MYH11','COL10A1','CFD','MMP2']

sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_individuals.mainlineage.pdf')
#
#
save_file = DATA_output_stem+'adata_Fibroblasts_v2.h5ad'
adata.write_h5ad(save_file)
#
#
#find degs between clusters
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges','pts']}).to_csv(DATA_output_stem + '/adata_rank_genes_groups.csv')
    
#
#
#import numpy as np
#malat1 = adata.var_names.str.startswith('MALAT1')
#mito_genes = adata.var_names.str.match("^MT-|^mt-")
#rp_genes = adata.var_names.str.contains("^RPS|^RPL|^Rps|^Rpl")
#rp_genes = adata.var_names.str.match("^RPS|^RPL|^Rps|^Rpl")
#remove = np.add(rp_genes,mito_genes, malat1)
#keep = np.invert(remove)
#adata_gfil = adata[:,keep]
#sc.tl.rank_genes_groups(adata_gfil, 'louvain', pts=True,method='wilcoxon',use_raw=False,key_added='rank_genes_groups_gfil')
###
###marker expression
###
###

marker_genes_dict = {
    'COL13A1+ matrix FBs': ['COL13A1','TCF21','ITGA8','CXCL14','NPNT'],
    'COL14A1+ matrix FBs': ['COL14A1','GSN','PI16','CYGB','PRRX1'],
    'Myofibroblasts':['ACTA2','MYH11','TAGLN','ACTG2','MYLK','MEF2C', 'ITGA7', 'PDGFA','MYLK','MCAM','MYL9','PLN','PTN','MT2A', 'FAM162B', 'KCNK3','RERGL',
                     'SORBS2'],
    'Smooth muscle cells':['CNN1','SYNPO2','CRYAB','DES'],
    'Mesothelial cells':['UPK3B','MSLN','CALB2','WT1'],
    'Pericytes':['RGS5','CSPG4','ABCC9','KCNJ8','PDGFRB','RRAS', 'RASL12','NOTCH3', 'PTP4A3', 'COL4A2', 'MMP9', 'DES','NDUFA4L2'],
    'Perivascular FB-like cells': ['CYP1B1','APOD'],
    'Lipofibroblast':['FABP4','FABP5','PPARG'],
    'Cancer-associated Fibroblast': ['COL3A1', 'COL1A1', 'COL1A2','COL4A1', 'POSTN', 'FAP', 'TGFB1',
                                    'CTHRC1', 'FN1', 'COL11A1', 'COL10A1'],
    'Immune-modulatory Fibroblast':['OGN', 'CFD', 'CXCL14','APOD','GSN',
                                    'MFAP5', 'PI16'],
    'CAF':['COMP','COL10A1','COL11A1','SERPINE1','PDGFRL','THY1','PDGFRA','MMP2','MMP3','PDPN', 'IGF1', 'WT1','CLDN1'],
    'Normal Fibroblast':['RGCC', 'MACF1', 'A2M', 'DST', 'NPNT']
}

#
sc.tl.dendrogram(adata, groupby='louvain')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain',swap_axes=False,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Fibroblasts.markers.pdf')
              
#
########### start sub clustering
###########
###########
#######
#######start cluster SMC
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2.h5ad')

#SMC
sc.tl.louvain(adata, restrict_to=('louvain', ['3']), 
              resolution=0.2, key_added='louvain_recluster')

adata.obs['louvain_recluster'].cat.categories

####rename the cluster
cluster_dict={'0':'0', '1':'1', 
              '2':'2', '3,0':'3',
              '3,1':'4', '4':'5'
        }
#####
adata.obs['louvain_recluster_renamed'] = (
    adata.obs["louvain_recluster"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)

with rc_context({'figure.figsize': (3, 4)}):
    sc.pl.umap(adata, color=['louvain','louvain_recluster_renamed'],
           hspace=0.5,ncols=1,
           legend_loc='on data',legend_fontsize='small',
           palette=color_CLASS
           ,save='louvain_recluster_renamed.pdf'
)

#find degs between clusters again
sc.tl.rank_genes_groups(adata, 'louvain_recluster_renamed', method='wilcoxon', pts=True,key_added='rank_genes_groups_louvain_recluster_renamed')
#sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
#
result = adata.uns['rank_genes_groups_louvain_recluster_renamed']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges','pts']}).to_csv(DATA_output_stem + '/adata_rank_genes_groups_louvain_recluster_renamed.csv')
    
#
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.tl.dendrogram(adata, groupby='louvain_recluster_renamed')

sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5,groupby= 'louvain_recluster_renamed',
                                   key='rank_genes_groups_louvain_recluster_renamed',
                                   use_raw=False, swap_axes=True,
                                   vmin=-3, vmax=3, cmap='bwr',min_logfoldchange=1,layer='scaled',
                                  save='_top5_genes_louvain_recluster_renamed.pdf')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain_recluster_renamed',swap_axes=False,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Fibroblasts.markers_louvain_recluster_renamed.pdf')

sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='louvain_recluster_renamed', 
                  dendrogram=True, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',
                 #vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='column scaled\nexpression'
                  ,save='_Fibroblasts.markers_louvain_recluster_renamed_scaled.pdf'
                 )

              
save_file = DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed.h5ad'
adata.write_h5ad(save_file)

####################################
#####check potential doublet group
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['scrubletScore_default'],
           hspace=0.5,ncols=1,color_map= 'plasma',
          save='_scrubletScore.pdf')
#########
import pandas as pd
from plotnine import *

base_plot= ggplot(adata.obs, aes(x='louvain_recluster_renamed', fill='scrublet_predict_default')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(4, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_scrublet_predict_default_vs_louvain_recluster_renamed.pdf',path= FIG_output_stem)
###
###check other markers and qc plot
marker_genes_dict = {
    'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    'Endothelial':['PECAM1', 'CLDN5', 'FLT1', 'RAMP2'],
    'Epithelial': ['EPCAM', 'KRT19', 'CDH1', 'KRT18'],
    'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1'],
    'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Myeloid': ['CD68', 'MARCO', 'LYZ','AIF1'],
    'NK': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
    'Tcell': ['CD3D', 'CD3E', 'CD3G', 'TRAC']#,
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}
sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='louvain_recluster_renamed', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_test.pdf'
                )


other_markers= ['CD3E','CD3G','CD3D']
sc.pl.umap(adata, color=other_markers,hspace=0.5,ncols=1,color_map= 'plasma'
           ,save='_marker_doublet_clusters.pdf')
           
#draw again
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    #'Endothelial':['MGP','RAMP2'],
    #'Epithelial': ['EPCAM','KRT19', 'KRT18'],
    'Tcell': ['CD3D', 'CD3E','CD3E'],
    'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1','ACTA2','MYH11'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    #'Myeloid': ['CD68', 'MARCO','LYZ','AIF1'],
    #'NK': [ 'KLRF1','XCL1','TYROBP'],
    #'DC':['LAMP3','FSCN1','CLEC4C','PTCRA','CD207', 'CD1A','CD1C']
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}
sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='louvain_recluster_renamed', dendrogram=True, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_reduced.pdf'
                )

################################
base_plot= ggplot(adata.obs, aes(x='louvain_recluster_renamed',y= 'nFeature_RNA', fill='louvain_recluster_renamed')) 
plots=[base_plot+
geom_boxplot(varwidth = False) + # vary boxes width according to n obs.
  #geom_jitter(alpha = 0.8, width = 0,size=2) + # adds random noise and limit its width
  #facet_wrap(~year) + # divide into 2 panels
  #theme(legend.position = "none") +# remove legend
  #scale_fill_manual(values = c("metastatic"="#74add1","non_metastatic"="#f46d43"))+
  #stat_compare_means(aes(group = group))+
  ggtitle("No. of genes detected")+
  #ylab("% of T cells")+
  theme(figure_size=(4, 3))
]
save_as_pdf_pages(plots,filename='boxplot_nFeature_RNA_vs_louvain_recluster_renamed.pdf',path= FIG_output_stem)

base_plot= ggplot(adata.obs, aes(x='louvain_recluster_renamed',y= 'nCount_RNA', fill='louvain_recluster_renamed')) 
plots=[base_plot+
geom_boxplot(varwidth = False) + # vary boxes width according to n obs.
  #geom_jitter(alpha = 0.8, width = 0,size=2) + # adds random noise and limit its width
  #facet_wrap(~year) + # divide into 2 panels
  #theme(legend.position = "none") +# remove legend
  #scale_fill_manual(values = c("metastatic"="#74add1","non_metastatic"="#f46d43"))+
  #stat_compare_means(aes(group = group))+
  ggtitle("No. of raw Counts")+
  #ylab("% of T cells")+
  theme(figure_size=(4, 3))
]
save_as_pdf_pages(plots,filename='boxplot_nCount_RNA_vs_louvain_recluster_renamed.pdf',path= FIG_output_stem)



############
############MARKER
marker_genes = ['RGS5','COL13A1','COL14A1','MYH11','COL3A1','CTHRC1','CD3E'

                
                #'MEF2C','ITGA7','RERGL','ACTG2','SORBS2',
                #'ITGA8','NPNT','GSN','COL10A1','CFD','MMP2'
               ]

sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_TEST.pdf')
#
##################################



marker_genes_dict_fun = {
'collagens': ['COL1A1', 'COL1A2','COL3A1', 'COL4A1', 'COL4A2','COL5A1', 'COL5A2','COL6A1', 'COL6A2',
           'COL7A1', 'COL8A1', 'COL10A1', 'COL11A1', 'COL12A1', 'COL13A1', 
           'COL14A1', 'COL15A1', 'COL18A1'],
'ECM': ['BGN', 'DCN', 'TAGLN', 'ELN', 'FN1'],

'MMPs':[ 'MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP10', 'MMP11', 'MMP14', 'MMP19'],
'TGFs':[ 'SERPINE1', 'CTHRC1', 'THBS2', 'SULF1', 'TGFB1', 'TGFB2', 'TGFB3', 'COMP', 'INHBA'],
'Angiogenesis':[ 'EGFL6', 'ANGPT2', 'PDGFA', 'PDGFC', 'VEGFA'],
'Contratile': ['ACTA2', 'MYL6', 'MYH9', 'MYH11', 'PLN', 'TPM1', 'TPM2', 'SORBS2'],
'RAS': ['RRAS', 'RASL12', 'RASGRP2'],
'Infla':['CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','CXCL8','CXCL9','CXCL10','CXCL12','CXCL11','CXCL13','CXCL14','CXCL16','CXCL17',
        'IL6','PDGFRA','CCL2'
        ],
    'Antigen-presenting':[ 'CD74','SLPI','HLA-DRA','HLA-DRB1','HLA-DOA','HLA-DQA1',
                         'HLA-DMA','HLA-DMB','HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DRB5'
        ]
}


sc.pl.dotplot(adata, marker_genes_dict_fun, groupby='louvain_recluster_renamed',swap_axes=False,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Fibroblasts.markers_FUNCTION_louvain_recluster_renamed.pdf')

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='louvain_recluster_renamed', 
                  dendrogram=True, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',
                 #vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='column scaled\nexpression'
                  ,save='_Fibroblasts.markers_FUNCTION_louvain_recluster_renamed_scaled.pdf'
                 )

################################
#2021-05-25
###DEFINE each cluster now

####rename the cluster
cluster_dict={'0':'iCAF',
              '1':'Normal-iCAF',
              '2':'CAF',
              '3':'Myofibroblasts',
              '4':'Pericyte',
              '5':'Doublet-Tcell'
             }

#####
adata.obs['Defined_Cell_Subtype'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)
#####
#####add markers
cluster_dict={'0': 'C2-CAF-COL14A1',
              '1': 'C1-CAF-COL13A1',
              '2': 'C3-CAF-CTHRC1',
              '3': 'C4-myCAF-MYH11',
              '4': 'C5-Pericyte-RGS5',
              '5': 'C6-Doublet-CD3E'
             }

adata.obs['Defined_Cell_Subtype_function'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)

#
marker_genes = ['RGS5','COL13A1','COL14A1','MYH11','COL3A1','CTHRC1','CD3E']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma'
          ,save='_marker_indi_refined.pdf'
          )
#
#
#
####reorder the subtype annotation
#
adata.obs['Defined_Cell_Subtype_function']=adata.obs['Defined_Cell_Subtype_function'].astype("category")

adata.obs['Defined_Cell_Subtype_function'].cat.reorder_categories(
    ['C1-CAF-COL13A1','C2-CAF-COL14A1',
              'C3-CAF-CTHRC1',
              'C4-myCAF-MYH11',
              'C5-Pericyte-RGS5',
              'C6-Doublet-CD3E'
    ], inplace=True)


save_file = DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad'
adata.write_h5ad(save_file)
######
######
#########show the annotated clusters
#########
with rc_context({'figure.figsize': (3, 4)}):
    sc.pl.umap(adata, color=['Defined_Cell_Subtype','Defined_Cell_Subtype_function'],
           hspace=0.5,ncols=1,
           #legend_loc='on data',legend_fontsize='small',
           palette=color_CLASS
           ,save='_cluster_anntated_v2.pdf'
          )
#####
#####filter out doublet and bad clusters
#####
adata_good = adata[~adata.obs['Defined_Cell_Subtype_function'].isin(['C6-Doublet-CD3E']),:]
save_file = DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm.h5ad'
adata_good.write_h5ad(save_file)
#pd.crosstab(adata.obs['Defined_Cell_Subtype'], adata.obs['louvain_recluster_renamed'])
######
######
######
###################################################
#####degs based on filtered clusters
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm.h5ad')

#####
sc.tl.rank_genes_groups(adata, 'Defined_Cell_Subtype_function', method='wilcoxon', pts=True,key_added ='Defined_Cell_Subtype_function_group')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key= 'Defined_Cell_Subtype_function_group' )


result = adata.uns['Defined_Cell_Subtype_function_group']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges','pts']}).to_csv(DATA_output_stem + '/adata_rank_genes_Defined_Cell_Subtype_function_group.csv')

#####
sc.tl.dendrogram(adata, groupby='Defined_Cell_Subtype_function')
#sc.pl.rank_genes_groups_matrixplot(adata, n_genes=20, use_raw=False, vmin=-3,  vmax=3, cmap='bwr', layer='scaled', key=          #                                    'Defined_Cell_Subtype_function_group',
#                                  save= '_rank_genes.pdf')
adata.obs['Defined_Cell_Subtype_function'].cat.reorder_categories(
    ['C1-CAF-COL13A1','C2-CAF-COL14A1',
              'C3-CAF-CTHRC1',
              'C4-myCAF-MYH11',
              'C5-Pericyte-RGS5',
    ], inplace=True)

sc.pl.rank_genes_groups_matrixplot(adata, n_genes=20, use_raw=False, swap_axes=True,
                                   #vmin=-3, vmax=3, 
                                   cmap='RdBu_r', layer='scaled',standard_scale='var',colorbar_title='column scaled\nexpression',
                                  key= 'Defined_Cell_Subtype_function_group',
                                  save= '_rank_genes.pdf')
#
marker_genes_dict_fun = {
'collagens': ['COL13A1','COL14A1', 'COL6A2','COL6A1',  
    'COL1A1', 'COL1A2','COL3A1', 'COL5A1', 'COL5A2',
           'COL7A1', 'COL8A1', 'COL10A1', 'COL11A1', 'COL12A1', 'COL15A1', 'COL4A1', 'COL4A2',
           'COL18A1'],
    
'ECM': ['DCN',  'ELN', 'FN1','TAGLN','BGN'],

'MMPs':[  'MMP2','MMP19', 'MMP1', 'MMP3', 'MMP9', 'MMP10', 'MMP11', 'MMP14'],
    
'TGFs':[  'TGFB3','SERPINE1', 'CTHRC1', 'THBS2', 'SULF1', 'TGFB1',  'COMP', 'INHBA','TGFB2'],
    
'Angiogenesis':[ 'VEGFA','PDGFC','PDGFA', 'EGFL6', 'ANGPT2'],
    
'Contratile': ['MYH11', 'PLN', 'SORBS2', 'ACTA2', 'MYL6','TPM2', 'TPM1'],
    
'RAS': ['RRAS', 'RASL12', 'RASGRP2'],
    
'Infla':['CCL2','CXCL14','PDGFRA','CXCL12','CXCL2','CXCL1','CXCL3','CXCL5','CXCL6','CXCL8','CXCL11','CXCL16', 'IL6',
        'CXCL10', 'CXCL9','CXCL13','CXCL17',
       
        ],
    'Antigen-presenting':[ 'CD74','SLPI','HLA-DRA','HLA-DRB1',
                         'HLA-DPA1'
        ]
}


#
marker_genes_dict2 = {
    'C1-CAF-COL13A1':['COL13A1','RGCC', 'MACF1', 'NPNT','CES1'],
    'C2-CAF-COL14A1':['COL14A1','GSN','APOD','OGN','PI16'],
    'C3-CAF-CTHRC1':['CTHRC1', 'POSTN','FAP','TGFB1','COMP'],
    'C4-myCAF-MYH11':['MYH11','ACTA2','ACTG2','MYL9','PLN'],
    'C5-Pericyte-RGS5':['RGS5','PDGFRB','NOTCH3','EGFL6','NDUFA4L2']
}

sc.pl.matrixplot(adata, var_names=marker_genes_dict2, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='row scaled\nexpression'
                  ,save='_markers_refined.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marker_genes_dict2, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='column scaled\nexpression'
                  ,save='_markers_refined_vertical.pdf'
                 )



sc.pl.matrixplot(adata, var_names=marker_genes_dict2, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_refined_zscore.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marker_genes_dict2, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_refined_zscore_vertical.pdf'
                 )

##

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='row scaled\nexpression'
                  ,save='_markers_refined_function.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='column scaled\nexpression'
                  ,save='_markers_refined_vertical_function.pdf'
                 )



sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_refined_zscore_function.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_refined_zscore_vertical_function.pdf'
                 )


###############################
###2021-05-26
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm.h5ad')


###
#define age group
cut_labels_3 = ['Young', 'Intermediated', 'Aged']
cut_bins = [20, 49, 60,100]
adata.obs['intage'] = adata.obs['age'].astype(int)
adata.obs['age_group'] = pd.cut(adata.obs['intage'], bins=cut_bins, labels=cut_labels_3)

####
####add group
adata.obs["aged_cluster"]= adata.obs["age_group"].str.cat(adata.obs['Defined_Cell_Subtype_function'], sep ="_")
adata.obs["aged_cluster"]=adata.obs['aged_cluster'].astype('category')

####
adata.obs['aged_cluster'].cat.reorder_categories(
    ['Aged_C1-CAF-COL13A1', 'Intermediated_C1-CAF-COL13A1','Young_C1-CAF-COL13A1',
     'Aged_C2-CAF-COL14A1' ,'Intermediated_C2-CAF-COL14A1','Young_C2-CAF-COL14A1',
     'Aged_C3-CAF-CTHRC1','Intermediated_C3-CAF-CTHRC1','Young_C3-CAF-CTHRC1',
     'Aged_C4-myCAF-MYH11','Intermediated_C4-myCAF-MYH11' ,'Young_C4-myCAF-MYH11',
     'Aged_C5-Pericyte-RGS5', 'Intermediated_C5-Pericyte-RGS5', 'Young_C5-Pericyte-RGS5'
], inplace=True)
##
##

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='aged_cluster', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_function_aged_cluster.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marker_genes_dict_fun, groupby='aged_cluster', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_function_vertical_aged_cluster.pdf'
                 )
####
sc.tl.rank_genes_groups(adata, 'Defined_Cell_Subtype_function', method='wilcoxon', pts=True,key_added ='Defined_Cell_Subtype_function_group')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key= 'Defined_Cell_Subtype_function_group' )


result = adata.uns['Defined_Cell_Subtype_function_group']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges','pts']}).to_csv(DATA_output_stem + '/adata_rank_genes_Defined_Cell_Subtype_function_group.csv')


save_file = DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked.h5ad'
adata.write_h5ad(save_file)
#####################################################
#####################################################

#####################################################
#####################################################
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked.h5ad')

####
#########proportion of cluster between age_group
import pandas as pd
from plotnine import *

base_plot= ggplot(adata.obs, aes(x='age_group', fill='Defined_Cell_Subtype_function')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(3, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_ageGroup_vs_Defined_Cell_Subtype_function.pdf',path= FIG_output_stem)


#proportion of  age_group in each cluster
import pandas as pd
from plotnine import *

base_plot= ggplot(adata.obs, aes(x='Defined_Cell_Subtype_function', fill='age_group')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = ["#F0E716","#47EFFA", "#EE2EE8"])+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(3, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_Defined_Cell_Subtype_function_vs_ageGroup.pdf',path= FIG_output_stem)


#############
###calculate mean z score between clusters
import pandas as pd
import numpy as np
#
def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names
    #
    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )
    #
    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out
    
#
#
mean_zscore = grouped_obs_mean(adata, group_key='Defined_Cell_Subtype_function', layer='scaled', gene_symbols=None)
mean_zscore.to_csv(DATA_output_stem+'mean_zscore_allgenes_among_Defined_Cell_Subtype_function_clusters.csv')
#
#
mean_zscore = grouped_obs_mean(adata, group_key='aged_cluster', layer='scaled', gene_symbols=None)
mean_zscore.to_csv(DATA_output_stem+'mean_zscore_allgenes_among_aged_cluster.csv')




#################################################################################################################
#################################################################################################################
#################################################################################################################
######2021-05-26
######RSS for TFs for each age
######
####
import re
import scanpy as sc
import pandas as pd
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.rss import regulon_specificity_scores
#
import os
DATA_output_stem = "/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/" + 'output_Fibroblasts'+'/'

print(DATA_output_stem)

# CREATE FIGURE DIRECTORY IF IT DOES NOT EXIST     
d = os.path.dirname(DATA_output_stem)
if not os.path.exists(d):
        os.makedirs(d)
#        
#load auc matrix
outputdata= '/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/output/'

auc_mtx = pd.read_csv(outputdata+ 'auc.csv', index_col=0)
#load Tcell data
#
adata_good= sc.read('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Fibroblasts_v2/adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked.h5ad')
#change cell type / with -
#adata_good.obs['Defined_Cell_Subtype_function']= re.sub("/", "-", adata_good.obs['Defined_Cell_Subtype_function'])
#use a loop to extract each cell type 

temp_list = adata_good.obs['Defined_Cell_Subtype_function'].unique()
for i in temp_list:
    #tissue,cell_type = analyte.split('.')
    adata_sel = adata_good[adata_good.obs['Defined_Cell_Subtype_function'].isin([i])]
    ##subset auc matrix with the same order as adata obs
    auc_mtx_sel= auc_mtx.iloc[pd.Index(auc_mtx.index).get_indexer(adata_sel.obs_names)]
    auc_mtx_sel

    ####RSS calculation
    rss = regulon_specificity_scores(auc_mtx_sel, adata_sel.obs.age)
    #rss.head()
    ####
    #cName= re.sub("/", "-", i)
    cName= i
    
    #datadir= DATA_output_stem+ cName +'/'
    # CREATE FIGURE DIRECTORY IF IT DOES NOT EXIST     
    #d = os.path.dirname(datadir)
    #if not os.path.exists(d):
    #    os.makedirs(d)
    #
    #
    rss.to_csv(DATA_output_stem+ cName + '_rss.csv')
    adata_sel.obs.to_csv(DATA_output_stem+ cName+'_adata.obs.csv')
    auc_mtx_sel.to_csv(DATA_output_stem+ cName+'_auc_sel.csv')


#
################################################################
################################################################
#2021-07-02
###########markers comparson specific
###########
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked.h5ad')
####
adata_sel= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C3-CAF-CTHRC1']),:]
marker_genes_dict_fun = {
'collagens': ['COL1A1', 'COL1A2','COL3A1', 'COL5A1', 'COL5A2','COL6A1', 'COL6A2',
           'COL7A1', 'COL8A1', 'COL10A1', 'COL11A1', 'COL12A1',  'COL15A1', 'COL18A1'],
    'MMPs':[ 'MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP10', 'MMP11', 'MMP14', 'MMP19'],
'TGFs':[ 'SERPINE1', 'CTHRC1', 'THBS2', 'SULF1', 'TGFB1', 'TGFB2', 'TGFB3', 'COMP', 'INHBA'],
}
sc.pl.dotplot(adata_sel, marker_genes_dict_fun, 'age_group', title='C3-CAF-CTHRC1', save= '_coll_MMP_TGF_C3-CAF-CTHRC1.pdf')

#
adata_sel= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C5-Pericyte-RGS5']),:]
marker_genes_dict_fun = {
'collagens': ['COL4A1', 'COL4A2', 'COL18A1'],
'Angiogenesis':[ 'EGFL6', 'ANGPT2', 'PDGFA', 'PDGFC', 'VEGFA']
}
sc.pl.dotplot(adata_sel, marker_genes_dict_fun, 'age_group', title='C5-Pericyte-RGS5', save= '_coll_angio_C5-Pericyte-RGS5.pdf')
#
adata_sel= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C4-myCAF-MYH11']),:]
marker_genes_dict_fun = {
'Contratile': ['ACTA2', 'MYL6', 'MYH9', 'MYH11', 'PLN', 'TPM1', 'TPM2', 'SORBS2'],
}
sc.pl.dotplot(adata_sel, marker_genes_dict_fun, 'age_group', title='C4-myCAF-MYH11', save= '_contratile_C4-myCAF-MYH11.pdf')
#
adata_sel= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C2-CAF-COL14A1']),:]
marker_genes_dict_fun = {
'Infla':['CCL2','CXCL14','PDGFRA','CXCL12','CXCL2','CXCL1','CXCL3','CXCL8','CXCL16', 'IL6'],
'Antigen-presenting':[ 'CD74','SLPI','HLA-DRA','HLA-DRB1', 'HLA-DPA1'
        ]
}
sc.pl.dotplot(adata_sel, marker_genes_dict_fun, 'age_group', title='C2-CAF-COL14A1', save= '_infl_antigenP_C2-CAF-COL14A1.pdf')
#
adata_sel= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C1-CAF-COL13A1']),:]
marker_genes_dict_fun = {
'Infla':['CCL2','CXCL14','PDGFRA','CXCL12','CXCL2','CXCL1','CXCL3','CXCL8','CXCL16', 'IL6'],
'Antigen-presenting':[ 'CD74','HLA-DRA','HLA-DRB1', 'HLA-DPA1'
        ]
}
sc.pl.dotplot(adata_sel, marker_genes_dict_fun, 'age_group', title='C1-CAF-COL13A1', save= '_infl_antigenP_C1-CAF-COL13A1.pdf')

#######################################


############################################
####2021-07-02
####
#########calculate signature scores from gsea
####fibroblast function signatures
#########
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked.h5ad')

pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/fibroblast_function_signature_gsea.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)

#
save_file = DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked_fibfunctionsignatureadded.h5ad'
adata.write_h5ad(save_file)

adata.obs.to_csv(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked_fibfunctionsignatureadded.csv')

####
###
###########################
###########################
#2021-09-08
####doublet stacked with defined cell type
######

import scanpy as sc
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
sc.settings.set_figure_params(dpi=100, frameon=True, figsize=(3, 4), facecolor='white')

# SPECIFY OUTPUT STEMS FOR FIGURES/PATHWAY ANALYSIS
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Fibroblasts_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Fibroblasts_v2' +'/'

sc.settings.figdir=FIG_output_stem

adata= sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad')

import pandas as pd
from plotnine import *

base_plot= ggplot(adata.obs, aes(x='Defined_Cell_Subtype_function', fill='scrublet_predict_default')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(4, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_scrublet_predict_default_vs_Defined_Cell_Subtype_function.pdf',path= FIG_output_stem)

#######marker matrix
#draw again
###check other markers and qc plot
#draw again
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    #'Endothelial':['MGP','RAMP2'],
    #'Epithelial': ['EPCAM','KRT19', 'KRT18'],
    'Tcell': ['CD3D', 'CD3E','CD3E'],
    'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1','ACTA2','MYH11'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    #'Myeloid': ['CD68', 'MARCO','LYZ','AIF1'],
    #'NK': [ 'KLRF1','XCL1','TYROBP'],
    #'DC':['LAMP3','FSCN1','CLEC4C','PTCRA','CD207', 'CD1A','CD1C']
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}

sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='Defined_Cell_Subtype_function', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_reduced_vs_Defined_Cell_Subtype_function.pdf'
                )
            

#
#####2021-12-01
#####new age group
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
#
sc.settings.set_figure_params(dpi=100, frameon=True, figsize=(3, 4), facecolor='white')
    
    
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Fibroblasts_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Fibroblasts_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem    
#
adata= sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked_fibfunctionsignatureadded.h5ad')
adata
#
sig_dict= {
    'signature':['HALLMARK_INFLAMMATORY_RESPONSE', 'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION', 'KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION', 'HALLMARK_ANGIOGENESIS', 'REACTOME_COLLAGEN_FORMATION']
}
#
sc.pl.matrixplot(adata, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression', 
                 save='Stroma_fib_signature_markers_vs_age20range.pdf')
with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata, ['HALLMARK_INFLAMMATORY_RESPONSE', 'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION', 'KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION', 'HALLMARK_ANGIOGENESIS', 'REACTOME_COLLAGEN_FORMATION'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_fib_signature_vs_age_group.pdf')


###################
###################ave_heatmap
######***********************avg_heatmap for each cluster
###
## Function for manual scaling and averaging expression per cell type in Heatmap
## Adapted heatmap function from scanpy
from scanpy.plotting._anndata import _plot_categories_as_colorblocks, _prepare_dataframe
from scanpy.plotting._anndata import _reorder_categories_after_dendrogram, _check_var_names_type
from scanpy.plotting._anndata import _plot_dendrogram
from scanpy.plotting._anndata import _plot_categories_as_colorblocks, _plot_colorbar
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import zscore

def avg_heatmap(adata, var_names, groupby = None, log = False, use_raw = False, num_categories = 7, order = None,
                gene_symbols = None, var_group_labels = None, var_group_positions = None, layer = None,
                show_gene_labels = None, scale = True, save = None, vmax = None, vmin = None, cmap = "viridis",
                dendrogram = True, colorbar_width = 0.2, categorical = True, figsize = (10, 5), ticksize = 10):

    categories, obs_tidy = _prepare_dataframe(adata, var_names, groupby, use_raw, log, num_categories,
                                                  gene_symbols = gene_symbols, layer = layer)
    var_group_labels = None
    var_group_positions = None
    var_names, var_group_labels, var_group_positions = _check_var_names_type(var_names,
                                                                             var_group_labels, var_group_positions)
    ## Inserted for mean per Cluster
    obs_means = obs_tidy.groupby(groupby).mean()
    
    if scale == True:
        for x in obs_means.columns:
            obs_means[x] = (obs_means[x] - obs_means[x].mean()) / obs_means[x].std()
        
    # get categories colors:
    groupby_colors = adata.uns[groupby + "_colors"]

    # obs_tidy = obs_tidy.sort_index()
    if order is not None:
        obs_means.index = obs_means.index.reorder_categories(order)
    else:
        obs_means = obs_means.sort_index()


    if show_gene_labels is None:
        if len(var_names) <= 50:
            show_gene_labels = True
        else:
            show_gene_labels = False
            #logg.warning(
            print('Gene labels are not shown when more than 50 genes are visualized. '
                 'To show gene labels set `show_gene_labels=True`')

    dendro_height = 0
    groupby_height = 0.13 if categorical else 0
    if figsize is None:
        if show_gene_labels:
            heatmap_height = len(var_names) * 0.18
        else:
            heatmap_height = 4
        width = 10
        height = heatmap_height + dendro_height + groupby_height  # +2 to account for labels
    else:
        width, height = figsize
        heatmap_height = height - (dendro_height + groupby_height)

    height_ratios = [dendro_height, heatmap_height, groupby_height]

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        width_ratios = [width, 0.14, colorbar_width]
    else:
        width_ratios = [width, 0, colorbar_width]

    fig = plt.figure(figsize=(width, height))
    axs = gridspec.GridSpec(nrows=3, ncols=3, wspace=0.25 / width,
                            hspace=0.3 / height,
                            width_ratios=width_ratios,
                            height_ratios=height_ratios)

    # plot heatmap
    heatmap_ax = fig.add_subplot(axs[1, 0])

    if order is not None:
        obs_means = obs_means.reindex(order, axis = 0)
    
    im = heatmap_ax.imshow(obs_means.T.values, aspect = 'auto', cmap = cmap, 
                           vmin = vmin, vmax = vmax)#, **kwds)              ## adapted for mean
    heatmap_ax.set_xlim(-0.5, obs_means.shape[0] - 0.5)                     ## adapted for mean
    heatmap_ax.set_ylim(obs_means.shape[1] - 0.5, - 0.5)                    ## adapted for mean

    heatmap_ax.tick_params(axis='x', bottom=False, labelbottom=False)
    heatmap_ax.set_xlabel('')
    heatmap_ax.grid(False)
    if show_gene_labels:
        heatmap_ax.tick_params(axis='y', labelsize='small', length = 1)
        heatmap_ax.set_yticks(np.arange(len(var_names)))
        heatmap_ax.set_yticklabels(var_names, rotation=0, size = ticksize)
        heatmap_ax.set
    else:
        heatmap_ax.tick_params(axis='y', labelleft=False, left=False)

    if categorical:
        groupby_ax = fig.add_subplot(axs[2, 0])
        groupby = obs_means.index.name
        groupby_cmap = ListedColormap(groupby_colors, groupby + '_cmap')
        norm = BoundaryNorm(np.arange(groupby_cmap.N+1)-.5, groupby_cmap.N)

        # determine groupby label positions such that they appear
        # centered next/below to the color code rectangle assigned to the category
        value_sum = -0.5  ## changed from 0
        ticks = []  # list of centered position of the labels
        labels = []
        label2code = {}  # dictionary of numerical values asigned to each label

        for code, (label, value) in enumerate(obs_means.index.value_counts(sort = False).iteritems()):
        
            ticks.append(value_sum + (value / 2))
            labels.append(label)
            value_sum += value
            label2code[label] = code
        groupby_ax.grid(False)

        groupby_ax.imshow(np.matrix([label2code[lab] for lab in obs_means.index]), aspect='auto', 
                          cmap = groupby_cmap, norm = norm)
        if len(labels) > 1:
            groupby_ax.set_xticks(ticks)
            if max([len(x) for x in labels]) < 3:
                # if the labels are small do not rotate them
                rotation = 0
            else:
                rotation = 90
            groupby_ax.set_xticklabels(labels, rotation=rotation)

        # remove x ticks, y ticks and labels
        groupby_ax.tick_params(axis='x', bottom = False, labelsize = ticksize)
        groupby_ax.tick_params(axis='y', left = False, labelleft = False)

        # remove surrounding lines
        groupby_ax.spines['right'].set_visible(False)
        groupby_ax.spines['top'].set_visible(False)
        groupby_ax.spines['left'].set_visible(False)
        groupby_ax.spines['bottom'].set_visible(False)

        groupby_ax.set_xlabel(groupby)

        # add lines to main heatmap
        line_positions = np.cumsum(obs_means.index.value_counts(sort=False))[:-1]
        heatmap_ax.vlines(line_positions - 0.505, -1, len(var_names) + 1, lw = 1)
    

    # plot colorbar
    _plot_colorbar(im, fig, axs[1, 2])
    
    #plt.show()
    plt.savefig(save)
    #plt.close()

    
    
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.pyplot import rc_context
mpl.rcParams['pdf.fonttype'] = 42
#
sc.settings.set_figure_params(dpi=100, frameon=True, figsize=(3, 4), facecolor='white')

#####
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Fibroblasts_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Fibroblasts_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem
#
adata=sc.read(DATA_output_stem+'adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped_generanked_fibfunctionsignatureadded.h5ad')
adata
#
#########
mg =  [ 'COL13A1','RGCC', 'MACF1', 'NPNT','CES1',
    'COL14A1','GSN','APOD','OGN','PI16',
    'CTHRC1', 'POSTN','FAP','TGFB1','COMP',
    'MYH11','ACTA2','ACTG2','MYL9','PLN',
    'RGS5','PDGFRB','NOTCH3','EGFL6','NDUFA4L2']   
avg_heatmap(adata, var_names = mg, groupby = "Defined_Cell_Subtype_function", ticksize = 6,
            figsize = (3,6),
            dendrogram = True,scale = True, use_raw = True, show_gene_labels = True, vmax = 2, vmin = -2, #cmap = "bwr", 
cmap='RdBu_r', save= FIG_output_stem+ '/heatmap_ave_test_markers_top_rank_genes.pdf')

#
#END
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

#
#import numpy as np
#import itertools
#a=np.array(list(marker_genes_fun.values()))
#print(list(itertools.chain.from_iterable(a)))
###
###################
###################


##########
##########

#########
########
#############################################
#############################################signature scores violin plot
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize = (5,7))

sns_plot= sns.violinplot(x=adata_good.obs["age_group"], y=adata_good.obs["naive_score"], palette=["#F0E716","#47EFFA", "#EE2EE8" ],saturation=1)
fig = sns_plot.get_figure()
fig.savefig(FIG_output_stem+"native_score_vs_ageGroup.pdf")

plt.close(fig)
#
#
plt.figure(figsize = (5,7))

sns_plot= sns.violinplot(x=adata_good.obs["age_group"], y=adata_good.obs["cytotoxic_score"], palette=["#F0E716","#47EFFA", "#EE2EE8" ],saturation=1)
fig = sns_plot.get_figure()
fig.savefig(FIG_output_stem+"cytotoxic_score_vs_ageGroup.pdf")

plt.close(fig)
##
plt.figure(figsize = (5,7))

sns_plot= sns.violinplot(x=adata_good.obs["age_group"], y=adata_good.obs["exhausted_score"], palette=["#F0E716","#47EFFA", "#EE2EE8" ],saturation=1)
fig = sns_plot.get_figure()
fig.savefig(FIG_output_stem+"exhausted_score_vs_ageGroup.pdf")

plt.close(fig)
##
plt.figure(figsize = (5,7))

sns_plot= sns.violinplot(x=adata_good.obs["age_group"], y=adata_good.obs["NK_score"], palette=["#F0E716","#47EFFA", "#EE2EE8" ],saturation=1)
fig = sns_plot.get_figure()
fig.savefig(FIG_output_stem+"NK_score_vs_ageGroup.pdf")

plt.close(fig)
########
#############################################
#########

              
              
              