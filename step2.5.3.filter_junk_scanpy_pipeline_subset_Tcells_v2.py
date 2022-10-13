import scanpy as sc
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
#
mpl.rcParams['pdf.fonttype'] = 42
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
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Tcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Tcell_v2' +'/'

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
adata= adataall[adataall.obs['Defined_Cell_Type'].isin(['T cells','NK cells'])]
adata_raw=adata.copy()

#
#prepeocess the data
####
####
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata.layers["counts"] = adata.X.copy()
#adata= adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['nCount_RNA', 'percent_mt_merge','CC_Difference'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#batch effect
sce.pp.harmony_integrate(adata, 'Sample')
#
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50,use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.tl.louvain(adata,resolution=0.8)
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
#
#density plot
sc.tl.embedding_density(adata, basis='umap', groupby='age_group')

sc.pl.embedding_density(adata, basis='umap', key='umap_density_age_group', group='all',hspace=0.5,ncols=1 ,save='_age_density.pdf')

#
#
marker_genes = ['CD8A','CD8B','CD4',
               'FCGR3A', 'GZMK', 'KLF2', 'GZMH', 'TRDC', 'XCL1', 'TYROBP','LAYN', 
'HSPH1', 'SLC4A10', 'IL7R', 'CD69', 'SELL', 'GZMA', 'IL32', 'MAL', 'CCR8', 'TNFRSF18','ISG15'
'MCM5', 'TOP2A']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_individuals.mainlineage.pdf')
#
#

save_file = DATA_output_stem+'adata_Tcell_v2.h5ad'
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
    'Tcell': ['CD3D', 'CD3E', 'CD3G'],
    'NK': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1','FCGR3A','FGFBP2','TYROBP','XCL1','XCL2','FCER1G'],
    'CD4-Tcell':['IL7R','CD4'],
    'CD8-Tcell':['CD8A','CD8B'],
    'CD4-T-Eff-Memory':['ZNF683','CD69','ANXA1','S100A4','CCR6','LTB','GPR183','CD40LG','CXCR3','CXCR6','CX3CR1','KLRG1'],
    'Naive': ['TCF7','SELL','LEF1','CCR7'],
    'Exhuasted':['IFI6','LAG3','TIGIT','PDCD1','HAVCR2','CTLA4','BTLA'],
    'Cytotoxic': ['IL2','GZMA','PRF1','GZMB','GZMK','IFNG','GZMH'],
    'Treg':['IL2RA','FOXP3','IKZF2'],
    'Tfh':['FTH1','MAF','CXCR5','CXCL13'],
    'MAIT':['SLC4A10','RORC', 'RORA','ZBTB16','CEBPD','KLRB1','IFNGR1','CCL20','NCR3'],
    #'Th17':['IRF4','CREM','NR4A2'],
    #'Th1':['STAT4','IL12RB2'],
    #'Th2': ['GATA3','STAT6','IL4'],
    'gamma delta T':['TRDC','TRGC2','TRGC1','TRDV1'],
    'Proliferation':['STMN1', 'MKI67', 'TOP2A', 'CDK1'],
    'other':['KLF2','LAYN', 'HSPH1', 'SLC4A10',  'IL32', 'MAL', 'CCR8', 'TNFRSF18','ISG15','ENTPD1','TOX']
}
#
sc.tl.dendrogram(adata, groupby='louvain')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain',swap_axes=True,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Tcell.markers.pdf')
              
#
#########calculate signature scores from Gueguen et al science advances
#########
#########calculate signature scores from Gueguen et al science advances
#########
pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/genesignature_Gueguen_SVds.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)

#
signatures = ['Effectors', 'Naive', 'Central_memory', 'Tfh',
       'Tregs', 'GDTcells', 'Interferon', 'Cycling', 'Stemness',
       'Negative_ICP ', 'Positive_ICP', 'Good_response_to_ICP',
       'Bad_response_to_ICP', 'Exhausted_progenitors', 'Terminal_exhausted',
       'DNA_repair', 'ER_stress']
sc.pl.umap(adata, color=signatures,hspace=0.5,ncols=1,color_map= 'plasma',save='_signature_scores.pdf')

#
save_file = DATA_output_stem+'adata_Tcell_v2_degs_signature_added.h5ad'
adata.write_h5ad(save_file)

#
########### start sub clustering CD8
###########
###########
#######
#######start cluster CD8
adata=sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_added.h5ad')

#CD8
sc.tl.louvain(adata, restrict_to=('louvain', ['1']), 
              resolution=0.2, key_added='louvain_recluster')
sc.tl.louvain(adata, restrict_to=('louvain_recluster', ['1,0']), 
              resolution=0.3, key_added='louvain_recluster')

sc.tl.louvain(adata, restrict_to=('louvain_recluster', ['2']), 
              resolution=0.2, key_added='louvain_recluster')
#CD4
sc.tl.louvain(adata, restrict_to=('louvain_recluster', ['7']), 
              resolution=0.3, key_added='louvain_recluster')
              
adata.obs['louvain_recluster'].cat.categories

####rename the cluster
cluster_dict={'0':'0', '1,0,0':'1', '1,0,1':'2', '1,0,2':'1', '1,1':'3', 
              '2,0':'4','2,1':'5', '3':'6', '4':'7', '5':'8', '6':'9', '7,0':'10',
              '7,1':'11', '8':'12', '9':'13', '10':'14', '11':'15', '12':'16', '13':'17', '14':'18'}
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
# scale and store results in layer
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.tl.dendrogram(adata, groupby='louvain_recluster_renamed')
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5,groupby= 'louvain_recluster_renamed', use_raw=False, swap_axes=True,
                                   vmin=-3, vmax=3, cmap='bwr',min_logfoldchange=1,layer='scaled',
                                  save='_top5_genes_louvain_recluster_renamed.pdf')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain_recluster_renamed',swap_axes=True,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Tcell.markers_louvain_recluster_renamed.pdf')
              
save_file = DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed.h5ad'
adata.write_h5ad(save_file)

#

marker_genes = ['CD8A','CD8B','CD4','FOXP3','XCL1','TYROBP','SLC4A10','KLRB1','RORC','RORA','TRDC','TRDV1','KLRC1',
               'FCGR3A','CX3CR1','KLRG1', 'GZMK', 'KLF2', 'GZMH',  'CXCR6','LAYN', 'ZNF683',
 'IL7R', 'CD69', 'SELL', 'GZMA', 'IL32', 'MAL', 'CCR8', 'TNFRSF18','ISG15','IFIT1','IFIT3',
'HSPH1','ENTPD1','CD40LG','HSPA1A','HSPA1B','CCR6','CXCR5','CXCL13',
'MCM5', 'TOP2A','MKI67']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_individuals.mainlineage_2.pdf')

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

other_markers= ['EPCAM','MARCO','CD68','RAMP2','DCN','MGP']
sc.pl.umap(adata, color=other_markers,hspace=0.5,ncols=1,color_map= 'plasma'
           ,save='_marker_doublet_clusters.pdf')
           
#draw again
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    'Endothelial':['MGP','RAMP2'],
    'Epithelial': ['EPCAM','CDH1'],
    #'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Myeloid': ['CD68', 'MARCO',],
    'NK': [ 'KLRF1','XCL1','TYROBP'],
    'Tcell': ['CD3D', 'CD3E']#,
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}
sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='louvain_recluster_renamed', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_reduced.pdf'
                )

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


#######
#######rename cluster based on function
cluster_dict={'0':'CD4-Naive-CCR7','1':'CD8-TRM-ZNF683','2':'CD8-gammadeltaT-TRDV1',
                '3':'CD8-TerminalllyExhausted-LAYN', '4':'CD8-TEM-GZMK', 
                '5':'CD8-TEM-CX3CR1', '6':'CD4-TEM-HSPA1A', '7':'CD4-Tregs-FOXP3', '8':'NK cell-FCGR3A','9':'CD4-TRM-CXCR6',
                '10':'CD4-MAIT-RORC','11':'CD8-MAIT-SLC4A10','12':'CD4-TfhExhuasted-CXCL13','13':'NK cell-XCL1',
                '14':'Doublet-EPCAM','15':'CD4/CD8-Cycling-MKI67', '16':'Doublet-CD68','17':'CD4/CD8-IFN-ISG15','18':'FewMarker-MGP'}
#####
adata.obs['Defined_Cell_Subtype'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)
#
#######rename cluster
cluster_dict={'0':'C1-CD4-CCR7',
                '1':'C10-CD8-ZNF683',
                '2':'C12-CD8-TRDV1',
                '3':'C11-CD8-LAYN', 
                '4':'C8-CD8-GZMK', 
                '5':'C9-CD8-CX3CR1',
                '6':'C2-CD4-HSPA1A', 
                '7':'C5-CD4-FOXP3', 
                '8':'C13-NK-FCGR3A',
                '9':'C3-CD4-CXCR6',
                '10':'C6-CD4-RORC',
                '11':'C7-CD8-SLC4A10',
                '12':'C4-CD4-CXCL13',
                '13':'C14-NK-XCL1',
                '14':'C17-Doublet-EPCAM',
                '15':'C16-CD4/CD8-MKI67', 
                '16':'C18-Doublet-CD68',
                '17':'C15-CD4/CD8-ISG15',
                '18':'C19-FewMarker-MGP'}
#####
adata.obs['Defined_Cell_Subtype_function'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)
#####
####reorder the subtype annotation
#
adata.obs['Defined_Cell_Subtype_function'].cat.reorder_categories(
    ['C1-CD4-CCR7','C2-CD4-HSPA1A', 'C3-CD4-CXCR6','C4-CD4-CXCL13','C5-CD4-FOXP3', 'C6-CD4-RORC',
     'C7-CD8-SLC4A10','C8-CD8-GZMK', 'C9-CD8-CX3CR1','C10-CD8-ZNF683','C11-CD8-LAYN', 'C12-CD8-TRDV1',
     'C13-NK-FCGR3A','C14-NK-XCL1','C15-CD4/CD8-ISG15','C16-CD4/CD8-MKI67', 
    'C17-Doublet-EPCAM','C18-Doublet-CD68','C19-FewMarker-MGP'
    ], inplace=True)

save_file = DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad'
adata.write_h5ad(save_file)
######
######
#####filter out doublet and bad clusters
#####
adata_good = adata[~adata.obs['Defined_Cell_Subtype_function'].isin(['C17-Doublet-EPCAM','C18-Doublet-CD68','C19-FewMarker-MGP']),:]
save_file = DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm.h5ad'
adata_good.write_h5ad(save_file)
#pd.crosstab(adata.obs['Defined_Cell_Subtype'], adata.obs['louvain_recluster_renamed'])
######
#define age group
cut_labels_3 = ['Young', 'Intermediated', 'Aged']
cut_bins = [20, 49, 60,100]
adata_good.obs['intage'] = adata_good.obs['age'].astype(int)
adata_good.obs['age_group'] = pd.cut(adata_good.obs['intage'], bins=cut_bins, labels=cut_labels_3)
#calculate cyto_exh_nk signature score
pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/cyto_exh_nk_signatures.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata_good,tmp,score_name=i,use_raw=True)
#
#
signatures = ['naive_score',  'cytotoxic_score',  'exhausted_score',  'NK_score']
sc.pl.umap(adata_good, color=signatures,hspace=0.5,ncols=1,color_map= 'plasma',save='_NaiveCytoExhNK_signature_scores.pdf')

adata_good.obs.to_csv(DATA_output_stem+'adata_good.obs.csv')
#
save_file = DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad'
adata_good.write_h5ad(save_file)

#
#########representatve markers for defined cell subtype
marker_genes = ['CD8A','CD8B','CD4','ILR7',
               'CCR7','ZNF683','TRDV1','LAYN','GZMK', 'CX3CR1',
               'HSPA1A','FOXP3','FCGR3A','CXCR6','RORC','SLC4A10',
               'CXCL13','XCL1','MKI67','ISG15']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_individuals.mainlineage_3.pdf')

#########show the annotated clusters
#########
with rc_context({'figure.figsize': (3, 4)}):
    sc.pl.umap(adata, color=['Defined_Cell_Subtype','Defined_Cell_Subtype_function'],
           hspace=0.5,ncols=1,
           #legend_loc='on data',legend_fontsize='small',
           palette=color_CLASS
           ,save='_cluster_anntated_v2.pdf'
          )
##########
##########

######need to selectively draw markers
marker_genes_dict_sel = {
     'Tcell': ['CD3G', 'CD4','CD8A','CD8B'],
     'CD4-Naive-CCR7':['CCR7','TCF7','SELL','LEF1'],
     'CD4-TEM-HSPA1A':['HSPA1A','HSPA1B','DNAJB1','HSPH1','HSPA2','DNAJA1','ANXA1'],
     'CD4-TRM-CXCR6':['CLU','GNA15','PLIN2','CPNE7','MSC','CXCR6','S100A4'],
     'CD4-TfhExhuasted-CXCL13':['MAF','CXCL13','TOX2','BTLA', 'PDCD1','TIGIT','CTLA4'],
     'CD4-Tregs-FOXP3':['IL2RA','FOXP3','TNFRSF4','IL1R2','CD177'],
     'CD4-MAIT-RORC':['IL4I1','CCR6','RORC', 'RORA','CCL20','CEBPD','IFNGR1'],
     'CD8-MAIT-SLC4A10':['SLC4A10','ZBTB16','TMIGD2','LST1','IFNGR1','NCR3'],
     'CD8-TEM-GZMK':['GZMK','CXCR4','CXCR3','CD44','EOMES','ENC1', 'CCL4L2', 'CCL4', 'ITM2C'], 
     'CD8-TEM-CX3CR1':['CX3CR1', 'ADGRG1', 'NKG7', 'FGFBP2','GZMH','KLRG1','PRF1','LILRB1','GNLY'], 
     'CD8-TRM-ZNF683':['ZNF683','CCL5','KLRC1','FXYD2','HOPX','ITGA1','ITM2C','SPRY1'],
     'CD8-TerminalllyExhausted-LAYN':['LAYN','VCAM1', 'ENTPD1','HAVCR2', 'LAG3','TNFRSF9','KRT86'], 
     'CD8-gammadeltaT-TRDV1':['TRDV1','TRDC','TRGC2','TRGC1','IKZF2'],
     'NK cell-FCGR3A':['FGFBP2','FCGR3A','KLRD1','TYROBP','KLRF1','CXCR2', 'CXCR1'],
     'NK cell-XCL1':['NCAM1', 'XCL1','XCL2','KLRC1','FCER1G'],
     'CD4/CD8-IFN-ISG15':['ISG15','IFIT1','MX1','IFI44L','IFI6'],
     'CD4/CD8-Cycling-MKI67':['STMN1', 'MKI67', 'TOP2A', 'CDK1','MCM5']   
}               
            
with rc_context({'figure.figsize': (3, 4)}):
     sc.pl.matrixplot(adata_good, var_names=marker_genes_dict_sel, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',#swap_axes=True,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var', 
                      colorbar_title='mean z-score'
                  ,save='_top_selective_markers_expression_Defined_Cell_Subtype_function.pdf'
                 )
sc.pl.matrixplot(adata_good, var_names=marker_genes_dict_sel, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_vertical_top_selective_markers_expression_Defined_Cell_Subtype_function.pdf'
                 )


sc.pl.matrixplot(adata_good, var_names=marker_genes_dict_sel, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='row scaled\nexpression'
                  ,save='_top_selective_markers_expression_Defined_Cell_Subtype_function_rowscaled.pdf'
                 )
#
#
######functional genesets
marker_genes_fun = {
    'Tcell': ['CD3G', 'CD4','CD8A','CD8B'],
     'Naive marker':['CCR7','TCF7','SELL','LEF1'],
     'Inhibitory receptor':['LAG3', 'PDCD1','TIGIT','CTLA4','HAVCR2','BTLA', 'KLRC1' ],
     'Cytotoxic effector':['IL2','GNLY','PRF1','GZMA','GZMK','GZMB','GZMH','NKG7','IFNG' ],
     'Memory effector':['ANXA1','ANKRD28','IL7R','CD69','CD40LG','S100A4','CXCR3','CXCR6'],
     'Co-stimulatory':['CD28','TNFRSF14','TNFRSF9','ICOS'],
     'Transcription factors':['EOMES','HOPX','TBX21','ZEB2','ZNF683','HIF1A','ID2','TOX','TOX2'],
     'Tregs':['IL2RA','FOXP3','TNFRSF4','IL1R2','CD177'],
     'MAIT cell':['RORC', 'RORA','CCL20','SLC4A10','CCR6','IL4I1'],
     'NK cell':['NCR1','NCAM1','FGFBP2','FCGR3A','KLRD1','TYROBP','KLRF1',
                   'CXCR2', 'CXCR1', 'XCL1','XCL2'],
     'gammadelta_T':['TRDV1','TRDC','TRGC2','TRGC1'],
     'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1','MCM5'],   
     'IFN':['ISG15','IFIT1','MX1','IFI44L','IFI6']
}               


sc.pl.matrixplot(adata_good, var_names=marker_genes_fun, groupby='Defined_Cell_Subtype_function', 
                 dendrogram=False, cmap='RdBu_r',swap_axes=True,
                     layer='scaled',vmin=-2, vmax=2,
                 #standard_scale='var', 
                 colorbar_title='mean z-score'
                 ,save='_vertical_functional_markers_Defined_Cell_Subtype_function.pdf'
                )

sc.pl.matrixplot(adata_good, var_names=marker_genes_fun, groupby='Defined_Cell_Subtype_function', 
                 dendrogram=False, cmap='RdBu_r',swap_axes=False,
                     layer='scaled',vmin=-2, vmax=2,
                 #standard_scale='var',
                 colorbar_title='mean z-score'
                 ,save='_horiz_functional_markers_Defined_Cell_Subtype_function.pdf'
                )

sc.pl.matrixplot(adata_good, var_names=marker_genes_fun, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='row scaled\nexpression'
                  ,save='_horiz_functional_markers_Defined_Cell_Subtype_function_rowscaled.pdf'
                 )

###################################
###################################
#####compare inhibitory and cyto gene exp adding age information

adata_good.obs["aged_cluster"]= adata_good.obs["age_group"].str.cat(adata_good.obs['Defined_Cell_Subtype_function'], sep ="_")


marker_genes_re = {
     'Inhibitory receptor':['LAG3', 'PDCD1','TIGIT','CTLA4','HAVCR2','BTLA', 'KLRC1' ],
     'Cytotoxic effector':['IL2','GZMA','GZMK','GZMB','GZMH','NKG7','PRF1','IFNG','GNLY' ]
}  

adata_good.obs['aged_cluster']=adata_good.obs['aged_cluster'].astype('category')

adata_good.obs['aged_cluster'].cat.reorder_categories(
    ['Aged_C1-CD4-CCR7','Intermediated_C1-CD4-CCR7','Young_C1-CD4-CCR7',
'Aged_C2-CD4-HSPA1A','Intermediated_C2-CD4-HSPA1A','Young_C2-CD4-HSPA1A',
'Aged_C3-CD4-CXCR6','Intermediated_C3-CD4-CXCR6','Young_C3-CD4-CXCR6',
'Aged_C4-CD4-CXCL13','Intermediated_C4-CD4-CXCL13','Young_C4-CD4-CXCL13',
'Aged_C5-CD4-FOXP3','Intermediated_C5-CD4-FOXP3','Young_C5-CD4-FOXP3',
'Aged_C6-CD4-RORC','Intermediated_C6-CD4-RORC','Young_C6-CD4-RORC',
'Aged_C7-CD8-SLC4A10','Intermediated_C7-CD8-SLC4A10','Young_C7-CD8-SLC4A10',
'Aged_C8-CD8-GZMK','Intermediated_C8-CD8-GZMK','Young_C8-CD8-GZMK',
'Aged_C9-CD8-CX3CR1','Intermediated_C9-CD8-CX3CR1','Young_C9-CD8-CX3CR1',
'Aged_C10-CD8-ZNF683','Intermediated_C10-CD8-ZNF683','Young_C10-CD8-ZNF683',
'Aged_C11-CD8-LAYN','Intermediated_C11-CD8-LAYN','Young_C11-CD8-LAYN',
'Aged_C12-CD8-TRDV1','Intermediated_C12-CD8-TRDV1','Young_C12-CD8-TRDV1',
'Aged_C13-NK-FCGR3A','Intermediated_C13-NK-FCGR3A','Young_C13-NK-FCGR3A',
'Aged_C14-NK-XCL1','Intermediated_C14-NK-XCL1','Young_C14-NK-XCL1',
'Aged_C15-CD4/CD8-ISG15','Intermediated_C15-CD4/CD8-ISG15','Young_C15-CD4/CD8-ISG15',
'Aged_C16-CD4/CD8-MKI67','Intermediated_C16-CD4/CD8-MKI67','Young_C16-CD4/CD8-MKI67'
    ], inplace=True)

sc.pl.matrixplot(adata_good, var_names=marker_genes_re, groupby='aged_cluster', 
                 dendrogram=False, cmap='RdBu_r',swap_axes=True,
                     layer='scaled',vmin=-2, vmax=2,
                 #standard_scale='var', 
                 colorbar_title='mean z-score'
                 ,save='_vertical_functional_markers_cytoexh_aged_clustered.pdf'
                )

sc.pl.matrixplot(adata_good, var_names=marker_genes_re, groupby='aged_cluster', 
                 dendrogram=False, cmap='RdBu_r',swap_axes=True,
                     layer='scaled',#vmin=-2, vmax=2,
                 standard_scale='var', 
                 colorbar_title='row scaled\nexpression'
                 ,save='_vertical_functional_markers_cytoexh_aged_clustered_rowscaled.pdf'
                )
############
#2021-06-28
############individually compare the marker expression
############
with rc_context({'figure.figsize': (5, 3)}):
    sc.pl.violin(adata_good, ['TIGIT'], groupby='age_group',
                 swap_axes=True, stripplot=False, inner='box'
                 #, save='_PDCD1_Tcells.compare.pdf'
                ) 

#sc.pl.stacked_violin(adata_good, marker_genes_re, groupby='age_group', dendrogram=False, log=True)


sc.pl.dotplot(adata_good, marker_genes_re, groupby='age_group',  save='_markers_Tcells.compare.age_three_groups.pdf')

#sc.pl.matrixplot(adata_good, marker_genes_re, groupby='age_group', cmap='viridis')





#import numpy as np
#import itertools
#a=np.array(list(marker_genes_fun.values()))
#print(list(itertools.chain.from_iterable(a)))
###
###################
###################
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
mean_zscore = grouped_obs_mean(adata_good, group_key='Defined_Cell_Subtype_function', layer='scaled', gene_symbols=None)
mean_zscore.to_csv(DATA_output_stem+'mean_zscore_allgenes_among_Defined_Cell_Subtype_function_clusters.csv')
#
#
mean_zscore = grouped_obs_mean(adata_good, group_key='aged_cluster', layer='scaled', gene_symbols=None)
mean_zscore.to_csv(DATA_output_stem+'mean_zscore_allgenes_among_aged_cluster.csv')

##########
##########
#########proportion of cluster between age_group
import pandas as pd
from plotnine import *

base_plot= ggplot(adata_good.obs, aes(x='age_group', fill='Defined_Cell_Subtype_function')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(3, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_ageGroup_vs_Defined_Cell_Subtype_function.pdf',path= FIG_output_stem)
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
#2021-05-14
########
adata_good=sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
adata_good
#####
#proportion of  age_group in each cluster
import pandas as pd
from plotnine import *

base_plot= ggplot(adata_good.obs, aes(x='Defined_Cell_Subtype_function', fill='age_group')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = ["#F0E716","#47EFFA", "#EE2EE8"])+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(3, 3))
]
save_as_pdf_pages(plots,filename='stacked_priportion_Defined_Cell_Subtype_function_vs_ageGroup.pdf',path= FIG_output_stem)
#####
#####
#################################################################################################################
#################################################################################################################
#################################################################################################################
#########
#2021-05-15
#########
#########aged degs between three group for each cell type
#########choose to use degs in seurat
#########
os.chdir('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/')
os.getcwd()
# SPECIFY OUTPUT STEMS FOR FIGURES/PATHWAY ANALYSIS
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Tcell_v2/'+'aged_degs' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Tcell_v2/'+'aged_degs' +'/'

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
datadir= '/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/'

adata_good=sc.read(datadir+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
adata_good
#
######
temp_list = adata_good.obs['Defined_Cell_Subtype_function'].unique()
for i in temp_list:
    #tissue,cell_type = analyte.split('.')
    adata_sel = adata_good[adata_good.obs['Defined_Cell_Subtype_function'].isin([i])]
    #find degs between clusters again
    sc.tl.rank_genes_groups(
    adata_sel,'age_group',method='wilcoxon',pts=False,key_added='rank_genes_groups_sub_aged')
    #sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
    #
    cName= re.sub("/", "-", i)
    result = adata_sel.uns['rank_genes_groups_sub_aged']
    groups = result['names'].dtype.names
    pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals_adj','logfoldchanges']}).to_csv(DATA_output_stem + '/degs_celltype_'+cName+'_.csv')
#
#
######
#################################################################################################################
#################################################################################################################
#################################################################################################################
######2021-05-18
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
DATA_output_stem = "/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/" + 'output_Tcell'+'/'

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
adata_good= sc.read('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
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
    cName= re.sub("/", "-", i)
    
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




###########################
###########################
#2021-09-06
####doublet stacked with defined cell type
######
adata= sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad')

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
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    'Endothelial':['MGP','RAMP2'],
    'Epithelial': ['EPCAM','CDH1'],
    #'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Myeloid': ['CD68', 'MARCO',],
    'NK': [ 'KLRF1','XCL1','TYROBP'],
    'Tcell': ['CD3D', 'CD3E']#,
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}
sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='Defined_Cell_Subtype_function', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_reduced_vs_Defined_Cell_Subtype_function.pdf'
                )


#####
#####
#####
#####2021-11-30
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

#####
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Tcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Tcell_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem
#
adata=sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
adata
#
adata.obs.age= adata.obs.age.astype('int')
#adata.obs['age_group'] = 'over80'
#adata.obs.loc[adata.obs['Age'] < 60, 'age_group'] = 'below60'
conditions = [
    (adata.obs['age'] < 40),
    (adata.obs['age'] >= 40) & (adata.obs['age'] < 60),
    (adata.obs['age'] >= 60) & (adata.obs['age'] < 80)
]
values = ['<40 years old', '40-59 years old', '60-79 years old']
adata.obs['age_group2'] = np.select(conditions, values)
adata.obs['age_group2']= adata.obs['age_group2'].astype('category')
#
adata.obs['age_group2'].cat.reorder_categories(
    [
      '<40 years old', '40-59 years old', '60-79 years old'
], inplace=True)
##
sig_dict= {
    'signature':['cytotoxic_score', 'exhausted_score'],
    'Exhuasted':['IFI6','LAG3','TIGIT','PDCD1','HAVCR2','CTLA4','BTLA'],
    'Cytotoxic': ['GZMA','PRF1','GZMB','GZMK','IFNG','GZMH','GNLY','NKG7','KLRK1','CTSW','CST7'],
}
##
##select just T cells
##
adata_tcell= adata[adata.obs['Defined_Cell_Type'].isin(['T cells'])].copy()
adata_tcell

sc.pl.matrixplot(adata_tcell, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='Tcell_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_tcell, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_Tcell_excluedNK_signature_vs_age_group.pdf')
#
#
############
adata_CD8= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C7-CD8-SLC4A10', 'C8-CD8-GZMK',
       'C9-CD8-CX3CR1', 'C10-CD8-ZNF683', 'C11-CD8-LAYN', 'C12-CD8-TRDV1'])].copy()
#
adata_CD4= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C1-CD4-CCR7', 'C2-CD4-HSPA1A', 'C3-CD4-CXCR6', 'C4-CD4-CXCL13',
       'C5-CD4-FOXP3'])].copy()


sc.pl.matrixplot(adata_CD8, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='Tcell_CD8_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_CD8, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_Tcell_CD8_excluedNK_signature_vs_age_group.pdf')
#   
sc.pl.matrixplot(adata_CD4, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='Tcell_CD4_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_CD4, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_Tcell_CD4_excluedNK_signature_vs_age_group.pdf')
#
###NK
##
adata_nkcell= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C13-NK-FCGR3A','C14-NK-XCL1'])].copy()
adata_nkcell

sc.pl.matrixplot(adata_nkcell, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='NKcell_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_nkcell, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_NKcell_excluedNK_signature_vs_age_group.pdf')
#
#
#####
#adata_nkcell.obs['age_ctype']= adata_nkcell.obs['Defined_Cell_Subtype_function'].astype(str)+adata_nkcell.obs['age_group'].astype(str)
adata_nkcell_c13= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C13-NK-FCGR3A'])].copy()
adata_nkcell_c14= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C14-NK-XCL1'])].copy()

sc.pl.matrixplot(adata_nkcell_c13, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='NKcell_indi_c13_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_nkcell_c13, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_NKcell_indi_c13_excluedNK_signature_vs_age_group.pdf')
#
sc.pl.matrixplot(adata_nkcell_c14, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression',save='NKcell_indi_c14_excluedNK_signature_markers_vs_age20range.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_nkcell_c14, ['cytotoxic_score', 
                               'exhausted_score'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_NKcell_indi_c14_excluedNK_signature_vs_age_group.pdf')

              
              
              
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
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Tcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Tcell_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem
#
adata=sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
adata
#
#########
mg =  ['CD3G', 'CD4','CD8A','CD8B',
     'CCR7','TCF7','SELL','LEF1',
     'HSPA1A','HSPA1B','DNAJB1','HSPH1','HSPA2','DNAJA1','ANXA1',
     'CLU','GNA15','PLIN2','CPNE7','MSC','CXCR6','S100A4',
     'MAF','CXCL13','TOX2','BTLA', 'PDCD1','TIGIT','CTLA4',
     'IL2RA','FOXP3','TNFRSF4','IL1R2','CD177',
    'IL4I1','CCR6','RORC', 'RORA','CCL20','CEBPD','IFNGR1',
    'SLC4A10','ZBTB16','TMIGD2','LST1','IFNGR1','NCR3',
    'GZMK','CXCR4','CXCR3','CD44','EOMES','ENC1', 'CCL4L2', 'CCL4', 'ITM2C', 
    'CX3CR1', 'ADGRG1', 'NKG7', 'FGFBP2','GZMH','KLRG1','PRF1','LILRB1','GNLY', 
    'ZNF683','CCL5','KLRC1','FXYD2','HOPX','ITGA1','ITM2C','SPRY1',
     'LAYN','VCAM1', 'ENTPD1','HAVCR2', 'LAG3','TNFRSF9','KRT86', 
     'TRDV1','TRDC','TRGC2','TRGC1','IKZF2',
     'FGFBP2','FCGR3A','KLRD1','TYROBP','KLRF1','CXCR2', 'CXCR1',
     'NCAM1', 'XCL1','XCL2','KLRC1','FCER1G',
     'ISG15','IFIT1','MX1','IFI44L','IFI6',
     'STMN1', 'MKI67', 'TOP2A', 'CDK1','MCM5']   
avg_heatmap(adata, var_names = mg, groupby = "Defined_Cell_Subtype_function", ticksize = 5, figsize = (4,8),
            dendrogram = True,scale = True, use_raw = True, show_gene_labels = True, vmax = 3, vmin = -3, #cmap = "bwr", 
cmap='RdBu_r', save= FIG_output_stem+ '/heatmap_ave_test_markers_top_rank_genes.pdf')
#
#
#
#--------------------------------------------------------------------

#degs among cellsubtype
#    
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
from operator import or_
from functools import reduce
#
sc.settings.set_figure_params(dpi=100, frameon=True, figsize=(3, 4), facecolor='white')

#####
#FIG_output_stem = "./output_step2.5.3.filter_junk_scanpy_pipeline/" + now +'/'
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Tcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Tcell_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem
#
adata=sc.read(DATA_output_stem+'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
adata


adata_rpremove= adata.copy()

prefixes = ["RPS",'RPL','MT-','MTRNR','IGKV']
exclude_genes = reduce(or_, [adata_rpremove.var_names.str.startswith(x) for x in prefixes])
adata_rpremove.var_names[exclude_genes]
#
adata_rpremove = adata_rpremove[:, ~exclude_genes].copy()
adata_rpremove
#
sc.tl.rank_genes_groups(adata_rpremove, 'Defined_Cell_Subtype_function', method='wilcoxon', pts=True, use_raw=False)
#sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
result = adata_rpremove.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:]: result[key][group]
    for group in groups for key in ['names','scores','pvals','pvals_adj','logfoldchanges','pts']}).to_csv(DATA_output_stem + 
                                                                                                          '/adata_rank_genes_groups_Defined_Cell_Subtype_function_T_RPremoved.csv')
#
pd.DataFrame(adata_rpremove.uns['rank_genes_groups']['names']).head(500).to_csv(DATA_output_stem + '/adata_rank_genes_groups_louvain_louvain_subtype_round2_T_RPremoved_top500list.csv')
######