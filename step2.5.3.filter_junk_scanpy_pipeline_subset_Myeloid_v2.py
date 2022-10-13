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
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Myeloidcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Myeloidcell_v2' +'/'

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
adata= adataall[adataall.obs['Defined_Cell_Type'].isin(['Myeloid cells'])]
#adata_raw=adata.copy()

#
#prepeocess the data
####
adata.layers["counts"] = adata.X.copy()
####
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
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
marker_genes = ['FABP4','CD14','FCGR3A','FCN1','CD1C','CD207','STMN1',
                'CLEC9A','LAMP3','LILRA4', 'CCL18','SPP1','VEGFA','SOX11']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_individuals.mainlineage.pdf')
#
#
save_file = DATA_output_stem+'adata_Myeloid_v2.h5ad'
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
    'Monocyte CD14': ['CTSS','FCN1','S100A8','S100A9','LYZ','VCAN'],
    'Non-classical CD16+ Monocyte': ['FCGR3A','CDKN1C','POU2F2','ZNF703'],
    'Classical CD14+ Monocyte':['CD14','FCN1','SERPINB2'],
    'Macrophage':['LGMN','CTSB'],
    'Alveolar-Mac':['MARCO','FABP4','MCEMP1'],
    'Perivascular resident Macrophage': ['LYVE1','PLTP','LILRB5','SELP','FOLR2','SLC40A1'],
    'Tumor-associated Macrophage':['SPP1','SLAMF9','LDHA','VEGFA','SLC2A3','CD81','CCL18','MMP9','CX3CR1'],
    'early stage Macrophage': ['CXCL10','GBP1','CXCL9','CXCL11','CCR2','CCL2'],
    'DC':['CLEC10A','CLEC4C','PTCRA'],
    'Migratory Conventional Dendritic cell':['FSCN1','CCR7','LAMP3','CCL17','CCL19','CCL22','CD40','BIRC3'],
    'CD1c+ DCs (LCs)':['CD1C','ITGAX'],
    'CD141+ DCs-cross-presenting DCs-Conventional DC type 1':['CLEC9A','XCR1'],
    'CD207+CD1a+ LCs': ['CD207','CD1A'],
    'Activated DCs':['CCR7','LAMP3'],
    'pDCs':['IL3RA','CLEC4C','LILRB4','IRF4','LILRA4','TCF4','MZB1'],
    'CD163+CD14+ DCs':['CD163'],
    'Langerhans Cells':['FCER1A','CD1C','CD1E','PKIB','CD1A'],
    'Monocytes Derived DC':['FCGR2B','CCL17','CLEC10A','LST1','IFITM2','LYPD2','CFP','LILRA5'],
    'Granulocytes':['S100A12','G0S2','FCGR3B'],
    'Proliferation':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}

#
sc.tl.dendrogram(adata, groupby='louvain')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain',swap_axes=False,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Myeloid.markers.pdf')
              
#
#########calculate signature scores from Gueguen et al science advances
#########
#########calculate signature scores from Gueguen et al science advances
#########
pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/M1_M2_macrophages_signatures.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)

#
signatures = ['M1','M2','Angiogenesis',	'Phagocytosis']
sc.pl.umap(adata, color=signatures,hspace=0.5,ncols=1,color_map= 'plasma',save='_signature_scores.pdf')

save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_added.h5ad'
save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature.4s_added.h5ad'
adata.write_h5ad(save_file)
#

########### start sub clustering
###########
###########
#######
#######start cluster DC
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_added.h5ad')

#DC
sc.tl.louvain(adata, restrict_to=('louvain', ['2']), 
              resolution=0.2, key_added='louvain_recluster')

adata.obs['louvain_recluster'].cat.categories

####rename the cluster
cluster_dict={'0':'0', '1':'1', 
              '2,0':'2', '2,1':'3',
              '3':'4', '4':'5', 
              '5':'6', '6':'7', '7':'8',
               '8':'9', '9':'10', '10':'11', 
              '11':'12', '12':'13', '13':'14', '14':'15'}
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
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5,groupby= 'louvain_recluster_renamed', use_raw=False, swap_axes=True,
                                   vmin=-3, vmax=3, cmap='bwr',min_logfoldchange=1,layer='scaled',
                                  save='_top5_genes_louvain_recluster_renamed.pdf')

sc.pl.dotplot(adata, marker_genes_dict, groupby='louvain_recluster_renamed',swap_axes=False,
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2,
              save='_Myeloid.markers_louvain_recluster_renamed.pdf')
              
save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed.h5ad'
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


other_markers= ['EPCAM','CD3G','CD3D']
sc.pl.umap(adata, color=other_markers,hspace=0.5,ncols=1,color_map= 'plasma'
           ,save='_marker_doublet_clusters.pdf')
           
#draw again
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    #'Endothelial':['MGP','RAMP2'],
    'Epithelial': ['EPCAM','KRT19', 'KRT18'],
    'Tcell': ['CD3D', 'CD3E'],
    #'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Myeloid': ['CD68', 'MARCO','LYZ','AIF1'],
    #'NK': [ 'KLRF1','XCL1','TYROBP'],
    'DC':['LAMP3','FSCN1','CLEC4C','PTCRA','CD207', 'CD1A','CD1C']
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
############MARKER from science advs lung paper and pan cell research

marsciadv= {
    'MAC-INHBA':['INHBA','IL1RN','CCL14'],
    'MAC-NLRP3':['NLRP3','EREG','IL1B'],
    'MAC-C1QC':['C1QC','C1QA','APOE'],
    'Alveolar-Mac':['MARCO','MCEMP1'],
    'MAC1-Alveolar': ['MRC1','EMP1' ,'EMP3', 'S100A6' ,'SH3BGRL3'],
'MAC2-Alveolar': ['SCD','RBP4' ,'IL17RB' ,'CAMP', 'GCHFR'] ,
'MAC3-Alveolar':['CES1' ,'FABP4', 'SLC19A3' ,'PPARG' ,'INHBA' ],
'MAC4-Pervascular':['LYVE1', 'LILRB5', 'SELP', 'FOLR2' ,'SLC40A1' ,'PLTP'],
'MAC5-anti-inflam':['MMP7', 'TIMP3', 'PLA2G7', 'CHI3L1', 'CTSB' ],
'MAC6-anti-inflam':['TNF', 'AXL', 'HS3ST1', 'RGS1','C3' ,'ADAMDEC1'],
'MAC7-TAM':['SPP1', 'SLAMF9' ,'LDHA', 'VEGFA' ,'SLC2A3', 'CCL18','MMP9','CX3CR1'],
'MAC8-Cycling':['STMN1', 'H2AFZ', 'TUBA1B' ,'PCNA' ],
'MAC9-early-stage':['CXCL9' ,'CXCL10', 'CXCL11', 'GBP1', 'CCL2','CCR2'],
'Mono1-classCD14':['CD14', 'FCN1' ,'S100A12','VCAN' ,'SERPINB2' ],
'Mono2-nonclassCD16':['FCGR3A', 'CDKN1C' ,'POU2F2', 'ZNF703' ],
'Mono-DC':['FCGR2B', 'CCL17', 'CLEC10A' ],
'DC1-CDC':['CLEC9A', 'BATF3', 'IRF8' ,'CPVL', 'CADM1'],
'DC2-CDC':['CD1C' ,'CD1A' ,'CD207' ],
'DC3-Migratory-CDC':['LAMP3' ,'FSCN1' ,'CCR7'] ,
'pDC':['IRF4', 'LILRA4' ,'TCF4' ]#
#Granu': ['G0S2', 'S100A8', 'S100A9','FCGR3B']
}

#sc.pl.matrixplot(adata, var_names=marsciadv, groupby='louvain_recluster_renamed', dendrogram=True, cmap='plasma',swap_axes=True,
 #standard_scale='var', colorbar_title='column scaled\nexpression'
                 #,save='_mainleage_marker_reduced.pdf'
  #              )

sc.pl.matrixplot(adata, var_names=marsciadv, groupby='louvain_recluster_renamed', 
                  dendrogram=True, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_from_sciadvpaper.pdf'
                 )
#

pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/M1_M2_macrophages_signatures.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)

#
signatures = ['M1','M2','Angiogenesis',	'Phagocytosis']
sc.pl.umap(adata, color=signatures,hspace=0.5,ncols=1,color_map= 'plasma',save='_signature_scores.pdf')

save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added.h5ad'
adata.write_h5ad(save_file)


#adata.var.to_csv(DATA_output_stem+'adta.var_gene_names.csv')

marMM= {
    'M1':['IL23A','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7'
],
    'M2':['IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA',
'VEGFB','VEGFC','EGF','CTSA',
'CTSB','CTSC','CTSD','TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B',
'FASLG','TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4','CD163','MRC1','IL10','FOLR2','STAB1','CD163L1','SELP','F13A1','MERTK','AXL','IL1RN','GPNMB','CHI3L1'
],
    'MMPS':['MMP1','MMP7','MMP12'],
    'Others':['PPARG','FABP4','INHBA','ALDH2','EGFL7','CD209','CH25H','LILRB5']

}

sc.pl.matrixplot(adata, var_names=marMM, groupby='louvain_recluster_renamed', 
                  dendrogram=True, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',vmin=-2, vmax=2,
                  #standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_MM.pdf'
                 )
#
#CELL RESEARCH description of LYVE1, and other cell type
#degs check in each clusters


marker_genes = ['FABP4','PPARG','LYVE1','CHI3L1','TNF','AXL', 'INHBA',
                'VEGFA','CCL18','CCL2','CD81','SPP1','MMP9','CX3CR1','CXCL9','CXCL10','TGFA','TGFB1', 'TGFB2','TGFB3',
                'CD86','CD68','TLR2','MSR1','MRC1','CD163'
                #'CD14','CD14','FCGR3A','FCN1','CD1C','CD207','STMN1','PCNA'
               # 'CLEC9A','LAMP3','LILRA4', 'SOX11'
               ]

sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma',save='_marker_TEST.pdf')
#
##################################

adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added.h5ad')
#2021-05-22
###DEFINE each cluster now

####rename the cluster
cluster_dict={'0':'Mac-Alveolar1',
              '1':'Mac-Pervascular',
              '2':'cDC2',
              '3':'DC-Langerhanslike',
              '4':'Mac-TAM1',
              '5':'Mac-Alveolar2',
              '6':'Monocyte-CD14',
              '7':'Mac-Alveolar3',
              '8':'Mac-TAM2',
              '9':'Cycling',
              '10':'Doublet-Tcell',
              '11':'Monocyte-CD16',
              '12':'cDC-Migratory',
              '13':'Doublet-Epi',
              '14': 'pDC',
              '15':'cDC1'
             }

#####
adata.obs['Defined_Cell_Subtype'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)
#####
#####add markers
cluster_dict={'0': 'C3-Mac-FABP4',
              '1': 'C6-Mac-CCL13',
              '2': 'C10-cDC2-CLEC10A',
              '3': 'C13-cDC-LCs-CD207',
              '4': 'C7-Mac-CX3CR1',
              '5': 'C4-Mac-RETN',
              '6': 'C1-Monocyte-CD14',
              '7': 'C5-Mac-MARCO',
              '8': 'C8-Mac-SPP1',
              '9': 'C14-Cycling-MKI67',
              '10': 'C15-Doublet-CD3D',
              '11':'C2-Monocyte-CD16',
              '12':'C11-cDCs-M-LAMP3',
              '13':'C16-Doublet-EPCAM',
              '14':'c12-pDC-LILRA4',
              '15':'C9-cDC1-CLEC9A'
             }

adata.obs['Defined_Cell_Subtype_function'] = (
    adata.obs["louvain_recluster_renamed"]
    .map(lambda x: cluster_dict.get(x, x))
    .astype("category")
)



marker_genes = [ 'FCGR3A','CD14','FABP4','CCL13','CX3CR1','RETN','MARCO',
                'SPP1','MKI67','CLEC9A','CLEC10A','CD207','CD1C',
                'LAMP3','LILRA4','CD3D','EPCAM'
               ]

sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1,color_map= 'plasma'
          ,save='_marker_indi_refined.pdf'
          )

#
#
####reorder the subtype annotation
#
adata.obs['Defined_Cell_Subtype_function'].cat.reorder_categories(
    ['C1-Monocyte-CD14','C2-Monocyte-CD16', 'C3-Mac-FABP4','C4-Mac-RETN','C5-Mac-MARCO',
        'C6-Mac-CCL13','C7-Mac-CX3CR1','C8-Mac-SPP1','C9-cDC1-CLEC9A','C10-cDC2-CLEC10A',
        'C11-cDCs-M-LAMP3','c12-pDC-LILRA4','C13-cDC-LCs-CD207','C14-Cycling-MKI67',
        'C15-Doublet-CD3D','C16-Doublet-EPCAM'
    ], inplace=True)

save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined.h5ad'
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
adata_good = adata[~adata.obs['Defined_Cell_Subtype_function'].isin(['C15-Doublet-CD3D','C16-Doublet-EPCAM']),:]
save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm.h5ad'
adata_good.write_h5ad(save_file)
#pd.crosstab(adata.obs['Defined_Cell_Subtype'], adata.obs['louvain_recluster_renamed'])
######
######
######
###################################################
#####degs based on filtered clusters
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm.h5ad')

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
    ['C1-Monocyte-CD14','C2-Monocyte-CD16', 'C3-Mac-FABP4','C4-Mac-RETN','C5-Mac-MARCO',
        'C6-Mac-CCL13','C7-Mac-CX3CR1','C8-Mac-SPP1','C9-cDC1-CLEC9A','C10-cDC2-CLEC10A',
        'C11-cDCs-M-LAMP3','c12-pDC-LILRA4','C13-cDC-LCs-CD207','C14-Cycling-MKI67'
    ], inplace=True)

sc.pl.rank_genes_groups_matrixplot(adata, n_genes=20, use_raw=False, swap_axes=True,
                                   #vmin=-3, vmax=3, 
                                   cmap='RdBu_r', layer='scaled',standard_scale='var',colorbar_title='column scaled\nexpression',
                                  key= 'Defined_Cell_Subtype_function_group',
                                  save= '_rank_genes.pdf')
#
#
marker_genes_dict2 = {
    'C1-Monocyte-CD14':['S100A8','S100A9','VCAN','CD14','SERPINB2'],
    'C2-Monocyte-CD16':['FCGR3A','CDKN1C','POU2F2','ZNF703','LYPD2'], 
    
    'C3-Mac-INHBA':['FABP4','INHBA','PPARG', 'ALDH2', 'TREM1'],
    'C4-Mac-PPARG':['LYZ','GRN','RETN','TNFSF13','PGD'],
    'C5-Mac-RBP4': ['MCEMP1','S100A6','MARCO','TSPO','FBP1'],
    'C6-Mac-CCL13':['CCL13','LILRB5','FOLR2','SLC40A1','GPNMB'],
    
    
    'C7-Mac-CX3CR1' :['CX3CR1','CCL3','CCL4','CCL4L2','KLF6'],
    'C8-Mac-SPP1': ['SPP1','LGALS1','S100A10','BRI3','AQP9'],
    'C9-cDC1-CLEC9A': ['CLEC9A', 'CPVL','XCR1','BATF3','SNX3'],
    'C10-cDC2-CLEC10A': ['CLEC10A','FCGR2B','CD1C','PEA15','CD1E'],
    'C11-cDCs-M-LAMP3' : ['LAMP3', 'FSCN1', 'CCR7' ,'CCL19','BIRC3'],
    'c12-pDC-LILRA4':['LILRA4','TCF4', 'IRF4','GZMB','IRF7'],
    'C13-cDC-LCs-CD207': ['CD207','CD1A','CD1E','LTB','FCGBP'],
    'C14-Cycling-MKI67':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
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
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_refined_vertical_zscore.pdf'
                 )


save_file=DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked.h5ad'
adata.write_h5ad(save_file)

####
###############################
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked.h5ad')

#
adata=adata[adata.obs['Defined_Cell_Subtype_function'].str.contains('Mac')]
adata
#
marMM= {
    'M1':['IL6','CXCL2','CXCL10','KYNU',
        'IL23A','IL1B','CXCL11',
          'CCL5','IRF1','CD86','CD40','IDO1','CCR7','CXCL9','SOCS3','TNF','IL1A','IRF5'
],
    'M2':['CCL20','EGF','CD276','IL1RN','CHI3L1',
          'MMP19','MMP7','CCL22','MMP12','MMP9','MMP14','MMP1','CLEC7A','FASLG','TNFSF12','IRF4',
       'IL10','AXL','TNFSF8',
        'IL4R','CCL4','CCL17','VEGFA','VEGFB','VEGFC','CTSD','CTSB','F13A1','MERTK','GPNMB','CD163','CCL13','CCL18','LYVE1','FOLR2','STAB1','CD163L1','SELP',
        'TGFB3',
        'WNT7B','VTCN1','CCL24','MRC1','MSR1','CTSC','TGFB1','TGFB2'
],

    'Others':['CD209','LILRB5','EGFL7','CH25H','ALDH2', 'PPARG','FABP4','INHBA']

}
#
sc.pl.matrixplot(adata, var_names=marMM, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_M1M2.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marMM, groupby='Defined_Cell_Subtype_function', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_M1M2_vertical.pdf'
                 )



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
    [ 'Aged_C3-Mac-FABP4', 'Intermediated_C3-Mac-FABP4', 'Young_C3-Mac-FABP4', 
      'Aged_C4-Mac-RETN','Intermediated_C4-Mac-RETN', 'Young_C4-Mac-RETN',
      'Aged_C5-Mac-MARCO', 'Intermediated_C5-Mac-MARCO','Young_C5-Mac-MARCO', 
      'Aged_C6-Mac-CCL13', 'Intermediated_C6-Mac-CCL13','Young_C6-Mac-CCL13', 
      'Aged_C7-Mac-CX3CR1','Intermediated_C7-Mac-CX3CR1', 'Young_C7-Mac-CX3CR1',
      'Aged_C8-Mac-SPP1','Intermediated_C8-Mac-SPP1','Young_C8-Mac-SPP1'
    ], inplace=True)
##
##
sc.pl.matrixplot(adata, var_names=marMM, groupby='aged_cluster', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=False,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_M1M2_aged_cluster.pdf'
                 )

sc.pl.matrixplot(adata, var_names=marMM, groupby='aged_cluster', 
                  dendrogram=False, cmap='RdBu_r',swap_axes=True,
                      layer='scaled',#vmin=-2, vmax=2,
                  standard_scale='var',
                  colorbar_title='mean z-score'
                  ,save='_markers_M1M2_vertical_aged_cluster.pdf'
                 )
####
#####################################################
#####################################################

#####################################################
#####################################################
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked.h5ad')

#define age group
cut_labels_3 = ['Young', 'Intermediated', 'Aged']
cut_bins = [20, 49, 60,100]
adata.obs['intage'] = adata.obs['age'].astype(int)
adata.obs['age_group'] = pd.cut(adata.obs['intage'], bins=cut_bins, labels=cut_labels_3)

####
####add group
adata.obs["aged_cluster"]= adata.obs["age_group"].str.cat(adata.obs['Defined_Cell_Subtype_function'], sep ="_")
####
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
######2021-05-24
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
DATA_output_stem = "/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/" + 'output_Myeloidcell'+'/'

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
adata_good= sc.read('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Myeloidcell_v2/adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked.h5ad')
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
############################################
####2021-07-01
####
#########calculate signature scores from pan cancer myloid paper
####DCs signatures
#########
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked.h5ad')

pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/DC_signature_pan.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)

#
save_file = DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked_DC_signatures_added.h5ad'
adata.write_h5ad(save_file)

adata.obs.to_csv(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked_DC_signatures_added.csv')

####
###

###########################
###########################
#2021-09-07
####doublet stacked with defined cell type
######
adata= sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined.h5ad')

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
marker_genes_dict = {
    #'Bcell': ['CD79A', 'IGHM', 'IGHG3', 'IGHA2'],
    #'Endothelial':['MGP','RAMP2'],
    'Epithelial': ['EPCAM','KRT19', 'KRT18'],
    'Tcell': ['CD3D', 'CD3E'],
    #'Fibroblasts': ['DCN', 'COL1A1', 'COL1A2', 'THY1'],
    #'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Myeloid': ['CD68', 'MARCO','LYZ','AIF1'],
    #'NK': [ 'KLRF1','XCL1','TYROBP'],
    'DC':['LAMP3','FSCN1','CLEC4C','PTCRA','CD207', 'CD1A','CD1C']
    #'Cycling':['STMN1', 'MKI67', 'TOP2A', 'CDK1']
}


sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='Defined_Cell_Subtype_function', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression'
                 ,save='_mainleage_marker_reduced_vs_Defined_Cell_Subtype_function.pdf'
                )
            


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
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Myeloidcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Myeloidcell_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem    
#
adata= sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked_DC_signatures_added.h5ad')
adata
#
sig_dict= {
    'signature':['M1','M2','Angiogenesis','Phagocytosis'],
     'M1':['IL6','CXCL2','CXCL10','KYNU',
        'IL23A','IL1B','CXCL11',
          'CCL5','IRF1','CD86','CD40','IDO1','CCR7','CXCL9','SOCS3','TNF','IL1A','IRF5'
],
    'M2':['CCL20','EGF','CD276','IL1RN','CHI3L1',
          'MMP19','MMP7','CCL22','MMP12','MMP9','MMP14','MMP1','CLEC7A','FASLG','TNFSF12','IRF4',
       'IL10','AXL','TNFSF8',
        'IL4R','CCL4','CCL17','VEGFA','VEGFB','VEGFC','CTSD','CTSB','F13A1','MERTK','GPNMB','CD163','CCL13','CCL18','LYVE1','FOLR2','STAB1','CD163L1','SELP',
        'TGFB3',
        'WNT7B','VTCN1','CCL24','MRC1','MSR1','CTSC','TGFB1','TGFB2'
]#,

  #  'Others':['CD209','LILRB5','EGFL7','CH25H','ALDH2', 'PPARG','FABP4','INHBA']
}

##select just Mac cells
##
adata_mac= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C3-Mac-FABP4', 'C4-Mac-RETN',
       'C5-Mac-MARCO', 'C6-Mac-CCL13', 'C7-Mac-CX3CR1', 'C8-Mac-SPP1'])].copy()
adata_mac

sc.pl.matrixplot(adata_mac, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression', 
                 save='Mac_signature_markers_vs_age_group.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_mac, ['M1','M2','Angiogenesis','Phagocytosis'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_Mac_signature_vs_age_group.pdf')
####
####select DC
sig_dict= {
    'signature':['ActivatedDC', 'MigratoryDC']
}
####
adata_dc= adata[adata.obs['Defined_Cell_Subtype_function'].isin(['C9-cDC1-CLEC9A', 'C10-cDC2-CLEC10A', 'C11-cDCs-M-LAMP3',
       'c12-pDC-LILRA4', 'C13-cDC-LCs-CD207'])].copy()
adata_dc

sc.pl.matrixplot(adata_dc, sig_dict, 'age_group', dendrogram=False, cmap='Blues', 
                 standard_scale='var', colorbar_title='column scaled\nexpression', 
                 save='dc_signature_markers_vs_age_group.pdf')

with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.violin(adata_dc, ['ActivatedDC', 'MigratoryDC'], 
                 groupby='age_group', stripplot=False, inner='box',
                save= 'violin_dc_signature_vs_age_group.pdf')
    
    
    
#
#
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
FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_Myeloidcell_v2' +'/'
DATA_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'data_Myeloidcell_v2' +'/'

print(FIG_output_stem)
print(DATA_output_stem)

sc.settings.figdir=FIG_output_stem
#
adata=sc.read(DATA_output_stem+'adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm_generanked_DC_signatures_added.h5ad')
adata
#
#########
mg = [ 'S100A8','S100A9','VCAN','CD14','SERPINB2',
    'FCGR3A','CDKN1C','POU2F2','ZNF703','LYPD2', 
    'FABP4','INHBA','PPARG', 'ALDH2', 'TREM1',
    'LYZ','GRN','RETN','TNFSF13','PGD',
    'MCEMP1','S100A6','MARCO','TSPO','FBP1',
    'CCL13','LILRB5','FOLR2','SLC40A1','GPNMB',
    'CX3CR1','CCL3','CCL4','CCL4L2','KLF6',
    'SPP1','LGALS1','S100A10','BRI3','AQP9',
    'CLEC9A', 'CPVL','XCR1','BATF3','SNX3',
    'CLEC10A','FCGR2B','CD1C','PEA15','CD1E',
    'LAMP3', 'FSCN1', 'CCR7' ,'CCL19','BIRC3',
    'LILRA4','TCF4', 'IRF4','GZMB','IRF7',
    'CD207','CD1A','CD1E','LTB','FCGBP',
    'STMN1', 'MKI67', 'TOP2A', 'CDK1']  
avg_heatmap(adata, var_names = mg, groupby = "Defined_Cell_Subtype_function", ticksize = 5, figsize = (4,8),
            dendrogram = True,scale = True, use_raw = True, show_gene_labels = True, vmax = 3, vmin = -3, #cmap = "bwr", 
cmap='RdBu_r', save= FIG_output_stem+ '/heatmap_ave_test_markers_top_rank_genes.pdf')
    
    
    

#end






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

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

              
              
              