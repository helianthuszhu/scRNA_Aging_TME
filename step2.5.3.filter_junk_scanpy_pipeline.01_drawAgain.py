####
####
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import scanpy.external as sce

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

sc.settings.set_figure_params(dpi=300, frameon=True, figsize=(3, 4), facecolor='white')

FIG_output_stem = "./output2_step2.5.3.filter_junk_scanpy_pipeline/" + 'figure_draw' +'/'
d = os.path.dirname(FIG_output_stem)
if not os.path.exists(d):
        os.makedirs(d)
        

sc.settings.figdir=FIG_output_stem

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
adata=sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data/adata.harmony.overclustered.filtered.CelltypeDefined.h5ad')
#
sc.pl.umap(adata, color=['louvain','Defined_Cell_Type','Sample'],palette=color_CLASS,hspace=0.5,ncols=1,save='_louvain_Defined_Cell_Type_Sample.pdf')

marker_genes = ['EPCAM','CD3D','CD79A','NKG7','CD68','KIT','DCN','RAMP2','MKI67']
sc.pl.umap(adata, color=marker_genes,hspace=0.5,ncols=1 ,color_map= 'plasma',save='_marker_individuals.mainlineage.pdf')

#
#define age group
cut_labels_3 = ['Young', 'Intermediated', 'Aged']
cut_bins = [20, 49, 60,100]
adata.obs['intage'] = adata.obs['age'].astype(int)
adata.obs['age_group'] = pd.cut(adata.obs['intage'], bins=cut_bins, labels=cut_labels_3)
#

sc.tl.embedding_density(adata, basis='umap', groupby='age_group')

sc.pl.embedding_density(adata, basis='umap', key='umap_density_age_group', group='all',hspace=0.5,ncols=1 ,save='_age_density.pdf')

#sc.pl.umap(adata, color=['intage', 'nCount_RNA', 'nFeature_RNA', 'percent_mt'],color_map= 'viridis_r', hspace=0.5,ncols=1,save='_intage.qc.pdf')
sc.pl.umap(adata, color=['intage', 'nCount_RNA', 'nFeature_RNA', 'percent_mt'],color_map= 'plasma', hspace=0.5,ncols=1,save='_intage.qc.pdf')

#
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
#
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
sc.pl.matrixplot(adata, marker_genes_dict, 'Defined_Cell_Type', dendrogram=True,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')

sc.pl.matrixplot(adata, var_names=marker_genes_dict, groupby='Defined_Cell_Type', dendrogram=False, cmap='plasma',swap_axes=True,
 standard_scale='var', colorbar_title='column scaled\nexpression',save='_mainleage.pdf')

#
#
#########
import pandas as pd
from plotnine import *
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

base_plot= ggplot(adata.obs, aes(x='age_group', fill='Defined_Cell_Type')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+
        theme_classic()+theme(figure_size=(3, 5))
]
save_as_pdf_pages(plots,filename='stacked_priportion_ageGroup_vs_celltype.pdf',path= FIG_output_stem)

#
base_plot= ggplot(adata.obs, aes(x='Sample', fill='Defined_Cell_Type')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(2, 5))
]
save_as_pdf_pages(plots,filename='stacked_priportion_Sample_vs_celltype.pdf',path= FIG_output_stem)

#
base_plot= ggplot(adata.obs, aes(x='Sample', fill='Defined_Cell_Type')) 
plots=[base_plot+
geom_bar()+scale_fill_manual(values = color_CLASS)+ylab('percentage')+coord_flip()+
        theme_classic()+theme(figure_size=(2, 5))
]
save_as_pdf_pages(plots,filename='stacked_abscounts_Sample_vs_celltype.pdf',path= FIG_output_stem)


base_plot= ggplot(adata.obs, aes(x='age_group', fill='Phase')) 
plots=[base_plot+
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+ylab('percentage')+
        theme_classic()+theme(figure_size=(3, 6))
]
save_as_pdf_pages(plots,filename='stacked_priportion_ageGroup_vs_phase.pdf',path= FIG_output_stem)

#
adata.obs.to_csv('output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/stacked_plot_ordered/adata.obs.for.stacked.plot.csv')