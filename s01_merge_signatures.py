#!/usr/bin/env python
# coding: utf-8

# In[134]:


mycolor=["#7F3C8DFF",'#11A579FF','#3969ACFF','#F2B701FF','#E73F74FF','#80BA5AFF','#E68310FF','#008695FF', 
             '#CF1C90FF','#F97B72FF','#4B4B8FFF','#A5AA99FF','#B2DF8AFF','#33A02CFF','#FB9A99FF','#E31A1CFF',
             '#FDBF6FFF','#FF7F00FF','#CAB2D6FF','#6A3D9AFF','#FFFF99FF','#B15928FF','#A6CEE3FF','#1F78B4FF']

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


# In[135]:


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
sc.settings.set_figure_params(dpi=80, frameon=True, figsize=(3, 4), facecolor='white')


# In[136]:


FIG_output_stem = "scanpy_out/"

print(FIG_output_stem)


# CREATE FIGURE DIRECTORY IF IT DOES NOT EXIST     
d = os.path.dirname(FIG_output_stem)
if not os.path.exists(d):
        os.makedirs(d)
#
sc.settings.figdir=FIG_output_stem


# In[16]:


##-----------------------------------------------------------------------
## nsclc inhouse
##-----------------------------------------------------------------------
adata= sc.read('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
##
adata_lc1= adata
##
obs_lc1= adata.obs[['Defined_Cell_Subtype_function','cytotoxic_score', 'exhausted_score','age','Sample']]
obs_lc1['study']= 'nsclc_inhouse_kim_etal'
obs_lc1.rename(columns = {'Defined_Cell_Subtype_function':'celltype',}, inplace = True)
obs_lc1.head()


# In[17]:


adata_lc1


# In[20]:


adata_lc1.obs['Defined_Cell_Subtype_function'].value_counts()


# In[22]:


####annotate clusters
####
my_anno = {
    'CD4T':['C1-CD4-CCR7', 'C2-CD4-HSPA1A','C5-CD4-FOXP3','C3-CD4-CXCR6', 'C4-CD4-CXCL13', 'C6-CD4-RORC',
           
           ],
    'CD8T':['C10-CD8-ZNF683', 'C8-CD8-GZMK','C9-CD8-CX3CR1', 'C11-CD8-LAYN', 'C7-CD8-SLC4A10', 
            'C15-CD4/CD8-ISG15','C12-CD8-TRDV1'],
    'NK':	['C13-NK-FCGR3A', 'C14-NK-XCL1' ],
    'cycling':['C16-CD4/CD8-MKI67']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata_lc1.obs['cell_type_defined_pool'] = (
    adata_lc1.obs["Defined_Cell_Subtype_function"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[26]:


with rc_context({'figure.figsize': (3, 5)}):
    sc.pl.umap(adata_lc1, color=['louvain','Sample','Defined_Cell_Subtype_function',
                             'Defined_Cell_Subtype','cell_type_defined_pool'],#legend_loc='on data',
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1,save= '_nsclc_inhouse_kim_annotation_Tcell.pdf'
          )


# In[31]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata_lc1, color=['cell_type_defined_pool'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1,save= '_nsclc_inhouse_kim_annotation_Tcell_pool.pdf'
          )


# In[ ]:





# In[32]:


##------------------------
## nsclc all data
##------------------------
ddi= '/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data/adata.harmony.overclustered.filtered.CelltypeDefined.h5ad'
adata= sc.read(ddi)


# In[33]:


adata.obs['Defined_Cell_Type'].value_counts()


# In[34]:


adata.obs['Defined_Cell_Type'].cat.reorder_categories(
    ['T cells', 'NK cells', 'Myeloid cells','B cells', 'Mast cells', 'Epithelial', 'Endothelial', 'Fibroblasts'
    ], inplace=True)


# In[36]:


sc.pl.umap(adata, color=['louvain','Defined_Cell_Type','Sample'],
           palette=mycolor,
           hspace=0.5,
           ncols=1,save='__nsclc_inhouse_kim_global.pdf'
          )


# In[ ]:





# In[10]:


#adata.obs['Defined_Cell_Type'].value_counts()
#adata.obs["Defined_Cell_Type"] = [str(x) for x in adata.obs["Defined_Cell_Type"]]
#adata.obs.loc[adata_lc1.obs_names, "Defined_Cell_Type"] = adata_lc1.obs["Defined_Cell_Subtype_function"]
#sc.pl.umap(adata, color=['Defined_Cell_Type'],palette=color_CLASS,hspace=0.5,ncols=1#,save='_louvain_Defined_Cell_Type_Sample.pdf')


# In[186]:


##-----------------------------------------------------------------------
## nsclc pubdata
##-----------------------------------------------------------------------
###--------------------------------------
### calculate signature score for T cells
###--------------------------------------
adata= sc.read('/data/Zhuxq/pub_sc_data/lung_Leader_etal_CancerCell_2021/Leader_et_al/mergedanalysisdata/adata_combined_normed_TNK_tumor.h5ad')
# calculate cyto_exh_nk signature score
pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/cyto_exh_nk_signatures.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)
#
adata_lc2= adata


# In[55]:


marker_genes = ['CD8A','CD4','IL7R','CCR7','NKG7','FCGR3A',
                'FOXP3','EPCAM','KIT','CD68','CD79A','CLDN5','DCN','STMN1']
#with rc_context({'figure.figsize': (4, 5)}):  
sc.pl.umap(adata_lc2, color=marker_genes,hspace=0.5,wspace=0.5,ncols=3 ,
           color_map= 'plasma')


# In[56]:


#with rc_context({'figure.figsize': (5, 5)}):
sc.pl.umap(adata_lc2, color=['louvain2','sub_lineage','sample_ID','patient_ID'],
               palette=color_CLASS,
               hspace=0.5,
               ncols=1
          )


# In[60]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_lc2, color=['louvain2'],legend_loc='on data',
           #palette= my_color2,wspace=0.5,
           hspace=0.5,ncols=2
          )


# In[61]:


dict_mainlineage = {
    'Bcell': ['MS4A1','CD79A', 'IGHM','TNFRSF17'],
    'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Mono':['CD14','FCGR3A'],
    'Myeloid': ['CD68', 'LYZ','AIF1', 'CD163'],
    'NK': [ 'NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
    'DC':['IDO1','CLEC10A','CD1C','CD1E'],
    'Tcell': ['CD3D', 'CD3E', 'CD3G', 'TRAC'],
    'CD4T':['CCR7','IL7R','CD4'],
    'CD8T':['CD8A','CD8B'],
    'cycling':['STMN1']
}

sc.tl.dendrogram(adata_lc2, groupby='louvain2')
sc.pl.dotplot(adata_lc2, dict_mainlineage, groupby='louvain2',
              swap_axes=False,color_map='Blues', figsize=(14,7),
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2)


# In[187]:


####annotate clusters
####
my_anno = {
    'CD4T':[12,3,13,6,9,0,2,1],
    'CD8T':['7,0','4,0','4,1',5],
    'Bcell':[8],
    'NK':	['7,1',10],
    'cycling':[11]
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata_lc2.obs['cell_type_defined'] = (
    adata_lc2.obs["louvain2"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[65]:


#with rc_context({'figure.figsize': (5, 5)}):
sc.pl.umap(adata_lc2, color=['louvain2','cell_type_defined','sample_ID','patient_ID'],#legend_loc='on data',
           palette= color_CLASS,
           #wspace=0.5,
           hspace=0.5,
           ncols=1
          )


# In[66]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_lc2, color=['louvain2'],legend_loc='on data',
           #palette= my_color2,wspace=0.5,
           hspace=0.5,ncols=2
          )


# In[130]:


obs_lc2= adata_lc2.obs[['cell_type_defined','cytotoxic_score', 'exhausted_score','Age','sample_ID']]
obs_lc2['study']= 'nsclc_leader_etal'
obs_lc2.rename(columns = {'cell_type_defined':'celltype','Age':'age','sample_ID':'Sample'}, inplace = True)
obs_lc2.head()


# In[188]:


adata_lc2_fil= adata_lc2[~adata_lc2.obs['louvain2'].isin(['8'])].copy()


# In[14]:


marker_genes = ['CD8A','CD4','IL7R','CCR7','NKG7','FCGR3A','FOXP3','STMN1']
with rc_context({'figure.figsize': (4, 5)}):  
    sc.pl.umap(adata_lc2_fil, color=marker_genes,
               hspace=0.5,wspace=0.5,ncols=3 ,
           color_map= 'plasma',save= '_nsclc_Leader_etal_marker_Tcell.pdf')


# In[37]:


sc.pl.umap(adata_lc2_fil, color=['louvain2','cell_type_defined','sub_lineage'],
           palette=mycolor,
           hspace=0.5,
           ncols=1,
           save='_nsclc_Leader_etal_annotation_Tcell.pdf'
          )


# In[189]:


save_file = 'scanpy_out/_adata_nsclc_Leader_etal_Tcell.h5ad'
adata_lc2_fil.write_h5ad(save_file)


# In[70]:


##---------------------
## Total cell
##---------------------
adata= sc.read('/data/Zhuxq/pub_sc_data/lung_Leader_etal_CancerCell_2021/Leader_et_al/mergedanalysisdata/adata_combined_normed.h5ad')
adata=adata[adata.obs['tissue'].isin(['Tumor'])].copy()


# In[71]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1)


# In[79]:


marker_genes = ['CD8A','CD4','IL7R','CCR7','NKG7','FCGR3A',
                'FOXP3','EPCAM','KIT','CD68','CD79A','CLDN5','DCN','STMN1']
#with rc_context({'figure.figsize': (4, 5)}):  
sc.pl.umap(adata, color=marker_genes,hspace=0.5,wspace=0.5,ncols=3 ,
           color_map= 'plasma')


# In[72]:


adata.obs['sub_lineage2']=  


# In[80]:


sc.pl.umap(adata, 
           color=['louvain','lineage','sub_lineage','norm_group','lig_rec_group'],
           palette=color_CLASS,
           hspace=0.5,
           ncols=1,
           save='_nsclc_Leader_etal_global.pdf'
          )


# In[81]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['louvain'],legend_loc='on data',
           #palette= my_color2
               wspace=0.5,
           hspace=0.5,ncols=1
          )


# In[74]:


save_file = 'scanpy_out/_adata_nsclc_Leader_etal_global.h5ad'
adata.write_h5ad(save_file)


# In[84]:


pd.crosstab(index= adata.obs['louvain'], columns= adata.obs['lig_rec_group'])


# In[41]:


####annotate clusters
####
my_anno = {
    'B':[1,17],	
    'DC':[12,21],
    'NK':[11],
    #'CD4T':[0,4,6,10],	
    #'CD8T':[7,8,20,22],
    'T':[0,4,6,10,7,8,20,22 ],
    'mac':[2,9],
    'mast':[13],
    'mono':[5],
    'plasma':[3],
    'cycling':[14],
    'Epithelial':[15,18,23],
    'Endo':[16],
    'Fib':[19]
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['cell_type_defined'] = (
    adata.obs["louvain"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[42]:


sc.pl.umap(adata, 
           color=['louvain','lineage','sub_lineage','norm_group','lig_rec_group','cell_type_defined'],
           palette=color_CLASS,
           hspace=0.5,
           ncols=1,
           save='_nsclc_Leader_etal_global.pdf'
          )


# In[43]:


adata.obs['cell_type_defined'].value_counts()


# In[45]:


adata.obs['cell_type_defined'].cat.reorder_categories(
    ['T', 'NK', 'mac','B', 'mast', 'Epithelial', 'Endo', 'Fib','plasma', 'mono','DC','cycling'
    ], inplace=True)


# In[48]:


with rc_context({'figure.figsize': (4, 4)}): 
    sc.pl.umap(adata, 
           color=['cell_type_defined'],
           palette=mycolor,
           hspace=0.5,
           ncols=1,
           save='_nsclc_Leader_etal_global_color.pdf'
          )


# In[49]:


save_file = 'scanpy_out/_adata_nsclc_Leader_etal_global.h5ad'
adata.write_h5ad(save_file)


# In[ ]:





# In[50]:


##-----------------------------------------------------------------------
## CRC pubdata
##-----------------------------------------------------------------------
adata_crc= sc.read('/data/Zhuxq/pub_sc_data/CRC_GSE178341_Cell2021/Totaldata/adata_total_doubleted_qc_doublet_filed_only_tumor_tissue_age_grouped_my_defined_cell_type2.h5ad')


# In[51]:


adata_crc


# In[52]:


adata_crc.obs['cell_type2'].cat.reorder_categories(
    ['Tcell', 'NK', 'Myeloid','Bcell', 'Mast', 'Epithelial', 'Endothelial', 'Fibroblasts','Plasma','cycling'
    ], inplace=True)


# In[100]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_crc, color=['louvain','cell_type2','PID'],#legend_loc='on data',
           #palette= my_color2,
               wspace=0.5,
           hspace=0.5,ncols=2
          )


# In[53]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata_crc, color=['cell_type2'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1, save= '_crc_global_color.pdf'
          )


# In[42]:


pd.crosstab(index= adata_crc.obs['louvain'], columns= adata_crc.obs['cell_type2'])


# In[61]:


####annotate clusters
####
my_anno = {
    'CD4T':[4,9],
    'CD8T':[2],
    'Bcell':[8],
    'Endothelial':[14],
    'Epithelial':[0,3,5,12,13,19,21]	,
    'Fibroblasts':	[15],
    'Mast':	[18],
    'Myeloid':	[1,7,17,20],
    'NK':	[10],
    'Plasma':	[6,11],
    'cycling':[16]
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata_crc.obs['cell_type_defined'] = (
    adata_crc.obs["louvain"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[45]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_crc, color=['louvain','cell_type_defined'],#legend_loc='on data',
           #palette= my_color2,
               wspace=0.5,
           hspace=0.5,ncols=2
          )


# In[101]:


obs_crc= adata_crc.obs[['cell_type_defined','cytotoxic_score (T cell)', 'exhausted_score (T cell)','Age','PID']]
obs_crc['study']= 'CRC_GSE178341'
obs_crc.rename(columns = {'cell_type_defined':'celltype',
                          'cytotoxic_score (T cell)':'cytotoxic_score'	,
                          'exhausted_score (T cell)':'exhausted_score',
                          'Age':'age',
                         'PID':'Sample'
                         
                         }, inplace = True)
obs_crc.head()


# ###---------------------
# ### reannotate TNK cell for CRC
# ###---------------------

# In[72]:


adata= adata_crc[adata_crc.obs['louvain'].isin(['2','4','9','10'])].copy()


# In[73]:


adata


# In[74]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'PatientTypeID', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[75]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)


# In[76]:


sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype')  # to check the doublet again, can set resolution bigger enough


# In[91]:


save_file = 'scanpy_out/_adata_CRC_GSE178341_Tcell.h5ad'
adata.write_h5ad(save_file)


# In[180]:


#adata=sc.read('scanpy_out/_adata_CRC_GSE178341_Tcell.h5ad')


# In[79]:


#with rc_context({'figure.figsize': (4, 4)}):
sc.pl.umap(adata, color=['louvain_subtype','cell_type_defined','clTopLevel',
                         'clMidwayPr','cl295v11SubShort','cl295v11SubFull'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1#, save= '_crc_global_color.pdf'
          )


# In[81]:


pd.crosstab(index= adata.obs['cl295v11SubFull'], columns = adata.obs['louvain_subtype'])


# In[181]:


adata_t= adata[adata.obs['cl295v11SubFull'].str.contains('cTN')].copy()


# In[94]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata_t, color=['cell_type_defined'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1, save= '_CRC_Tcell_color.pdf'
          )
##
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_t, color=['cl295v11SubFull'],#legend_loc='on data',
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1, save= '_CRC_Tcell_subtype_color.pdf'
          )


# In[182]:


save_file = 'scanpy_out/_adata_CRC_GSE178341_Tcell_clean.h5ad'
adata_t.write_h5ad(save_file)


# In[ ]:





# In[ ]:





# In[95]:


##-----------------------------------------------------------------------
## GC pubdata
##-----------------------------------------------------------------------
# total
#
adata= sc.read('/data/Zhuxq/pub_sc_data/GC_GSE183904_CD2021/copy_node4/Totaldata/adata_GC_GSE183904_26samples_PT_doubleted_qc_filtered_clustered_defined_sigscoreCaled.h5ad')


# In[96]:


adata


# In[99]:


adata.obs['cell_type'].cat.reorder_categories(
    ['Tcell', 'NK', 'Myeloid','Bcell', 'Mast', 'Epithelial', 'Endothelial', 'Fibroblasts','Plasma'
    ], inplace=True)


# In[101]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cell_type'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1, 
               save= '_GC_global_color.pdf'
          )


# In[176]:


## Tcell
adata_gc= sc.read('/data/Zhuxq/pub_sc_data/GC_GSE183904_CD2021/copy_node4/TNKdata/_adata_cluster_doubletADD_louvain_inital_annotated_refined_age_group2.h5ad')


# In[177]:


adata_gc= adata_gc[~adata_gc.obs['cell_type'].isin(['TNK_others'])].copy()


# In[178]:


adata_gc.obs['cell_type'].value_counts()


# In[111]:


####annotate clusters
####
my_anno = {
    'CD4T':['CD4-FOXP3','CD4-CCR7', 'CD4-ANXA1', 'CD4-PRF1', 'CD4-CXCL13' ],
    'CD8T':['CD8-GZMK', 'CD8-FTH1', 'CD8-ZNF683','CD8-GZMH','CD8-CD52',
            'CD8-CXCL13','CD8-FOS','CD8-STMN1','CD8-ISG15'],
    'NK':	['NK']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata_gc.obs['cell_type_pool'] = (
    adata_gc.obs["cell_type"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[112]:


#with rc_context({'figure.figsize': (5, 5)}):
sc.pl.umap(adata_gc, color=['louvain','cell_type','cell_type_pool','sample'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
          )


# In[114]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata_gc, color=['cell_type_pool'],#legend_loc='on data',
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1, save= '_GC_Tcell_color.pdf'
          )
##
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata_gc, color=['cell_type'],#legend_loc='on data',
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1, save= '_GC_Tcell_subtype_color.pdf'
          )


# In[179]:


save_file = 'scanpy_out/_adata_GC_GSE183904_Tcell.h5ad'
adata_gc.write_h5ad(save_file)


# In[113]:


obs_gc= adata_gc.obs[['cell_type','cytotoxic_score (T cell)', 'exhausted_score (T cell)','Age','sample']]
obs_gc['study']= 'GC_GSE183904'
obs_gc.rename(columns = {'cell_type':'celltype',
                          'cytotoxic_score (T cell)':'cytotoxic_score'	,
                          'exhausted_score (T cell)':'exhausted_score',
                          'Age':'age',
                        'sample':'Sample'
                        }, inplace = True)
obs_gc.head()


# In[ ]:





# In[ ]:





# In[160]:


##-----------------------------------------------------------------------
## HCC pubdata
##-----------------------------------------------------------------------
adata_hcc=  sc.read('/data/Zhuxq/pub_sc_data/liver_Zhangzemin_etal_GSE140228_Cell_2019/out/Totalonly2019data/adata_qcfil_normed_signatured_renamed.h5ad')


# In[161]:


adata_hcc= adata_hcc[~adata_hcc.obs['celltype_global'].isin(['Myeloid-Liver-doublets'])].copy()
adata_hcc= adata_hcc[adata_hcc.obs['celltype_global'].notna()].copy()


# In[139]:


#with rc_context({'figure.figsize': (4, 5)}):
sc.pl.umap(adata_hcc, color=['louvain','celltype_sub','celltype_global','Sample','Donor'],
           #palette= my_color2,
               wspace=0.5,
           hspace=0.5,ncols=1#,save='_test_louvain.pdf'
          )
#with rc_context({'figure.figsize': (5, 5)}):
sc.pl.umap(adata_hcc, color=['louvain'],legend_loc='on data',
           #palette= my_color2,
               wspace=0.5,
           hspace=0.5,ncols=1
          )


# In[145]:


adata_hcc.obs['celltype_global'].value_counts()


# In[147]:


adata_hcc.obs['celltype_global'].cat.reorder_categories(
    ['Lymphoid-T', 'Lymphoid-NK', 'Myeloid','Lymphoid-B', 
     'Myeloid-Mast', 'Lymphoid-B-Plasma','Lymphoid-T-NK-Cycling','ILCs'
    ], inplace=True)


# In[149]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata_hcc, color=['celltype_global'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1,save='_HCC_global_color.pdf'
          )


# In[140]:


pd.crosstab(index= adata_hcc.obs['celltype_sub'], columns= adata_hcc.obs['louvain'])


# In[141]:


marker_genes = ['CD8A','CD4','IL7R','CCR7','NKG7','FCGR3A',
                'FOXP3','EPCAM','KIT','CD68','CD79A','CLDN5','DCN','STMN1']
with rc_context({'figure.figsize': (4, 5)}):  
    sc.pl.umap(adata_hcc, color=marker_genes,hspace=0.5,wspace=0.5,ncols=3 ,
           color_map= 'plasma',save='_markers_HCC_global.pdf')


# In[142]:


dict_mainlineage = {
    'Bcell': ['MS4A1','CD79A', 'IGHM','TNFRSF17'],
    'Mast': ['KIT', 'MS4A2','GATA2','TPSB2'],
    'Mono':['CD14','FCGR3A'],
    'Myeloid': ['CD68', 'LYZ','AIF1', 'CD163'],
    'NK': [ 'NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
    'DC':['IDO1','CLEC10A','CD1C','CD1E'],
    'Tcell': ['CD3D', 'CD3E', 'CD3G', 'TRAC'],
    'CD4T':['CCR7','IL7R','CD4'],
    'CD8T':['CD8A','CD8B'],
    'cycling':['STMN1']
}

sc.tl.dendrogram(adata_hcc, groupby='louvain')
sc.pl.dotplot(adata_hcc, dict_mainlineage, groupby='louvain',
              swap_axes=False,color_map='Blues', figsize=(14,7),
              dendrogram=True, use_raw=True,standard_scale='var',
              dot_max=0.5,dot_min=0.2)


# In[162]:


#####------
##### only for T cells
#####------
adata_hcc=  adata_hcc[adata_hcc.obs['celltype_sub'].str.startswith(('CD','NK'))].copy()
adata= adata_hcc


# In[163]:


adata_hcc


# In[164]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)g
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'Sample', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[165]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)


# In[166]:


sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype')  # to check the doublet again, can set resolution bigger enough


# In[168]:


adata.obs['celltype_sub'].value_counts()


# In[170]:


my_anno = {
    'NK':	['NK-C7-CD160','NK-C5-CD69','NK-C4-IFNG-HSPA1A','NK-C8-CD160-HSPA1A','NK-C3-IFNG','NK-C6-IL7R',
            'NK-C1-FCGR3A','NK-C9-MKI67','NK-C2-SELL'
            ]
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['celltype_sub2'] = (
    adata.obs["celltype_sub"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[172]:


adata.obs['celltype_sub3'] =  adata.obs['celltype_sub'].str.slice(0, 3)


# In[175]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['celltype_sub'],
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1,save='_HCC_Tcell_subtype_color2.pdf'
          )


# In[173]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['celltype_sub3'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1,save='_HCC_Tcell_pool_color.pdf'
          )

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['celltype_sub2'],
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1,save='_HCC_Tcell_subtype_color.pdf'
          )


# In[174]:


save_file = 'scanpy_out/_adata_HCC_GSE140228_Tcell.h5ad'
adata.write_h5ad(save_file)


# In[ ]:





# In[ ]:





# In[70]:


my_anno = {
    'CD4T':[5,12,0],
    'CD8T':[1,7,16],
    'Bcell':[8,10],
    'Mast':	[17],
    'Myeloid':	[2,6,9,13,14,15],
    'NK':	[3,4],
    'cycling':[11]
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata_hcc.obs['cell_type_defined'] = (
    adata_hcc.obs["louvain"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[71]:


with rc_context({'figure.figsize': (4, 5)}):
    sc.pl.umap(adata_hcc, color=['louvain','cell_type_defined'],
           #palette= my_color2,
               wspace=0.5,
           hspace=0.5,ncols=2#,save='_test_louvain.pdf'
          )


# In[72]:


adata_hcc


# In[114]:


obs_hcc= adata_hcc.obs[['cell_type_defined','cytotoxic_score (T cell)', 'exhausted_score (T cell)','Age','Donor']]
obs_hcc['study']= 'HCC_GSE140228'
obs_hcc.rename(columns = {'cell_type_defined':'celltype',
                          'cytotoxic_score (T cell)':'cytotoxic_score'	,
                          'exhausted_score (T cell)':'exhausted_score',
                          'Age':'age',
                         'Donor':'Sample'
                         }, inplace = True)
obs_hcc.head()


# In[131]:


obs_merge= pd.concat([obs_lc1, obs_lc2, obs_crc, obs_hcc, obs_gc], axis=0)


# In[132]:


obs_merge


# In[133]:


obs_merge.to_csv('merged_obs_Tcell_sigscore_5_cohorts.csv')


# ###---------
# ### PD1 scRNA-seq
# ###---------

# In[272]:


adata= sc.read('/data/Zhuxq/pub_sc_data/ccRCC_Krishna_etal_CancerCell_2021/out/Totaldata/adata_tumourTissue_cluster_name2.h5ad')


# In[273]:


adata.obsm['X_umap']=adata.obs[['UMAP1','UMAP2']].to_numpy()


# In[274]:


adata= adata[~adata.obs['cluster_name2'].isin(['Ambiguous'])].copy()


# In[275]:


adata


# In[276]:


adata.obs['cluster_name2'].value_counts()


# In[277]:


adata.obs['cluster_name2'].cat.reorder_categories(
    ['T cell', 'NK cell', 'TAM-HLA+','B cell', 'Mast', 'Epithelial', 'Endothelial', 'Fibroblast','TAM-ISG+',
     'DCs','pDC','Monocyte','Megakaryocyte'
    ], inplace=True)


# In[206]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name2'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_global_color.pdf'
          )


# In[290]:


adata.obs['treatment']= adata.obs['treatment'].str.replace(' ','_')


# In[291]:


adata.obs['treatment'].value_counts()


# In[292]:


my_anno = {
    'Response':['_Ipi.Nivo_mixed_response', '_Ipi.Nivo_complete_response'],
    'Non-response':['_Ipi.Nivo_resistant'
           ],
    'Exposed':['_Nivo.exposed'],
    'Treatment.naive':['_Untreated_1','_Untreated_2']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['treatmentGroup'] = (
    adata.obs["treatment"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[293]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['treatmentGroup'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
              # ,save='_ccRCC_global_color.pdf'
          )


# In[294]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata[adata.obs['treatmentGroup'].isin(['Response','Non-response'])], 
               color=['treatmentGroup'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
              # ,save='_ccRCC_global_color.pdf'
          )


# In[296]:


adata=adata[adata.obs['treatmentGroup'].isin(['Response','Non-response'])].copy()


# In[297]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'Sample_name', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[298]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)


# In[ ]:


sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[316]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name2'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_global_color_response_celltype.pdf'
          )


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, 
               color=['treatmentGroup'],
           palette= ["#E69F00","#56B4E9"],
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_global_color_response.pdf'
          )


# In[301]:


save_file = 'scanpy_out/_adata_ccRCC_global_responsesample.h5ad'
adata.write_h5ad(save_file)


# In[ ]:





# In[ ]:





# In[ ]:





# In[216]:


####
### focus on T cells
####
adata=  adata[adata.obs['cluster_name'].str.startswith(('CD4+','CD8A+','NK'))].copy()


# In[ ]:





# In[223]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'Sample_name', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[224]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)


# In[225]:


sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[ ]:





# In[222]:


adata.obs['Sample_name'].value_counts()


# In[218]:


adata


# In[228]:


adata.obs['cluster_name'].value_counts()


# In[231]:


my_anno = {
    'CD4T':['CD4+ Activated IEG', 'CD4+ Treg','CD4+ Effector', 'CD4+ Proliferating', 'CD4+ Naive' ],
    'CD8T':['CD8A+ Exhausted', 'CD8A+ Tissue-resident', 'CD8A+ Proliferating','CD8A+ NK-like',
            'CD8A+ Exhausted IEG',
           ],
    'NK':['NK HSP+']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['cluster_name3'] = (
    adata.obs["cluster_name"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[233]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name3'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_Tcell_pool_color.pdf'
          )
    
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name'],
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_Tcell_subtype_color.pdf'
          ) 
    
    


# In[234]:


save_file = 'scanpy_out/_adata_ccRCC_Tcell.h5ad'
adata.write_h5ad(save_file)


# In[ ]:


##----------
## Tcell response samples
##--------------


# In[302]:


adata=  adata[adata.obs['cluster_name'].str.startswith(('CD4+','CD8A+','NK'))].copy()
adata.obs['treatment']= adata.obs['treatment'].str.replace(' ','_')

my_anno = {
    'Response':['_Ipi.Nivo_mixed_response', '_Ipi.Nivo_complete_response'],
    'Non-response':['_Ipi.Nivo_resistant'
           ],
    'Exposed':['_Nivo.exposed'],
    'Treatment.naive':['_Untreated_1','_Untreated_2']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['treatmentGroup'] = (
    adata.obs["treatment"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)



# In[303]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'Sample_name', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[304]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)


# In[305]:


sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[307]:


my_anno = {
    'CD4T':['CD4+ Activated IEG', 'CD4+ Treg','CD4+ Effector', 'CD4+ Proliferating', 'CD4+ Naive' ],
    'CD8T':['CD8A+ Exhausted', 'CD8A+ Tissue-resident', 'CD8A+ Proliferating','CD8A+ NK-like',
            'CD8A+ Exhausted IEG',
           ],
    'NK':['NK HSP+']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['cluster_name3'] = (
    adata.obs["cluster_name"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[310]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name3'],
           palette= mycolor,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_Tcell_pool_color_response.pdf'
          )
    
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['cluster_name'],
           palette= color_CLASS,
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_Tcell_subtype_color_response.pdf'
          ) 
    
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['treatmentGroup'],
           palette= ["#E69F00","#56B4E9"],
               wspace=0.5,
           hspace=0.5,ncols=1
               ,save='_ccRCC_Tcell_subtype_color_responsetreat.pdf'
          ) 


# In[ ]:





# In[ ]:





# In[348]:


###------------
### BCC scRNA pd1
###------------
adata= sc.read('/data/Zhuxq/pub_sc_data/melanoma_Yost_etal_GSE123813_NM_2019/out/Totaldata/adataall_bcc_renamed.h5ad')


# In[349]:


adata


# In[350]:


adata.obsm['X_umap']=adata.obs[['UMAP1','UMAP2']].to_numpy()

adata.obs['cluster2'].cat.reorder_categories(
    [
'T cell','NK_cells', 'Macrophages','B cell', 'pDCs','Tumor', 
        'Endothelial','Fibroblasts', 'Plasma_cells','DCs','Melanocytes'
    ], inplace=True)


# In[351]:


adata.obs['treatment'].value_counts()


# In[352]:


pd.crosstab(index= adata.obs['treatment'], columns=  adata.obs['Response'])


# In[245]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster2'#,'cluster'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_global_cell_color.pdf'
          )


# In[246]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           #save='_BCC_global_cell_color.pdf'
          )


# In[271]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['Response'
                        ],
           palette= ["#E69F00","#56B4E9"],wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_global_cell_color_response.pdf'
          )


# In[344]:


adata


# In[327]:


adata.obs['mygroup']=adata.obs['treatment'].astype(str) + adata.obs['Response'].astype(str)


# In[329]:


# calculate cyto_exh_nk signature score
pws = pd.read_table('/data/Zhuxq/young_LC_analysis/gseaDB/cyto_exh_nk_signatures.txt', delimiter='\t')
col=pws.columns.values
for i in col: 
    tmp = pws[i].dropna()
    sc.tl.score_genes(adata,tmp,score_name=i,use_raw=True)


# In[330]:


adata


# In[331]:


with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, ['cytotoxic_score', 'exhausted_score'], 
                 groupby='mygroup', stripplot=False, inner='box')  # use stripplot=False to remove the internal dots, inner='box' adds a boxplot inside violins


# In[345]:


###-------------
### only use pretreated samples
###-----------------
adata= adata[adata.obs['treatment'].isin(['pre'])].copy()

# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[346]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster2'#,'cluster'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
          save='_BCC_global_cell_color_response_pre.pdf'
          )
    
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['Response'
                        ],
           palette= ["#E69F00","#56B4E9"],wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
          save='_BCC_global_cell_color_responsegroup_pre.pdf'
          )


# In[347]:


save_file = 'scanpy_out/_adata_BCC_global_cell.pre.h5ad'
adata.write_h5ad(save_file)


# In[ ]:





# In[252]:


###-----
###focus on T cell 
###-----
adata=  adata[adata.obs['cluster'].str.startswith(('CD','Tcell','Treg','NK'))].copy()


# In[253]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           #save='_BCC_global_cell_color.pdf'
          )


# In[254]:


adata


# In[255]:


# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T


# In[256]:


sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[258]:


adata.obs['cluster'].value_counts()


# In[260]:


my_anno = {
    'CD4T':['CD4_T_cells','Tregs'],
    'CD8T':['CD8_mem_T_cells','CD8_act_T_cells', 'CD8_ex_T_cells'
           ],
    'NK':['NK_cells'],
    'cycling':['Tcell_prolif']
}
annot_dict = {
    str(c): ct for ct, clusters in my_anno.items() for c in clusters
}
annot_dict

adata.obs['cluster2'] = (
    adata.obs["cluster"]
    .map(lambda x: annot_dict.get(x, x))
    .astype("category")
)


# In[261]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster2'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_pool_color.pdf'
          )
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster'
                        ],
           palette= color_CLASS,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_color.pdf'
          )
    
    


# In[262]:


save_file = 'scanpy_out/_adata_BCC_Tcell.h5ad'
adata.write_h5ad(save_file)


# In[263]:


adata


# In[266]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['Response'
                        ],
           palette= ["#E69F00","#56B4E9"],wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_color_response.pdf'
          )


# In[336]:


###-----
### only pre T cells
###-----

adata=sc.read('scanpy_out/_adata_BCC_Tcell.h5ad')


###-------------
### only use pretreated samples
###-----------------
adata= adata[adata.obs['treatment'].isin(['pre'])].copy()

# start to normalize and cluster
#
####normalzation
# use raw counts to reanalyze
adata.X= adata.layers["counts"]
# normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.highly_variable_genes(adata, n_top_genes = 4000)
#adata= adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt','CC.Difference'])
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
#sce.pp.harmony_integrate(adata, 'sample')
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
#
import harmonypy as hm
ho=hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient', max_iter_harmony=50)
#adjusted_pcs = pd.DataFrame(ho.Z_corr)
#adata.obsm['X_pca'] = adjusted_pcs
#adata.obsm['X_pca_harmonized'] = adjusted_pcs
adata.obsm['X_pca_harmonized'] =ho.Z_corr.T

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50,use_rep='X_pca_harmonized')
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1, key_added='louvain_subtype') 


# In[338]:


with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster2'
                        ],
           palette= mycolor,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_pool_color_pre.pdf'
          )
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['cluster'
                        ],
           palette= color_CLASS,wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_color_pre.pdf'
          )

with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color=['Response'
                        ],
           palette= ["#E69F00","#56B4E9"],wspace=0.5,size= 5,title= 'BCC_GSE123813',
           hspace=0.5,ncols=1,
           save='_BCC_T_cell_color_response_pre.pdf'
          )


# In[339]:


save_file = 'scanpy_out/_adata_BCC_Tcell_pre.h5ad'
adata.write_h5ad(save_file)


# In[ ]:




