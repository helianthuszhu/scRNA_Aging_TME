import scanpy as sc
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import os
from matplotlib.pyplot import rc_context
import re
#
adataall= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data/adata.harmony.overclustered.filtered.CelltypeDefined.withrawcounts.h5ad')
adataall.obs.Defined_Cell_Type.value_counts()
#
adata_tcell= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad')
#
adata_myecell= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Myeloidcell_v2/adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined.h5ad')
#
adata_fib= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Fibroblasts_v2/adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined.h5ad')
#
adata_endo= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Endothelial_v2/adata_Endothelial_v2_degs_subtype_Defined.h5ad')
adata_bcell= sc.read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Bcell_v2/adata_Bcell_v2_degs_subtype_Defined.h5ad')
#
adata_mast= adataall[adataall.obs['Defined_Cell_Type'].isin(['Mast cells'])]
adata_epi= adataall[adataall.obs['Defined_Cell_Type'].isin(['Epithelial'])]

#####
obs_tcell=adata_tcell.obs[['Index','Defined_Cell_Subtype_function']]
obs_myecell= adata_myecell.obs[['Index','Defined_Cell_Subtype_function']]
obs_fib= adata_fib.obs[['Index','Defined_Cell_Subtype_function']]
obs_endo= adata_endo.obs[['Index','Defined_Cell_Subtype_function']]
obs_bcell=adata_bcell.obs[['Index','Defined_Cell_Subtype_function']]
obs_mast= adata_mast.obs[['Index','Defined_Cell_Type']]
obs_mast = obs_mast.rename(columns={'Defined_Cell_Type': 'Defined_Cell_Subtype_function'})
obs_epi= adata_epi.obs[['Index','Defined_Cell_Type']]
obs_epi = obs_epi.rename(columns={'Defined_Cell_Type': 'Defined_Cell_Subtype_function'})
#####
frames = [obs_tcell, obs_myecell, obs_fib, obs_endo, obs_bcell, obs_mast, obs_epi]
obs_merged = pd.concat(frames)
obs_merged.to_csv('output6_merge_adataobs/obs_merged.csv')
#####
#####


adata_sel= adataall[adataall.obs.index.isin(obs_merged.index)]
adata_sel.obs['Defined_Cell_subType']= obs_merged['Defined_Cell_Subtype_function']
adata_sel.obs['Defined_Cell_subType'].value_counts()

save_file = 'output6_merge_adataobs/'+'adata_defined_58031cells.h5ad'
adata_sel.write_h5ad(save_file)