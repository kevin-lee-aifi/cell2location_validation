import scanpy as sc

def preprocess_vis_adata(adata):
    # Convert spatial coordinates to float
    adata.obsm['spatial'] = adata.obsm['spatial'].astype('float')
    
    # Ensure variable names are unique
    adata.var_names_make_unique()
    
    # Add a new column 'mt' in the var DataFrame to identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    # Calculate quality control metrics and add them to the .obs and .var fields
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    return adata
