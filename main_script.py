import os

import numpy as np
import scanpy as sc
import scvelo as scv
import detecting_bifurcations_in_scRNA_seq_data as dbsc
from multiprocessing import freeze_support


def process_and_save_data(filename):
    """
        Pre-process the data from a h5ad file using standard Scanpy procedure and saves the pre-processed
        data to a new h5ad file.

        Parameters:
            filename (str): The base filename (without extension) of the input h5ad file.

        Returns:
            None
    """
    # Check if the input file exists with .h5ad extension
    filename_with_extension = f'{filename}.h5ad'
    if not os.path.isfile(filename_with_extension):
        print(f"The input file '{filename_with_extension}' does not exist.")
        return

    # Load or create the processed file name with .h5ad extension
    processed_filename = f'processed_{filename}.h5ad'

    # Check if the processed file already exists
    if os.path.isfile(processed_filename):
        print(f"The processed file '{processed_filename}' already exists.")
        return

    # Continue with processing the data
    adata = sc.read_h5ad(filename_with_extension)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 3500, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]

    sc.pp.normalize_total(adata, target_sum=1e4)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    adata.uns['iroot'] = np.flatnonzero(adata.obs['time'] == '0d')[0]
    sc.tl.dpt(adata)

    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)

    scv.pp.filter_and_normalize(adata, n_top_genes=50)

    # Save the processed file with .h5ad extension
    adata.write_h5ad(processed_filename)


def main():
    freeze_support()
    filename1 = 'OVCA420_TGFB1'
    # Move the below function as a main script for tutorial.
    # process_and_save_data(filename1)

    #  adata = sc.read_h5ad(f'processed_{filename1}.h5ad')
    #  dbsc.tl.traj_inference(adata) #run before plot_stability and spectrum_full
    #  dbsc.pl.plot_stability(adata, filename1) #creates regression file (pos and max eig values)
    #  dbsc.pl.spectrum_full(adata, filename1) #creates spectrum file (eigen spectrum)
    #  dbsc.tl.create_weights_geneexpress(adata, filename1)

    # Below for spliceJAC parameter variation
    #  dbsc.tl.par_width_var(adata)
    #  dbsc.pl.plot_traj_widthvar(adata, filename1)

    # CD adata file generation
    #  adata = sc.read_h5ad(f'processed_{filename1}.h5ad')
    #  dbsc.tl.create_weights_geneexpress(adata, filename1)  # create adata with gene exp in windows

    # Community detection
    #  adata = sc.read_h5ad(f'gene_exp_{filename1}.h5ad')
    #  jac_list = adata.uns['Jacobian']['jacobians']  # stores range of Jacobians calculated
    #  geneexp_list = adata.uns['Cells_Captured']['gene_expression']  # stores range of mean expressions calculated
    #  G_list=dbsc.tl.G_listgen_a(adata, jac_list)
    #  Community_list1, Num_com_list1=dbsc.tl.cd_g_w(G_list)
    #  dbsc.pl.plot_comm(adata, Num_com_list1, filename1)

    #  G_list_2 = dbsc.tl.G_listgen_geneexp(adata, jac_list, geneexp_list)
    #  Community_list2, Num_com_list2=dbsc.tl.cd_grM_w(G_list)
    #  dbsc.pl.plot_comm(adata, Num_com_list2, filename1, method_key='greedy_modularity')


if __name__ == '__main__':
    main()
