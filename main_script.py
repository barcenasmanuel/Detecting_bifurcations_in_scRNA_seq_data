import scanpy as sc
import Generalized_Method_W24 as gm
from multiprocessing import freeze_support


def main():
    freeze_support()
    filename1 = 'OVCA420_TGFB1'
    # gm.tl.process_and_save_data(filename1) #adapted load_preprocess function

    adata = sc.read_h5ad(f'processed_{filename1}.h5ad')
    # gm.tl.traj_inference(adata) #run before plot_stability and spectrum_full
    # gm.pl.plot_stability(adata, filename1) #creates regression file (pos and max eig values)
    # gm.pl.spectrum_full(adata, filename1) #creates spectrum file (eigen spectrum)
    # gm.tl.create_weights_geneexpress(adata, filename1)
    gm.tl.par_width_var(adata)
    gm.pl.plot_traj_widthvar(adata, filename1)


if __name__ == '__main__':
    main()
