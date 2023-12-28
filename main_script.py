import scanpy as sc
import Generalized_Method_W24 as gm
from multiprocessing import freeze_support


def main():
    freeze_support()
    filename1 = 'OVCA420_TGFB1'
    gm.tl.process_and_save_data(filename1)

    adata = sc.read_h5ad(f'processed_{filename1}.h5ad')
    gm.tl.traj_inference(adata)
    gm.pl.plot_stability(adata, filename1)


if __name__ == '__main__':
    main()
