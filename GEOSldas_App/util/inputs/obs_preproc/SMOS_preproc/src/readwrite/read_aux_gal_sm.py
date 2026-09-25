import numpy as np
    
def read_aux_gal_sm(filename=None):
    N_row=721
    N_col=1441
    tmp = np.fromfile(filename, dtype='<f4')
    TB_Sky_H = tmp[0:int(tmp.size/2),].reshape((N_row, N_col))
    TB_Sky_V = tmp[int(tmp.size/2):, ].reshape((N_row, N_col))
    
    return TB_Sky_H, TB_Sky_V

