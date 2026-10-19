import numpy as np
from os import path
import sys
import struct
from datetime import datetime

def read_bin_SMOS_reg(fname=None,N_out_fields=None,read_ind_latlon=None,data_product=None,*args,**kwargs):

# Read various types of gridded files, related to SMOS.
    
# SMOS_data  = data, corresponding to the output fields [array(N_out_fields,N_grid,N_angle)]
# inc_angle  = incidence angles [array(N_angle)]
# col_ind    = EASE(v2) M36 column indices, or longitudes [array(N_grid)]
# row_ind    = EASE(v2) M36 row indices, or latitudes [array(N_grid)]
# asc_flag   = ascending or descending flag: 1=ascending, 0=descending [int]
# version    = SMOS processor version number [int]
# prep_version = gmao's pre-processor version number [int]
# start_time = start data recording in this swath [structure]
# end_time   = end data recording in this swath [structure]
# N_grid     = number of grid cells [int]
# tile_id    = optional LDASsa tile_ids [array(N_grid)]
    
# ------------------------------------------------------------------
#                 # 1 - Col-index, 0-based;
                  # 2 - Row-index, 0-based; 
                  #OR (for nearest neighbour)
                  # 1 - Lon;
                  # 2 - Lat;
    
#Tb-files:
#------------------------------------------------------------------
#N_out_fields     # 1 - Tbh 
                  # 2 - Tbv
    
    		  # 3 - heterogeneity index Tbh
                  # 4 - heterogeneity index Tbv
    
    		  # 5 - # SMOS pixels in EASE grid pixel Tbh
                  # 6 - # SMOS pixels in EASE grid pixel Tbv
    
  		  # 7 - RA Tbh
                  # 8 - RA Tbv
                  #=> repeated for T3 and T4 (9-16)
    
#SM-files:
#------------------------------------------------------------------
#N_out_fields     # 1  - SM [N float] 
                  # 2  - ST [N float] 
                  # 3  - tau [N float] 
                  # 4  - Tbh 42.5^o simulated on antenna reference frame [N float]
                  # 5  - Tbv 42.5^o simulated on antenna reference frame [N float]
                  # 6  - RSTDSM [N float]
                  # 7  - RSTDST [N float]
                  # 8  - RSTDtau [N float]
                  # 9  - stdv in SM EASE pixel (heterogeneity index) [N float]
	          # 10 - number of SMOS SM pixels per EASE grid cell [N int] 
	          # 11 - science flag [N int]
   
# ------------------------------------------------------------------

    print('reading from '+fname)
    
    if path.isfile(fname):
        with open(fname,'rb') as fin:
            din=fin.read()
        fin.close()
    else:
        sys.exit('file does not exist')
    
    # fortran tag before and after each record   
    byte_beg = 0; byte_end = 4
    
    # byte size for record
    byte_beg = byte_end
    byte_end = byte_beg + 4*3
    asc_flag,version,prep_version = struct.unpack('>'+'i'*3,din[byte_beg:byte_end])
    
    # 2 fortran tags between records
    byte_beg = byte_end; byte_end = byte_beg + 4*2
    
    byte_beg = byte_end
    byte_end = byte_beg + 4*5
    
    tmp = struct.unpack('>'+'i'*5, din[byte_beg:byte_end])
    start_time = datetime(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],0)

    # 2 fortran tags
    byte_beg = byte_end; byte_end = byte_beg + 4*2
    
    byte_beg = byte_end
    byte_end = byte_beg + 4*5
    
    tmp = struct.unpack('>'+'i'*5, din[byte_beg:byte_end])
    end_time = datetime(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],0)

    # 2 fortran tags
    byte_beg = byte_end; byte_end = byte_beg + 4*2
    
    byte_beg  = byte_end
    if data_product == 'scaling':
        byte_end = byte_beg  + 4*3
        N_grid,N_angle,N_tile = struct.uppack('>'+'i'*3, din[byte_beg:byte_end])
    else:
        byte_end = byte_beg  + 4*2
        N_grid,N_angle = struct.unpack('>'+'i'*2, din[byte_beg:byte_end])
        N_tile=1
    
    print('N_grid and N_angle : '+ str(N_grid) +', '+ str(N_angle))

    # read all records
    
    # 2 fortran tags
    byte_beg = byte_end; byte_end = byte_beg + 4*2
    
    byte_beg  = byte_end   
    if (N_grid > 1):
        byte_end = byte_beg + 4*N_angle
        inc_angle  = struct.unpack('>'+'f'*N_angle,din[byte_beg:byte_end])

        # 2 fortran tags
        byte_beg = byte_end; byte_end = byte_beg + 4*2
        
        byte_beg = byte_end 
        byte_end = byte_beg + 4*N_grid
        if read_ind_latlon == 'ind':
            col_ind=struct.unpack('>'+'i'*N_grid,din[byte_beg:byte_end])

            # 2 fortran tags
            byte_beg = byte_end; byte_end = byte_beg + 4*2
            
            byte_beg = byte_end 
            byte_end = byte_beg + 4*N_grid  
           
            row_ind=struct.unpack('>'+'i'*N_grid,din[byte_beg:byte_end])
        else:
            col_ind=struct.unpack('>'+'f'*N_grid,din[byte_beg:byte_end])

            # 2 fortran tags
            byte_beg = byte_end; byte_end = byte_beg + 4*2
            
            byte_beg = byte_end 
            byte_end = byte_beg + 4*N_grid  
            
            row_ind=struct.unpack('>'+'f'*N_grid,din[byte_beg:byte_end])
            
            if read_ind_latlon == 'latlon_id':
                    byte_beg = byte_end; byte_end = byte_beg + 4*2
                    
                    byte_beg = byte_end 
                    byte_end = byte_beg + 4*N_grid  
                    tile_id = struct.unpack('>'+'f'*N_grid,din[byte_beg:byte_end])
                   

        SMOS_data=np.full([N_out_fields,N_grid,N_angle],np.nan)
        
        for i in range(N_out_fields):
            for j in range(N_angle):
                byte_beg = byte_end; byte_end = byte_beg + 4*2
                
                byte_beg = byte_end 
                byte_end = byte_beg + 4*N_grid  
                if (((i == 4 or i == 5 or i == 11 or i == 13) and \
                     (data_product =='SCLF1C' or data_product == 'BWLF1C')) or \
                     (i == 9 and data_product == 'SMUDP2') or \
                     ((i == 4 or i == 9) and data_product == 'scaling')):
                    
                    tmp_data=struct.unpack('>'+'i'*N_grid,din[byte_beg:byte_end])
                    
                else:
                    
                    tmp_data=struct.unpack('>'+'f'*N_grid,din[byte_beg:byte_end])

                SMOS_data[i,:,j]=tmp_data
                
    else:
        SMOS_data=np.full([N_out_fields],np.nan)
        col_ind=np.nan
        row_ind=np.nan
        inc_angle=np.nan
    
    return np.array(SMOS_data), np.array(inc_angle), np.array(col_ind), np.array(row_ind), asc_flag, \
        version, prep_version, start_time, end_time, N_grid
    
