import numpy as np
import sys
from os import path
import struct

def write_bin_SMOS_reg(fname=None,colind=None,rowind=None,av_angle_bin=None, \
                          data=None,asc_flag=None,version=None,version_prep=None,\
                              start_time=None,end_time=None,overwrite=False,N_out_fields=None,\
                                  write_ind_latlon=None,data_product=None,tile_id=None,*args,**kwargs):

# write "fortran sequential" tile tavg files (identical to LDASsa output)
    
# optional input:
#   overwrite = 0 -- do NOT overwrite existing files, print warning
#                    message, return
#   overwrite = 1 -- overwrite existing files, print warning message
    
# ------------------------------------------------------------------
    
    #N_out_fields     # 1 - Col-index, 0-based;
                      # 2 - Row-index, 0-based; 
                      #OR (for nearest neighbout)
                      # 1 - Lon;
                      # 2 - Lat;
    
    #N_out_fields     # 1 - Tbh; 
                      # 2 - Tbv;
    
                      # 3 - heterogeneity index Tbh
                      # 4 - heterogeneity index Tbv
    
                      # 5 - # SMOS pixels in EASE grid pixel Tbh
                      # 6 - # SMOS pixels in EASE grid pixel Tbv
    
                      # 7 - RA Tbh
                      # 8 - RA Tbv
    
                      #=> repeated for T3 and T4 (9-16)
    
    #OR FOR SMUDP2:
    
                      # 1 - SM
                      # 2 - ST
                      # 3 - opacity      
                      # 4 - Tbh; 
                      # 5 - Tbv;
    
                      # 6 - SM RSTD
                      # 7 - ST RSTD
                      # 8 - opac RSTD
    
                      # 9 - stdv in SM (grid cell averaging)
                      # 10 - # small SMOS pixels inside 1 EASE grid cell;
    
                      # 11 - accumulated flag
    
                  #NOT 11 - omega; scattering albedo
                  #NOT 12 - diff_albedos (om_H-om_V) 
                  #NOT 13 - max_roughness
                  #NOT 14 - RSTD omega
                  #NOT 15 - RSTD diff_omega
                  #NOT 16 - RSTD max_roughness
    
   
    # check dimensions
    if data.shape[0] != N_out_fields:
        sys.exit('ERROR: size of data incompatible with N_out_fields')

    # check if file exists
    if path.isfile(fname):
        if not overwrite:
            sys.exit('RETURNING!!! -- NOT OVERWRITING EXISTING FILE '+fname)
            return
        else:
            print('OVERWRITING '+fname)
    else:
        print('writing '+fname)    
   
    N_grid= data.shape[1]
    N_angle=1
    if (len(data.shape) == 3):
        N_angle=data.shape[2]
        data_org=data
        if (N_angle != av_angle_bin.size):
            sys.exit('ERROR in N_angle')
                
    if ( write_ind_latlon =='latlon_id' and len(args) == 14):
        if tile_id.shape[0] != N_grid:
            sys.exit('tile_id dimensions ??')
        if tile_id.shape[1] > 1:
            print('# subgridcells per gridcell: '+ str(tile_id.shape[1]))
    
    # open file    
    fout=open(fname,'wb')
    # determine number of grid cells ; further check dimensions

    # write all records
    # length of each reoord in bytes
    fortran_tag= 3*4
    
    # write output, '>' lead the format str force big endian order
    fout.write(struct.pack('>i',fortran_tag))
    fout.write(struct.pack('>'+'i'*3,asc_flag,version,version_prep)) 
    fout.write(struct.pack('>i',fortran_tag))   
    
    fortran_tag=5*4
    fout.write(struct.pack('>i',fortran_tag))
    fout.write(struct.pack('>'+'i'*5, int(start_time.year),int(start_time.month),int(start_time.day),int(start_time.hour),int(start_time.minute)))
    fout.write(struct.pack('>i',fortran_tag))
    
    fortran_tag=5*4
    fout.write(struct.pack('>i',fortran_tag))
    fout.write(struct.pack('>'+'i'*5, int(end_time.year),int(end_time.month),int(end_time.day),int(end_time.hour),int(end_time.minute)))
    fout.write(struct.pack('>i',fortran_tag))
               
    if not (data_product == 'scaling' and write_ind_latlon =='latlon_id' and len(args) == 14):
        fortran_tag=2*4
        fout.write(struct.pack('>i',fortran_tag))
        fout.write(struct.pack('>'+'i'*2,N_grid,N_angle))
        fout.write(struct.pack('>i',fortran_tag))       
    else:
        fortran_tag=3*4
        fout.write(struct.pack('>i',fortran_tag))
        fout.write(struct.pack('>'+'i'*3,N_grid,N_angle,tile_id.shape[1]))
        fout.write(struct.pack('>i',fortran_tag))       
    
    if N_grid >= 1:
        fortran_tag=N_angle*4
        fout.write(struct.pack('>i',fortran_tag))
        fout.write(struct.pack('>'+'f'*N_angle,*av_angle_bin))
        fout.write(struct.pack('>i',fortran_tag))       
 
        fortran_tag=N_grid*4
        if write_ind_latlon == 'ind':
            fout.write(struct.pack('>i',fortran_tag))
            fout.write(struct.pack('>'+'i'*N_grid,np.round(colind)))
            fout.write(struct.pack('>i',fortran_tag)) 

            fout.write(struct.pack('>i',fortran_tag))
            fout.write(struct.pack('>'+'i'*N_grid,np.round(rowind)))
            fout.write(struct.pack('>i',fortran_tag)) 
            
        else:
            if write_ind_latlon == 'latlon':
                fout.write(struct.pack('>i',fortran_tag))
                fout.write(struct.pack('>'+'f'*N_grid, *colind))
                fout.write(struct.pack('>i',fortran_tag)) 
    
                fout.write(struct.pack('>i',fortran_tag))
                fout.write(struct.pack('>'+'f'*N_grid, *rowind))
                fout.write(struct.pack('>i',fortran_tag)) 
                
            else:
                if write_ind_latlon == 'latlon_id' and len(args) == 14:
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>'+'f'*N_grid, *colind))
                    fout.write(struct.pack('>i',fortran_tag)) 
        
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>'+'f'*N_grid, *rowind))
                    fout.write(struct.pack('>i',fortran_tag)) 
                    
                    for i in range(tile_id.shape[1]):
                        fout.write(struct.pack('>i',fortran_tag))
                        fout.write(struct.pack('>'+'i'*N_grid, *tile_id[:][i]))
                        fout.write(struct.pack('>i',fortran_tag)) 
                else:
                    sys.exit('output-arguments do not line up')
                    
        fortran_tag=N_grid*4
        for i in range(N_out_fields):
            for jj in range(N_angle):
                if N_angle > 1:
                    data=np.squeeze(data_org[:,:,jj])

                if ((i == 4 or i == 5 or i == 12 or i == 13) and (data_product == 'SCLF1C' or data_product == 'BWLF1C')) \
                    or (i == 9 and data_product == 'SMUDP2'):
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>'+'i'*N_grid, *data[i,:].round().astype('int')))
                    fout.write(struct.pack('>i',fortran_tag)) 

                else:
                    if (i == 4 or i == 9) and data_product == 'scaling':
                        fout.write(struct.pack('>i',fortran_tag))
                        fout.write(struct.pack('>'+'i'*N_grid, *data[i,:].round().astype('int')))
                        fout.write(struct.pack('>i',fortran_tag)) 
                    else:
                        fout.write(struct.pack('>i',fortran_tag))
                        fout.write(struct.pack('>'+'f'*N_grid, *data[i,:]))
                        fout.write(struct.pack('>i',fortran_tag)) 

    else:
        fortran_tag=N_angle*4
        fout.write(struct.pack('>i',fortran_tag))
        fout.write(struct.pack('>'+'f'*N_angle,*av_angle_bin))
        fout.write(struct.pack('>i',fortran_tag))       

        fortran_tag=4
        if write_ind_latlon == 'ind':
            fout.write(struct.pack('>i',fortran_tag))
            fout.write(struct.pack('>i',0))
            fout.write(struct.pack('>i',fortran_tag)) 

            fout.write(struct.pack('>i',fortran_tag))
            fout.write(struct.pack('>i',0))
            fout.write(struct.pack('>i',fortran_tag)) 
             
        else:
            if write_ind_latlon == 'latlon':
                fout.write(struct.pack('>i',fortran_tag))
                fout.write(struct.pack('>f', 0.))
                fout.write(struct.pack('>i',fortran_tag)) 
    
                fout.write(struct.pack('>i',fortran_tag))
                fout.write(struct.pack('>f', 0.))
                fout.write(struct.pack('>i',fortran_tag)) 

            else:
                if write_ind_latlon == 'latlon_id' and len(args) == 14:
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>f', 0.))
                    fout.write(struct.pack('>i',fortran_tag)) 
        
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>f', 0.))
                    fout.write(struct.pack('>i',fortran_tag)) 
                    
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>i', 0))
                    fout.write(struct.pack('>i',fortran_tag)) 
        
                else:
                        sys.exit('output-arguments do not line up')
            
        for i in range(N_out_fields):
            for jj in range(N_angle):
                if ((i == 4 or i == 5 or i == 12 or i == 13) and ( data_product == 'SCLF1C' or data_product == 'BWLF1C')) \
                      or (i == 9 and data_product == 'SMUDP2'):
                    fout.write(struct.pack('>i',fortran_tag))
                    fout.write(struct.pack('>i', 0))
                    fout.write(struct.pack('>i',fortran_tag)) 

                else:
                    if (i == 4 or i == 9) and data_product == 'scaling':
                        fout.write(struct.pack('>i',fortran_tag))
                        fout.write(struct.pack('>i', 0))
                        fout.write(struct.pack('>i',fortran_tag)) 
                    else:
                        fout.write(struct.pack('>i',fortran_tag))
                        fout.write(struct.pack('>f', -999.0))
                        fout.write(struct.pack('>i',fortran_tag))     
    
    fout.close()
