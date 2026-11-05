 function [ asc_flag, N_sclprm, N_angle, inc_angle,                                   ...
            lon, lat,                                                                 ...
            tbh_mean_obs, tbh_std_obs, tbh_mean_mod, tbh_std_mod, tbh_N_data,         ...
            tbv_mean_obs, tbv_std_obs, tbv_mean_mod, tbv_std_mod, tbv_N_data          ...
            ]                                                                         ...
            = read_Tb_scale_param( fname );

% Matlab reader for binary files with GEOSldas Tb obs scaling params used in SMAP L4_SM algorithm
%   (see GEOSldas Fortran subroutine scale_obs_Tb_zscore()).
%
% Example file name:  L4SM_v015_L1C_zscore_stats_D_p04.bin
%
%    Stores version 015 (v015) scaling parameters (used in L4_SM Version 8) for SMAP Tb
%    observations from descending (D) half-orbits for pentad 4 (p04; 16-20 January). 
%  
% The L4_SM algorithm runs on the 9-km EASEv2_M09 grid and assimilates SMAP L1C_Tb observations,
%    which are are on the EASEv2_M36 grid.
% The Tb observations are "scaled" to the climatology of the simulated (36-km) Tbs by subtracting the
%    the climatological mean of the SMAP Tb observations and adding the climatological mean of the
%    corresponding simulated Tbs.  The climatological mean values are seasonally varying and stored
%    in pentad scaling files.
%  
% outputs:
%
%   asc_flag      : orbit direction: 1=ascending (6pm), 0=descending (6am) [matching (A/D) in file name]
%   N_sclprm      : number of elements (essentially, L1C_TB EASEv2_M36 grid cells) in scaling parameter vectors
%   N_angle       : = 1 (number of Tb incidence angles)
%   inc_angle     : =40 (Tb incidence angle)
%
%   lon           : center longitude of 9-km (L4_SM EASEv2_M09) grid cell associated with Tb scaling parameters 
%   lat           : center latitude  of 9-km (L4_SM EASEv2_M09) grid cell associated with Tb scaling parameters
%
%   tbh_mean_obs  : mean of observed  H-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell
%   tbh_std_obs   : stdv of observed  H-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell [not used in L4_SM]
%   tbh_mean_mod  : mean of simulated H-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell
%   tbh_std_mod   : stdv of simulated H-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell [not used in L4_SM]
%   tbh_N_data    : number of data used in computation of H-pol mean and stdv values
%
%   tbv_mean_obs  : mean of observed  V-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell
%   tbv_std_obs   : stdv of observed  V-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell [not used in L4_SM]
%   tbv_mean_mod  : mean of simulated V-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell
%   tbv_std_mod   : stdv of simulated V-pol Tb for pentad and L1C_TB EASEv2_M36 grid cell [not used in L4_SM]
%   tbv_N_data    : number of data used in computation of V-pol mean and stdv values
%
% - reichle, 24 Aug 2026
%
% ------------------------------------------------------------------

int_precision   = 'int32';        % precision of fortran tag and integer data
float_precision = 'float32';      % precision of data in input file

nodata_stats    = -9999;          % no-data-value for stats data

disp(' ')
disp(['reading from ', fname])

ifp = fopen( fname, 'r', 'b' );   % *big* endian

% read header lines:

% header line 1: asc_flag (orbit direction), N_data_min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 2], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

asc_flag   = tmp_data(1);

N_data_min = tmp_data(2);                              % ignore, see "tb[X]_N_data"   

disp(' ')
disp(['asc_flag                     : ', num2str(asc_flag   )]);
disp(' ')
disp(['[ignore] N_data_min          : ', num2str(N_data_min )]);

% header line 2: start year/month/day/hour/min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 5], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

start_time.year  = tmp_data(1);
start_time.month = tmp_data(2);
start_time.day   = tmp_data(3);
start_time.hour  = tmp_data(4);
start_time.min   = tmp_data(5);

disp(['[ignore] start_time.year     : ', num2str(start_time.year )])
disp(['[ignore] start_time.month    : ', num2str(start_time.month)])
disp(['[ignore] start_time.day      : ', num2str(start_time.day  )])
disp(['[ignore] start_time.hour     : ', num2str(start_time.hour )])
disp(['[ignore] start_time.min      : ', num2str(start_time.min  )])

% header line 3: *end* year/month/day/hour/min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 5], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

end_time.year  = tmp_data(1);
end_time.month = tmp_data(2);
end_time.day   = tmp_data(3);
end_time.hour  = tmp_data(4);
end_time.min   = tmp_data(5);

disp(['[ignore] end_time.year       : ', num2str(end_time.year )])
disp(['[ignore] end_time.month      : ', num2str(end_time.month)])
disp(['[ignore] end_time.day        : ', num2str(end_time.day  )])
disp(['[ignore] end_time.hour       : ', num2str(end_time.hour )])
disp(['[ignore] end_time.min        : ', num2str(end_time.min  )])
disp(' ')

% header line 4: N_sclprm, N_angle, "N_tile"? (not read in F90]

fortran_tag  = fread( ifp, 1, int_precision );
tmp_data     = fread( ifp, [1 3], int_precision );
fortran_tag  = fread( ifp, 1, int_precision );

N_sclprm  = tmp_data(1);    % = N_tile for tile-space scaling params
N_angle   = tmp_data(2);
N_tile    = tmp_data(3);

disp(['N_sclprm                     : ', num2str(N_sclprm )]);
disp(['N_angle                      : ', num2str(N_angle)]);
disp(' ')
disp(['[ignore] N_tile              : ', num2str(N_tile )]);       
disp(' ')

if (N_angle ~= 1), error('reader only works for N_angle=1'), end

% header line 5: incidence angle(s)

fortran_tag = fread( ifp, 1,           int_precision );
inc_angle   = fread( ifp, [1 N_angle], float_precision );
fortran_tag = fread( ifp, 1,           int_precision );

disp(['inc_angle(s)                 : ', num2str(inc_angle)]);

% coord info:

% records 6-7: longitude and latitude

fortran_tag = fread( ifp, 1,            int_precision   );
lon         = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag = fread( ifp, 1,            int_precision   );

fortran_tag = fread( ifp, 1,            int_precision   );
lat         = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag = fread( ifp, 1,            int_precision   );

disp(['min max lon                  : ', num2str([min(lon), max(lon)]) ])
disp(['min max lat                  : ', num2str([min(lat), max(lat)]) ])

% record 8: tile ID  (EASEv2_M09 tile ID associated with Tb scaling parameter)


fortran_tag = fread( ifp, 1,            int_precision );
tile_id     = fread( ifp, [1 N_sclprm], int_precision );
fortran_tag = fread( ifp, 1,            int_precision );       

disp(' ')
disp(['[ignore] min max tile_id     : ', num2str([min(tile_id), max(tile_id)]) ])
disp(' ')

% records  9-13: Tb*H* mean_obs, std_obs, mean_mod, std_mod, N_data

fortran_tag  = fread( ifp, 1,            int_precision   );
tbh_mean_obs = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbh_std_obs  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbh_mean_mod = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbh_std_mod  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbh_N_data   = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

% records 14-18: Tb*V* mean_obs, std_obs, mean_mod, std_mod, N_data

fortran_tag  = fread( ifp, 1,            int_precision   );
tbv_mean_obs = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbv_std_obs  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbv_mean_mod = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbv_std_mod  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tbv_N_data   = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

% records 19-22: unused data of length N_sclprm

fortran_tag  = fread( ifp, 1,            int_precision   );
tmprecord19  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tmprecord20  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tmprecord21  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       

fortran_tag  = fread( ifp, 1,            int_precision   );
tmprecord22  = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag  = fread( ifp, 1,            int_precision   );       


% replace no-data-values with NaNs

tbh_mean_obs( tbh_mean_obs==nodata_stats ) = NaN;
tbh_std_obs ( tbh_std_obs ==nodata_stats ) = NaN;
tbh_mean_mod( tbh_mean_mod==nodata_stats ) = NaN;
tbh_std_mod ( tbh_std_mod ==nodata_stats ) = NaN;
tbh_N_data  ( tbh_N_data  ==nodata_stats ) = NaN;

tbv_mean_obs( tbv_mean_obs==nodata_stats ) = NaN;
tbv_std_obs ( tbv_std_obs ==nodata_stats ) = NaN;
tbv_mean_mod( tbv_mean_mod==nodata_stats ) = NaN;
tbv_std_mod ( tbv_std_mod ==nodata_stats ) = NaN;
tbv_N_data  ( tbv_N_data  ==nodata_stats ) = NaN;

tmprecord19 ( tmprecord19 ==nodata_stats ) = NaN;
tmprecord20 ( tmprecord20 ==nodata_stats ) = NaN;
tmprecord21 ( tmprecord21 ==nodata_stats ) = NaN;
tmprecord22 ( tmprecord22 ==nodata_stats ) = NaN;

disp(['min max tbh_mean_obs         : ', num2str([min(tbh_mean_obs), max(tbh_mean_obs)]) ])
disp(['min max tbh_std_obs          : ', num2str([min(tbh_std_obs ), max(tbh_std_obs )]) ])
disp(['min max tbh_mean_mod         : ', num2str([min(tbh_mean_mod), max(tbh_mean_mod)]) ])
disp(['min max tbh_std_mod          : ', num2str([min(tbh_std_mod ), max(tbh_std_mod )]) ])
disp(['min max tbh_N_data           : ', num2str([min(tbh_N_data  ), max(tbh_N_data  )]) ])
disp(' ')

disp(['min max tbv_mean_obs         : ', num2str([min(tbv_mean_obs), max(tbv_mean_obs)]) ])
disp(['min max tbv_std_obs          : ', num2str([min(tbv_std_obs ), max(tbv_std_obs )]) ])
disp(['min max tbv_mean_mod         : ', num2str([min(tbv_mean_mod), max(tbv_mean_mod)]) ])
disp(['min max tbv_std_mod          : ', num2str([min(tbv_std_mod ), max(tbv_std_mod )]) ])
disp(['min max tbv_N_data           : ', num2str([min(tbv_N_data  ), max(tbv_N_data  )]) ])
disp(' ')

disp(['[ignore] min max tmprecord19 : ', num2str([min(tmprecord19 ), max(tmprecord19 )]) ])
disp(['[ignore] min max tmprecord20 : ', num2str([min(tmprecord20 ), max(tmprecord20 )]) ])
disp(['[ignore] min max tmprecord21 : ', num2str([min(tmprecord21 ), max(tmprecord21 )]) ])
disp(['[ignore] min max tmprecord22 : ', num2str([min(tmprecord22 ), max(tmprecord22 )]) ])
disp(' ')


fclose(ifp);
  
% ======= EOF ==========================================================
