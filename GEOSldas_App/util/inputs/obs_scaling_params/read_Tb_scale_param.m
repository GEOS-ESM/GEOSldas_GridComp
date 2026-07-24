 function [ asc_flag, N_data_min,                                                     ...
            start_time, end_time,                                                     ...
            N_sclprm, N_angle, inc_angle,                                             ...
            lon, lat, tile_id,                                                        ...
            tbh_mean_obs, tbh_std_obs, tbh_mean_mod, tbh_std_mod, tbh_N_data,         ...
            tbv_mean_obs, tbv_std_obs, tbv_mean_mod, tbv_std_mod, tbv_N_data          ...
            ]                                                                         ...
            = read_Tb_scale_param( fname );

% read binary files with GEOSldas Tb obs scaling params
%   (see GEOSldas Fortran subroutine scale_obs_Tb_zscore())
%  
% - reichle  22 Jun 2026
%
% ------------------------------------------------------------------

int_precision   = 'int32';        % precision of fortran tag and integer data
float_precision = 'float32';      % precision of data in input file

nodata_stats    = -9999;          % no-data-value for stats data

disp(['reading from ', fname])

ifp = fopen( fname, 'r', 'b' );   % *big* endian

% read header lines:

% header line 1: asc_flag (orbit direction), N_data_min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 2], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

asc_flag   = tmp_data(1);
N_data_min = tmp_data(2);   

disp(['asc_flag:    ', num2str(asc_flag   )]);
disp(['N_data_min:  ', num2str(N_data_min )]);

% header line 2: start year/month/day/hour/min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 5], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

start_time.year  = tmp_data(1);
start_time.month = tmp_data(2);
start_time.day   = tmp_data(3);
start_time.hour  = tmp_data(4);
start_time.min   = tmp_data(5);

disp(['start_time.year : ', num2str(start_time.year )])
disp(['start_time.month: ', num2str(start_time.month)])
disp(['start_time.day  : ', num2str(start_time.day  )])
disp(['start_time.hour : ', num2str(start_time.hour )])
disp(['start_time.min  : ', num2str(start_time.min  )])

% header line 3: *end* year/month/day/hour/min

fortran_tag = fread( ifp, 1,     int_precision );
tmp_data    = fread( ifp, [1 5], int_precision );
fortran_tag = fread( ifp, 1,     int_precision );

end_time.year  = tmp_data(1);
end_time.month = tmp_data(2);
end_time.day   = tmp_data(3);
end_time.hour  = tmp_data(4);
end_time.min   = tmp_data(5);

disp(['end_time.year : ', num2str(end_time.year )])
disp(['end_time.month: ', num2str(end_time.month)])
disp(['end_time.day  : ', num2str(end_time.day  )])
disp(['end_time.hour : ', num2str(end_time.hour )])
disp(['end_time.min  : ', num2str(end_time.min  )])

% header line 4: N_sclprm, N_angle, "N_tile"? (not read in F90]

fortran_tag  = fread( ifp, 1, int_precision );
tmp_data     = fread( ifp, [1 3], int_precision );
fortran_tag  = fread( ifp, 1, int_precision );

N_sclprm  = tmp_data(1);    % = N_tile for tile-space scaling params
N_angle   = tmp_data(2);
N_tile    = tmp_data(3);

disp(['N_sclprm:  ', num2str(N_sclprm )]);
disp(['N_angle:   ', num2str(N_angle)]);
disp(['N_tile:    ', num2str(N_tile )]);        % ignore...

if (N_angle ~= 1), error('reader only works for N_angle=1'), end

% header line 5: incidence angle(s)

fortran_tag = fread( ifp, 1,           int_precision );
inc_angle   = fread( ifp, [1 N_angle], float_precision );
fortran_tag = fread( ifp, 1,           int_precision );

disp(['inc_angle(s):  ', num2str(inc_angle)]);

% coord info:

% records 6-7: longitude and latitude

fortran_tag = fread( ifp, 1,            int_precision   );
lon         = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag = fread( ifp, 1,            int_precision   );

fortran_tag = fread( ifp, 1,            int_precision   );
lat         = fread( ifp, [1 N_sclprm], float_precision );
fortran_tag = fread( ifp, 1,            int_precision   );

disp(['min max lon: ', num2str([min(lon), max(lon)]) ])
disp(['min max lat: ', num2str([min(lat), max(lat)]) ])

% record 8: tile ID

fortran_tag = fread( ifp, 1,            int_precision );
tile_id     = fread( ifp, [1 N_sclprm], int_precision );
fortran_tag = fread( ifp, 1,            int_precision );       

disp(['min max tile_id: ', num2str([min(tile_id), max(tile_id)]) ])

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


disp(['min max tbh_mean_obs: ', num2str([min(tbh_mean_obs), max(tbh_mean_obs)]) ])
disp(['min max tbh_std_obs : ', num2str([min(tbh_std_obs ), max(tbh_std_obs )]) ])
disp(['min max tbh_mean_mod: ', num2str([min(tbh_mean_mod), max(tbh_mean_mod)]) ])
disp(['min max tbh_std_mod : ', num2str([min(tbh_std_mod ), max(tbh_std_mod )]) ])
disp(['min max tbh_N_data  : ', num2str([min(tbh_N_data  ), max(tbh_N_data  )]) ])
                                                                               
disp(['min max tbv_mean_obs: ', num2str([min(tbv_mean_obs), max(tbv_mean_obs)]) ])
disp(['min max tbv_std_obs : ', num2str([min(tbv_std_obs ), max(tbv_std_obs )]) ])
disp(['min max tbv_mean_mod: ', num2str([min(tbv_mean_mod), max(tbv_mean_mod)]) ])
disp(['min max tbv_std_mod : ', num2str([min(tbv_std_mod ), max(tbv_std_mod )]) ])
disp(['min max tbv_N_data  : ', num2str([min(tbv_N_data  ), max(tbv_N_data  )]) ])


% ------------------------------------------------
%
% records 15-22: unknown data of length N_sclprm
%
% ------------------------------------------------

fclose(ifp);
  
% ======= EOF ==========================================================
