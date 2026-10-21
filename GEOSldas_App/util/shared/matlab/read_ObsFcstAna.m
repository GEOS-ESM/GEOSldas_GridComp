function [date_time,              ...
	  obs_assim,              ...
	  obs_species,            ...
	  obs_tilenum,            ...
	  obs_lon,                ...
	  obs_lat,                ...
	  obs_obs,                ...
 	  obs_obsvar,             ...
	  obs_fcst,               ...
 	  obs_fcstvar,            ...
	  obs_ana,                ...
 	  obs_anavar              ...
 	 ] =                      ...
    read_ObsFcstAna( fname )

%
% read_ObsFcstAna.m can be used to read "ObsFcstAna" files that 
%  contain Observations and observation-space model forecasts and 
%  analysis data
%
% data format:
%  binary: see f90 subroutine output_ObsFcstAna()    in GEOSldas
%  nc4:    see f90 subroutine write_ObsFcstAna_nc4() in GEOSldas
%
% reichle,  4 Oct 2011
%
% amfox, reichle, 21 Oct 2026 - added nc4 read functionality, removed isLDASsa flag
%
% ------------------------------------------------------------------

% initialize outputs in case file does not exist or is empty

nodata = -9999;

date_time   = struct('year',   nodata, ...
		     'month',  nodata, ...
		     'day',    nodata, ...
		     'hour',   nodata, ...
		     'min',    nodata, ...
		     'sec',    nodata, ...
		     'dofyr',  nodata, ...
		     'pentad', nodata );

obs_assim               = [];
obs_species             = [];
obs_tilenum             = [];
obs_lon                 = [];
obs_lat                 = [];
obs_obs                 = [];
obs_obsvar              = [];
obs_fcst                = [];
obs_fcstvar             = [];
obs_ana                 = [];
obs_anavar              = [];


% determine file format from file name

fname = strtrim(fname);

is_nc4 = strcmp('.nc4',fname(end-3:end));


% read file

if exist(fname)~=2

  disp(['file does not exist: ', fname])

else

  disp(['reading from ', fname  ])

  if is_nc4

    % time stamp is not stored in nc4 file, get it from file name
    
    tmp_str = regexp( fname, '(\d{8})_(\d{4})z', 'tokens', 'once' );
    
    if isempty(tmp_str)
      error(['cannot get time stamp from file name: ', fname])
    end
    
    date_time.year   = str2double( tmp_str{1}(1:4) );
    date_time.month  = str2double( tmp_str{1}(5:6) );
    date_time.day    = str2double( tmp_str{1}(7:8) );
    date_time.hour   = str2double( tmp_str{2}(1:2) );
    date_time.min    = str2double( tmp_str{2}(3:4) );
    date_time.sec    = 0;
    
    date_time        = get_dofyr_pentad( date_time );
    
    % read observation assim flag
    
    tmp_data     = double( ncread( fname, 'assim_flag' ) );
    
    obs_assim    = zeros( length(tmp_data), 1);
    obs_assim( tmp_data~= 0 ) = 1;
    
    % read obs data (convert to double)
    
    obs_species  = double( ncread( fname, 'species'    ) );
    obs_tilenum  = double( ncread( fname, 'tilenum'    ) );
    obs_lon      = double( ncread( fname, 'lon'        ) );
    obs_lat      = double( ncread( fname, 'lat'        ) );
    obs_obs      = double( ncread( fname, 'obs'        ) );
    obs_obsvar   = double( ncread( fname, 'obsvar'     ) );
    obs_fcst     = double( ncread( fname, 'fcst'       ) );
    obs_fcstvar  = double( ncread( fname, 'fcstvar'    ) );
    obs_ana      = double( ncread( fname, 'ana'        ) );
    obs_anavar   = double( ncread( fname, 'anavar'     ) );
    
  else

    % read binary file

    int_precision     = 'int32';         % precision of fortran tag
    float_precision   = 'float32';       % precision of data in input file
    logical_precision = 'int32';         % precision of data in input file

    machfmt           = 'l';             % little-endian, GEOSldas
      
    ifp = fopen( fname, 'r', machfmt );       
    
    % read N_obs and time stamp entry
    
    fortran_tag  = fread( ifp, 1, int_precision );
    N_obs        = fread( ifp, 1, int_precision );
    year         = fread( ifp, 1, int_precision );
    month        = fread( ifp, 1, int_precision );
    day          = fread( ifp, 1, int_precision );
    hour         = fread( ifp, 1, int_precision );
    minute       = fread( ifp, 1, int_precision );
    second       = fread( ifp, 1, int_precision );
    dofyr        = fread( ifp, 1, int_precision );
    pentad       = fread( ifp, 1, int_precision );
    fortran_tag  = fread( ifp, 1, int_precision );
    
    date_time.year   = year;
    date_time.month  = month;
    date_time.day    = day;
    date_time.hour   = hour;
    date_time.min    = minute;
    date_time.sec    = second;
    date_time.dofyr  = dofyr;
    date_time.pentad = pentad;
    
    % read observation assim flag
    
    fortran_tag  = fread( ifp, 1, int_precision );
    tmp_data     = fread( ifp, [N_obs 1], logical_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    obs_assim = zeros( N_obs, 1);
    obs_assim( tmp_data~= 0 ) = 1;
    
    % read species information
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_species  = fread( ifp, [N_obs 1], int_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
      
    % read tile number information
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_tilenum  = fread( ifp, [N_obs 1], int_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    % read longitude
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_lon      = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    % read latitude
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_lat      = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    
    % read observation value
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_obs      = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
      
    % read observation variance
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_obsvar   = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    
    % read observation-space model forecast value
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_fcst     = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
      
    % read observation-space model forecast variance
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_fcstvar  = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    
    % read observation-space analysis value
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_ana      = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
      
    % read observation-space analysis variance
    
    fortran_tag  = fread( ifp, 1, int_precision );
    obs_anavar   = fread( ifp, [N_obs 1], float_precision );
    fortran_tag  = fread( ifp, 1, int_precision );	
    
    % close file  
    
    fclose(ifp);

  end      % is_nc4?

  
  % no-data check:
  %  - nc4: ncread() already returns NaN for _FillValue
  %  - single ensemble member integrations yield obs_obsvar==nodata)
  %  - in some cases obs_fcst (a.k.a. Obs_pred) is no-data-value,
  %     eg. SMOS Tb when snow is present)
  
  obs_obs(     obs_obs       == nodata ) = NaN;
  obs_obsvar(  obs_obsvar    == nodata ) = NaN;
  
  obs_fcst(    obs_fcst      == nodata ) = NaN;
  obs_fcstvar( obs_fcstvar   == nodata ) = NaN;
  
  obs_ana(     obs_ana       == nodata ) = NaN;
  obs_anavar(  obs_anavar    == nodata ) = NaN;
  
end        % file exists?


% ======= EOF ==========================================================
