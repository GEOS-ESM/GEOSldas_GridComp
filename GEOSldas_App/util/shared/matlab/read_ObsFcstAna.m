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
%  nc4, see f90 subroutine write_ObsFcstAna_nc4() in GEOSldas
%
% reichle,  4 Oct 2011
% amfox,   15 Oct 2026 - read nc4 instead of Fortran sequential binary
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

% read file if it exists

if exist(fname)==2

  disp(['reading from ', fname  ])

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

  % no-data check (ncread() already returns NaN for _FillValue)
  %  - single ensemble member integrations yield obs_obsvar==nodata)
  %  - in some cases obs_fcst (a.k.a. Obs_pred) is no-data-value,
  %     eg. SMOS Tb when snow is present)

  obs_obs(     obs_obs       == nodata ) = NaN;
  obs_obsvar(  obs_obsvar    == nodata ) = NaN;

  obs_fcst(    obs_fcst      == nodata ) = NaN;
  obs_fcstvar( obs_fcstvar   == nodata ) = NaN;

  obs_ana(     obs_ana       == nodata ) = NaN;
  obs_anavar(  obs_anavar    == nodata ) = NaN;

else            % if exist(fname)==2

  disp(['file does not exist: ', fname])

end


% ======= EOF ==========================================================
