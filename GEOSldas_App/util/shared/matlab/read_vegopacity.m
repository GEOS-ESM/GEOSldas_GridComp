
function [ veg_opacity, N_tile, avg_period ] = read_vegopacity( fname, N_tile );

% Reader for climatological mwRTM vegopacity file (little-endian, Fortran sequential binary).
%
% vegopacity files are created with the script
%   GEOSldas_GridComp/GEOSldas_App/util/inputs/mwRTM_params/Create_vegopacity_8day_clim.m
%
% The format of the vegopacity file is compatible with MAPL_ReadForcing().  There are 96 Fortran
%   records for 48 pairs of header/vegopacity entries, stored in a temporal wrap-around manner:
%
% record  1:   [FortranTag] header     [FortranTag]                      % last  record of year
% record  2:   [FortranTag] vegopacity [FortranTag] 
% record  3:   [FortranTag] header     [FortranTag]                      % first record of year
% record  4:   [FortranTag] vegopacity [FortranTag]
% ...
% record 93:   [FortranTag] header     [FortranTag]                      % last  record of year
% record 94:   [FortranTag] vegopacity [FortranTag] 
% record 95:   [FortranTag] header     [FortranTag]                      % first record of year
% record 96:   [FortranTag] vegopacity [FortranTag]
%
% header: [FortranTag] Year Month Day Hour Minute Second    ...     % start of 8-day avg period  (float)
%                      Year Month Day Hour Minute Second    ...     % end   of 8-day avg period  (float)
%                      N_data_vegopacity                    ...     % number of tiles            (float)
%                      1                                    ...                                  (float)
%         [FortranTag]
%
% The [FortranTag] is an integer that indicates how many bytes of data are stored in the
%   corresponding Fortran sequential binary record.  Integers and reals are stored in
%   simple (4-byte) precision. 
%   Odd  ("header")     records: FortranTag = 14     * 4
%   Even ("vegopactiy") records: FortranTag = N_tile * 4
%
% N_tile is *required* as input because some vegopacity files that were interpolated to cube-sphere
%   have the wrong N_tile info in their header lines.
%
% - reichle, 29 Sep 2025
%
% ------------------------------------------------------------------

MAPL_UNDEF      = 1.e15;

machfmt         = 'l';          % little-endian, GEOSldas

int_precision   = 'int32';      % precision of fortran tag and integer data
float_precision = 'float32';    % precision of real data 

disp(['read_vegopacity.m: reading from ', fname])

ifp = fopen( fname, 'r', machfmt );

N_head = 14; 
N_time = 48;

for ii=1:N_time

  disp(['reading avg period ', num2str(ii)])

  % odd records ("header")

  fortran_tag = fread( ifp,      1, int_precision   );  if fortran_tag~=N_head*4, error('wrong fortran_tag A'), end
  
  tmp         = fread( ifp, N_head, float_precision );
  
  fortran_tag = fread( ifp,      1, int_precision   );  if fortran_tag~=N_head*4, error('wrong fortran_tag B'), end

  % populate avg_period
  
  avg_period(ii,:) = tmp(1:12);
  
  % verify N_tile

  N_tile_tmp = tmp(13);
  
  if N_tile_tmp~=N_tile, disp(['WARNING: wrong N_tile in header lines: ', num2str(N_tile_tmp)]), end
    

  % even records ("vegopacity")
  
  fortran_tag = fread( ifp, 1,      int_precision   );  if fortran_tag~=N_tile*4, disp(['WARNING: wrong N_tile in vegopacity file ', num2str(fortran_tag/4)]), end
  
  tmp         = fread( ifp, N_tile, float_precision );

  fortran_tag = fread( ifp, 1,      int_precision   );  if fortran_tag~=N_tile*4, disp(['WARNING: wrong N_tile in vegopacity file ', num2str(fortran_tag/4)]), end

  % populate veg_opacity
  
  veg_opacity(ii,:) = tmp;

end

fclose(ifp);

% change no-data-values to NaN

veg_opacity( veg_opacity>MAPL_UNDEF*0.9999 ) = NaN;


disp(['done reading ', fname])

% =========== EOF ===========================================

