function [ tile_coord ] = read_til_file( fname_til, fname_catchmentdef )

% read tile coordinates from tile file ("til") file (nc4 or ASCII, EASEv2 or cube-sphere) and
% put data into tile_coord structure (similar to that defined in GEOSldas)
%
% inputs:
%           fname_til          : Name (with full path) of nc4 or ASCII tile file (see examples below).
%                               
%           fname_catchmentdef : Name of matching catchment.def file.
%                                OPTIONAL, for use with ASCII tile file from bcs in legacy layout 
%                                (assembled from fname_til if tile file is from bcs in new layout).
%
% fname_til examples: 
%   '/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/geometry/EASEv2_M36/EASEv2_M36_964x406.til'
%   '/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/geometry/CF0360x6C_DE1440xPE0720/CF0360x6C_DE1440xPE0720-Pfafstetter.til'
%
% IMPORTANT: 
% This reader is designed to work for bcs versions NL5, v11, v12, and newer (beginning ~2022).  
% It should be backward-compatible for most older bcs versions, but there is no 
%   expectation that it will work universally across all old and new versions. 
%
% reichle, 27 Jan 2025
%
% -------------------------------------------------------------

disp(['reading from ', fname_til])

% detect nc4 vs ASCII

is_nc4 = 1;

if ~strcmp(fname_til(end-2:end),'nc4'), is_nc4=0; end

% -------------------------------------------------------------

if is_nc4

  % read nc4 file

  s = ncinfo( fname_til );

  % dimensions

  tile_coord.N_tile = s.Dimensions.Length;
  
  % ---------------------------------
  
  % read tile space attributes
  
  attnames = {s.Attributes.Name};
  
  for kk=1:length(attnames)
  
    this_attname = attnames{kk};
  
    % skip select attributes
    
    if strcmp(this_attname,'NCO'    ), continue, end
    if strcmp(this_attname,'history'), continue, end
    
    tmpdata = ncreadatt( fname_til, '/', this_attname );
    
    cmd = ['tile_coord.', this_attname, ' = tmpdata;'];
    
    %disp(cmd)
    eval(cmd)
    
  end

  % ---------------------------------
  
  % read tile variables
  
  varnames = {s.Variables.Name};
  
  for kk=1:length(varnames)
  
    this_varname = varnames{kk};
  
    tmpdata = ncread( fname_til, this_varname );
    
    cmd = ['tile_coord.', this_varname, ' = tmpdata;'];

    %disp(cmd)
    eval(cmd)
  
  end

else

  % -------------------------------------------------------------
  %
  % read ASCII file

  % determine fname_catchmentdef
  
  if ~exist('fname_catchmentdef','var')
  
    ind=strfind(fname_til,'/');
  
    % verify that tile file is from bcs in new layout and naming convention
  
    if ~strcmp( fname_til(ind(end-2)+1:ind(end-1)-1), 'geometry' )
    
      error('cannot derive fname_catchmentdef; try using optional input argument')
  
    end
  
    resolution = fname_til(ind(end-1)+1:ind(end)-1);
    
    fname_catchmentdef = [ fname_til(1:ind(end-2)), 'land/', resolution, '/clsm/catchment.def' ]; 
  
  end
  
  % -------------------------------------------------------------
  %
  % read tile file
  
  disp(['reading from ', fname_til])
  
  ifp = fopen( fname_til, 'r' );
  
  % -------------------------------------------------------------
  %
  % EASE grid tile space? 
  
  isEASE = 0;
  
  if ~isempty(findstr('EASE',fname_til)), isEASE=1; end
  
  % --------------
  %
  % determine number of data columns
  
  if isEASE
  
    % NOTE: very old EASE versions have N_data_col=8
    %
    % between SMAP bcs v001 (through Vv3030) and SMAP bcs v003[a] (from Vv4030), 
    %   a new (9th) column was added in the *.til file
  
    N_data_col =  9;      
  
  else
  
    N_data_col = 12;
  
  end
  
  % -------------------------------------------------------------
  %
  % Number of integers in first header line depends on bcs version:
  
  N_headerline1 = 4;   % for bcs versions after ~2022
  
  old_versions = { 'Fortuna',  ...
                   'Ganymed',  ...
                   'Heracles', ...
                   'Icarus',   ...
                   'Jason',    ...
                   'DYAMOND2', ...
                   'FV3-1_0',  ...
                   'NL3',      ...
                   'NL4',      ...
                   };
  
  for kk=1:length(old_versions) 
  
    if ~isempty( strfind( fname_til, old_versions{kk} ))
  
      N_headerline1 = 3;
  
      disp([ 'detected old bcs version ', old_versions{kk} ])
      
    end
  end
  
  % make exception for select "new" versions (stored in legacy bcs directory) that would otherwise
  %   be classified as "old" per the loop above
  
  new_versions = { 'Icarus-NLv5', ...
                   'DYAMOND',     ...
                   };
  
  for kk=1:length(new_versions) 
  
    if ~isempty( strfind( fname_til, new_versions{kk} ))
  
      N_headerline1 = 4;
  
      disp([ 'detected exception for new bcs version ', new_versions{kk} ])
      
    end
  end
  
  % echo final value of N_headerline1
  
  disp(['using N_headerline1 = ', num2str(N_headerline1) ])
  
  % -------------------------------------------------------------
  %
  % read header
  
  % header line 1: 
  
  tmpdata = fscanf( ifp, '%f', N_headerline1 );
  
  N_tile_tmp  = tmpdata(1);

  if     N_headerline1==4
  
    tile_coord.N_PfafCat = tmpdata(2);
    tile_coord.raster_nx = tmpdata(3);
    tile_coord.raster_ny = tmpdata(4);
  
  elseif N_headerline1==3
  
    tile_coord.N_PfafCat = NaN;
    tile_coord.raster_nx = tmpdata(2);
    tile_coord.raster_ny = tmpdata(3);
  
  end
  
  % --------------
  %
  % header line 2: 
  
  tile_coord.N_Grids = fscanf( ifp, '%f', 1 );     
  
  % verify N_Grids (should be 1 for EASE and 2 for non-EASE)
  
  if  isEASE & tile_coord.N_Grids~=1, error(['unexpected N_Grids for EASE tile space: ',     num2str(tile_coord.N_Grids)]), end 
  
  if ~isEASE & tile_coord.N_Grids~=2, error(['unexpected N_Grids for non-EASE tile space: ', num2str(tile_coord.N_Grids)]), end 
  
  % Deal with older EASE bcs versions having (useless) header lines for grid 2, 
  %   despite having the correct value of N_Grids=1 in header line 2.
  
  if isEASE & ( ~isempty(findstr('NL3',fname_til)) | ~isempty(findstr('NL4',fname_til)) )
  
    tmp_n_grids = 2; 
  
  else 
  
    tmp_n_grids = tile_coord.N_Grids; 
  
  end
  
  % --------------
  %
  % header lines 3-5:
  
  tile_coord.Grid_Name = fscanf( ifp, '%s', 1 );
  tile_coord.IM        = fscanf( ifp, '%f', 1 );
  tile_coord.JM        = fscanf( ifp, '%f', 1 );
  
  % --------------
  %
  % header lines 6-8 (if present):
  
  if tmp_n_grids==2
  
    % NOTE: EASE NL3 and NL4 contain additional header lines for second grid, despite (correct) N_Grids=1
    %       For these old versions, read (and then ignore) additional header lines.
  
    tile_coord.Grid_ocn_Name = fscanf( ifp, '%s', 1 );
    tile_coord.IM_ocn        = fscanf( ifp, '%f', 1 );
    tile_coord.JM_ocn        = fscanf( ifp, '%f', 1 );

  end
  
  % --------------
  %
  % read data
  
  tmpdata = fscanf( ifp, '%f' );
  
  fclose(ifp);
  
  disp('done reading "til" file')
  
  % verify N_tile*N_data_col against number of data read
  
  if length(tmpdata)~=N_tile_tmp*N_data_col, error('something wrong with N_tile*N_data_col'), end
  
  % convert into 2d array
  
  tmpdata = transpose( reshape(tmpdata, [N_data_col, N_tile_tmp]) );
  
  % --------------------------------------------------
  %
  % assign tile_id
  %
  % by convention, tile_id is equal to index number of tile in *.til file
  
  tmp_tileid = transpose(1:N_tile_tmp);
  
  % --------------------------------------------------
  %
  % subset (if requested)
  
  %ind = find(tmpdata(:,1)==100);   % keep only land
  %ind = find(tmpdata(:,1)== 19);   % keep only lakes
  %ind = find(tmpdata(:,1)== 20);   % keep only ice
  %ind = find(tmpdata(:,1)~=  0);   % keep land, lake, landice
  
  %tmpdata    = tmpdata(    ind,:);
  
  %tmp_tileid = tmp_tileid( ind, 1);
  
  % --------------------------------------------------
  %
  % copy data into tile_coord structure 
  
  tile_coord.N_tile          = size(tmpdata,1);        % number of tiles (assign here in case of subsetting above)
  
  % the following are universal (for EASE and non-EASE, all tile types)
  
  tile_coord.tile_id         = tmp_tileid;             % tile ID
  
  tile_coord.typ             = tmpdata(:, 1);          % tile type
  
  tile_coord.com_lon         = tmpdata(:, 3);          % center-of-mass longitude of tile
  tile_coord.com_lat         = tmpdata(:, 4);          % center-of-mass latitude  of tile
  
  tile_coord.i_indg          = tmpdata(:, 5);          % i index of tile on global "atm" (or EASE) grid 
  tile_coord.j_indg          = tmpdata(:, 6);          % j index of tile on global "atm" (or EASE) grid 
  tile_coord.frac_cell       = tmpdata(:, 7);          % area fraction of "atm" (or EASE) grid cell 
  
  % initialize remaining fields (to be filled below)
  
  tmpNaN = NaN*ones(tile_coord.N_tile,1);
  
  tile_coord.min_lon         = tmpNaN;                 % min longitude of tile
  tile_coord.max_lon         = tmpNaN;                 % max longitude of tile
  tile_coord.min_lat         = tmpNaN;                 % min latitude  of tile
  tile_coord.max_lat         = tmpNaN;                 % max latitude  of tile
                                                     
  tile_coord.elev            = tmpNaN;                 % elevation of tile
                                                     
  tile_coord.area            = tmpNaN;                 % area of "atm" grid cell 
  tile_coord.pfaf_index      = tmpNaN;                 % index of (hydrological) Pfafstetter catchment
  tile_coord.frac_pfaf       = tmpNaN;                 % area fraction of Pfafstetter catchment
                     
  if ~isEASE
                                 
    tile_coord.i_indg_ocn    = tmpNaN;                 % i index of tile on global "ocean" grid (grid 2)
    tile_coord.j_indg_ocn    = tmpNaN;                 % j index of tile on global "ocean" grid (grid 2)
    tile_coord.frac_cell_ocn = tmpNaN;                 % area fraction of "ocean" grid cell (grid 2)
  
  end
  
  % --------------------------------------------------
  %
  % get indices for tile types
  
  ind_land     = find( tile_coord.typ == 100);   
  ind_lake     = find( tile_coord.typ ==  19);   
  ind_landice  = find( tile_coord.typ ==  20);  
  ind_ocean    = find( tile_coord.typ ==   0);  
  
  ind_NOTland  = find( tile_coord.typ ~= 100);
  ind_NOTocean = find( tile_coord.typ ~=   0);  
  
  if isEASE
  
    col_area      = NaN;
    col_pfafindex = 2;      % land tiles only
    col_fracpfaf  = NaN;    % land tiles only 
  
  else
  
    col_area      = 2;
    col_pfafindex = 9;      % land tiles only
    col_fracpfaf  = 11; 
    
  end 
  
  
  if ~isnan(col_area)
  
    tile_coord.area        = tmpdata(:, col_area);   % tile area
    
  end
  
  tile_coord.pfaf_index( ind_NOTocean) = tmpdata(ind_NOTocean, col_pfafindex);   
  
  % -------------------------------------------------
  %
  % EASE grid tile file:
  %
  %   column  8: Pfafstetter index (same as column 2)
  %   column  9: Pfafstetter ID
  %
  % non-EASE grid tile file:
  %
  %   column  9: i_indg2    [for ocean tiles] *OR* pfaf_index [for non-ocean tiles]
  %   column 10: j_indg2
  %   column 11: frac_cell2 [for ocean tiles] *OR* frac_pfaf  [for non-ocean tiles]
  %   column 12: dummy_index2
  
  if ~isEASE
  
    tile_coord.frac_pfaf( ind_NOTocean) = tmpdata(ind_NOTocean, col_fracpfaf);
  
    tile_coord.i_indg_ocn(   ind_ocean)    = tmpdata(ind_ocean,     9);
    tile_coord.j_indg_ocn(   ind_ocean)    = tmpdata(ind_ocean,    10);
    tile_coord.frac_cell_ocn(ind_ocean)    = tmpdata(ind_ocean,    11);
  
  end
  
  % -------------------------------------------------
  %
  % get additional info from (land-only) "catchment.def" file 
  
  tctmp = read_catchmentdef( fname_catchmentdef );
  
  % verify consistency of catchment.def file and tile file
  
  if tctmp.N_tile~=length(ind_land)
  
    error('mismatch between tile file and catchment.def file: N_tile')
  
  end
  
  if ( any(tile_coord.tile_id(   ind_land) - tctmp.tile_id   )  |   ...
       any(tile_coord.pfaf_index(ind_land) - tctmp.pfaf_index)      ...
       )
  
    error('mismatch between tile file and catchment.def file: tile_id or pfaf_index')
  
  end
  
  % fill tile_coord with relevant info from catchment.def
  
  tile_coord.min_lon(ind_land) = tctmp.min_lon;
  tile_coord.max_lon(ind_land) = tctmp.max_lon;
  tile_coord.min_lat(ind_land) = tctmp.min_lat;
  tile_coord.max_lat(ind_land) = tctmp.max_lat;
  
  tile_coord.elev(   ind_land) = tctmp.elev;

end    % if is_nc4

% ========== EOF ===================================================

