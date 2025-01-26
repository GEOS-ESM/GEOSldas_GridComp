function [ tile_coord ] = read_til_file_nc4( fname )

% read tile coordinates from nc4 "til" file (EASEv2 and cube-sphere)
%
% reichle, 24 Jan 2025
%
% -------------------------------------------------------------

disp(['reading from ', fname])

s = ncinfo( fname );

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

  tmpdata = ncreadatt( fname, '/', this_attname );

  cmd = ['tile_coord.', this_attname, ' = tmpdata;'];

  %disp(cmd)
  eval(cmd)

end

% ---------------------------------

% read tile variables

varnames = {s.Variables.Name};

for kk=1:length(varnames)

  this_varname = varnames{kk};

  tmpdata = ncread( fname, this_varname );

  % rename some fields 
  
  if strcmp(this_varname,'pfaf_index'  ), this_varname = 'pfaf';        end
  if strcmp(this_varname,'i_indg1'     ), this_varname = 'i_indg';      end
  if strcmp(this_varname,'j_indg1'     ), this_varname = 'j_indg';      end
  if strcmp(this_varname,'frac_cell1'  ), this_varname = 'frac_cell';   end
  if strcmp(this_varname,'dummy_index1'), this_varname = 'dummy_index'; end

  cmd = ['tile_coord.', this_varname, ' = tmpdata;'];

  %disp(cmd)
  eval(cmd)

end

% ========== EOF ===================================================
