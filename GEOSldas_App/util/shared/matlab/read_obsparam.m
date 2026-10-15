function [ N_obs_param, obs_param ] = read_obs_param( fname )
 
% Get observation parameters
% Format as in module enkf_types, subroutine write_obs_param
%
% Gabrielle De Lannoy  - 26 Oct 2011
%
%  1 Dec 2011 - reichle: minor modifications and check-in to CVS
%
%  8 Jun 2017 - reichle: added "flistpath" and "flistname"
%
% 17 Apr 2026 - reichle: added "fcstvarname" and "fcstunits"
%
% 15 Oct 2026 - amfox: read string fields line by line (strings can contain
%                blanks, e.g. units 'm3 m-3')
%
% ------------------------------------------------------------------     

fid = fopen(fname);

disp(['Reading ',fname]);

N_obs_param = fscanf(fid, '%d ', 1);

for i=1:N_obs_param

    obs_param(i).descr           = read_str(fid);
    obs_param(i).species         = fscanf(fid, '%f ', 1);    
    obs_param(i).orbit           = fscanf(fid, '%f ', 1);  %1=A, 2=D
    obs_param(i).pol             = fscanf(fid, '%f ', 1);  %1=H, 2=V
    obs_param(i).N_ang           = fscanf(fid, '%f ', 1);
    
    obs_param(i).ang             = fscanf(fid, '%f ', obs_param(i).N_ang);

    obs_param(i).freq            = fscanf(fid, '%f ', 1);
    obs_param(i).FOV             = fscanf(fid, '%f ', 1);
    obs_param(i).FOV_units       = read_str(fid);
    obs_param(i).assim           = read_str(fid);
    obs_param(i).scale           = read_str(fid);
    obs_param(i).getinnov        = read_str(fid);
    obs_param(i).RTM_ID          = fscanf(fid, '%f ', 1);
    obs_param(i).bias_Npar       = fscanf(fid, '%f ', 1);
    obs_param(i).bias_trel       = fscanf(fid, '%f ', 1);
    obs_param(i).bias_tcut       = fscanf(fid, '%f ', 1);
    obs_param(i).nodata          = fscanf(fid, '%f ', 1);    
    obs_param(i).varname         = read_str(fid);
    obs_param(i).units           = read_str(fid);
    obs_param(i).fcstvarname     = read_str(fid);
    obs_param(i).fcstunits       = read_str(fid);
    obs_param(i).path            = read_str(fid);
    obs_param(i).name            = read_str(fid);
    obs_param(i).maskpath        = read_str(fid);
    obs_param(i).maskname        = read_str(fid);
    obs_param(i).scalepath       = read_str(fid);
    obs_param(i).scalename       = read_str(fid);                     
    obs_param(i).flistpath       = read_str(fid);
    obs_param(i).flistname       = read_str(fid);                     
    obs_param(i).errstd          = fscanf(fid, '%f ', 1);
    obs_param(i).std_normal_max  = fscanf(fid, '%f ', 1);
    obs_param(i).zeromean        = read_str(fid);
    obs_param(i).coarsen_pert    = read_str(fid);
    obs_param(i).xcorr           = fscanf(fid, '%f ', 1);
    obs_param(i).ycorr           = fscanf(fid, '%f ', 1);
    obs_param(i).adapt           = fscanf(fid, '%f ', 1);
        
    
end

fclose(fid);

disp(['Done reading obs_param for ',num2str(N_obs_param),' species']);

% ===================================================================

function [ str ] = read_str( fid )

% read the next non-empty line and remove leading/trailing blanks and
% quotes (logicals, e.g. T or F, are not quoted)

str = '';

while isempty(str)
  str = strtrim( fgetl(fid) );
end

if str(1)=='''' && str(end)==''''
  str = str(2:end-1);
end

% =========================== EOF ====================================
