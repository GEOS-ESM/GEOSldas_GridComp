function [tile_coord ] = read_catchmentdef( fname )

% read (land) tile properties from "catchment.def" file (EASEv2 and cube-sphere)
%
% reichle, 24 Jan 2025
%
% -------------------------------------------------------------

% read file

disp(['reading from ', fname])

ifp = fopen( fname, 'r' );

% read header line

tmpdata = fscanf( ifp, '%f', 1 );

tile_coord.N_tile        = tmpdata(1);

% read rest of data

tmpdata = fscanf( ifp, '%f' );

fclose(ifp);

disp('done reading "catchment.def" file')

% --------------------------------------------------

% process data

tmpdata = reshape(tmpdata, [7, tile_coord.N_tile])';

tile_coord.tile_id       = tmpdata(:,1);
tile_coord.pfaf_index    = tmpdata(:,2); 

tile_coord.min_lon       = tmpdata(:,3);   
tile_coord.max_lon       = tmpdata(:,4);   
tile_coord.min_lat       = tmpdata(:,5);   
tile_coord.max_lat       = tmpdata(:,6);

tile_coord.elev          = tmpdata(:,7);
             
% ========== EOF ===================================================

