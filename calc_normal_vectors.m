function [ panel_normals ] = calc_normal_vectors( x_points, y_points )
%CALC_NORMAL_VECTORS Takes in the airfoil panels, and calculates the
%normal of each panel.
%   inputs:
%   _____________________________________________________________
%   x_points -- an array containing the x-components of the panel start and
%                   end points. Length n+1
%   y_points -- an array containing the y-components of the panel start and
%                   end points. Length n+1
%
%   outputs:
%   ______________________________________________________________
%   panel_normals -- an nx2 array where panel_normals( i, : ) is the
%                       normal vector of that panel.
%
%
% Written by Gerard Boberg and Trevor Buck
% 18 Nov 2014

n = length( x_points ) - 1; % The number of panels on the airfoil.

if ( n ~= (length( y_points) - 1) )
    error( 'x length and y length do not match' )
end

panel_normals = ones( n, 2 ); % pre-allocate for efficieny

for ii = 1:n   % for each panel,
   delta_x = x_points( ii+1 ) - x_points( ii ); % find the current slopes 
   delta_y = y_points( ii+1 ) - y_points( ii );
   
   n_hat = [ delta_x, ( -1 / delta_y ) ]; % the normal vector is < x, -1/y >
   n_hat = n_hat / norm( n_hat );         % and is of magnitude 1
   
   panel_normals( ii, : ) = n_hat; % add to the array
end

end % End of File


