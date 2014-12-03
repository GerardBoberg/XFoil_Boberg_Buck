% Gerard Boberg and Trevor Buck and Zane Patterson
% 22 Nov 2014
function [ B ] = calc_b_matrix_line_vortex( colocate_x, colocate_y,...
                                            vortex_x, vortex_y )
%CALC_BETA_MATRIX calculates the normalized induced velocity of a number of
%line vorticies on a number of colocation points.
%
%   colocate_x and colocate_y are 1xM vectors with the x and y locations of
%       each of the colocation points
%
%   vortex_x and vortex_y are 1xN vectors with the x and y locations of
%       each of the colocation points
%
%   B is a N x M x 2 matrix containing the normalized induced velocity
%       vectrors of each voretx panel on
%

m = length( colocate_x );
n = length( vortex_x   ) - 1;

if ( m ~= length( colocate_y ) )
    error( 'colocate_x and colocate_y must have same dimensions' );
end
if ( (n + 1) ~= length( vortex_y ) )
    error( 'vortex_x and vortex_x must have same dimensions' );
end

B = zeros( n, m, 2 );  % Pre-allocate for speed

for ii = 1:m           % Calcualtes the normalized velocity vector
    for jj = 1:n       % of vortex i on colocation point j
         [dx, dy]      = line_vortex_constant(   1, ...
                               vortex_x(ii:ii+1)  , vortex_y(ii:ii+1)  ,...
                                  colocate_x(jj)  , colocate_y(jj)   );
        B( ii, jj, 1 ) = dx;
        B( ii, jj, 2 ) = dy;
    end
end

end % End of File

