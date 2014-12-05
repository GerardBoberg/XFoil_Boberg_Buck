% Gerard Boberg and Trevor Buck and Zane Patterson
% 18 Nov 2014
function [ A ] = dot_coefficient_matrix( b , normals )
%CREATE_COEFFICIENT_MATRIX Dots each element of the beta matrix with the 
%corrosponding normal vector.
%   inputs
%       b      : an NxMx2 matrix, where beta( i, j, : ) is a vector-2
%       normals: an Nx2 matrix, where normals( i, : ) is a vector-2 that
%           represents the surface normal for panel i.
%   outputs:
%       A   : an NxM matrix that is simply every spot in beta dotted with
%           the normal for that row.

n = size( b, 1 ); % number of rows
m = size( b, 2 ); % number of columns

if ( size( b, 3 ) ~= 2 )
    error( 'beta must be an N x M x 2 matrix' )
end

if ( size( normals, 1 ) ~= n )
    error( 'input arrays must have the same number of rows' )
end

if ( size( normals, 2 ) ~= 2 )
    error( 'normals must have 2 columns' )
end

A = zeros( n, m ); % pre-allocate

for row = 1:n % for each row
    for col = 1:m % for each column
        
        vect(1, : ) = b( row, col, : );
        norm(1, : ) = normals( row, : );
        % dot the vector at beta(i,j) with the vector normal(i)
        A( row, col ) = dot( vect, norm );
    end
end

end % end of file


