function [ u_bar ] = calc_freestream_matrix( alpha, n )
%CALC_FREESTREAM_MATRIX creates the freestream matrix 
%   creates a  n x 1 x 2 matrix that represents the induced velocity of
%   the freestream flow at each colocation point, based on the angle of
%   attack.
%
%   alpha -- angle of attack
%   n     -- the number of colocation points 


u_bar = ones( n, 1, 2 ); % preallocate

u_bar( :, 1, 1 ) = cos( alpha ) * u_bar( :, 1, 1 );
u_bar( :, 1, 2 ) = sin( alpha ) * u_bar( :, 1, 2 );

end % End of File

