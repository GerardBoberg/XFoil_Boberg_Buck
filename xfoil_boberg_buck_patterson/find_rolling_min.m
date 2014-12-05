function [ avg ] = find_rolling_min( data, num_to_average )
%FIND_ROLLING_AVERAGE Calulates a rolling minimum of points.
%
% Takes in an array of data, and returns min which is a rolling average
%   with of a given number of points to roll.
%
% num_to_average must be an odd number. If it is even, one will be added.

if( mod( num_to_average, 2 ) == 0 )
    num_to_average = num_to_average + 1;
end

r = (num_to_average - 1) / 2;

L = length( data );
avg = 1:L; % pre-allocate

for ii = 1:L
   nums_low  = max( 1, ii-r );
   nums_high = min( L, ii+r );
   
   avg( ii ) = min( data( nums_low:nums_high ) );
end


end

