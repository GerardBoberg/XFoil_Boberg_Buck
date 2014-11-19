function [ x_locations, y_locations ] = find_colocation_points( ...
                                             x_panels, y_panels, location )
%FIND_COLOCATION_POINTS_3QUARTER returns the location of the colocation
%points at the three-quarters panel point. For use with Lump Vortex.
%   inputs:
%   x_panels: array of the locations of the x-components of the start and
%               end points of the panels. Length n+1.
%   y_panels: array of the locations of the y-components of the start and
%               end points of the panels. Length n+1.
%   location: A number between 0.0 and 1.0, the spot on each panel to find 
%               the colocation point at. Can also be an array of length n+1
%               to specify a different colocation point for each panel.
%
%   outputs:
%   x_locations: array of the x-components of the colocation points. 
%                   Length n-1
%   y_locations: array of the y-components of the colocation points.
%                   Length n-1
%
%  Gerard Boberg and Trevor Buck
%  18 Nov 2014

n = length( x_panels ) - 1;

if ( n ~= ( length( y_panels ) - 1 ) )
    error( 'arrays inputs must be the same length' )
end

if ( length( location ) ~= 1 ) 
    if ( length( location ) ~= (n+1) )
        error( 'location array must be of length 1, or n+1' )
    end
end

x_locations = x_panels(1:n) .* (1 - location) + x_panels(2:(n+1)) .* location;
y_locations = y_panels(1:n) .* (1 - location) + y_panels(2:(n+1)) .* location;

end % End of File

% the above is equivelent to this code block:
%{
for ii = 1:n   % for each panel
    dx = x_panels( ii+1 ) - x_panels( ii );   % find the slope
    dy = y_panels( ii+1 ) - y_panels( ii );
    
    if( different_co )         % lerp an amount equal to the location value
        dx = location(ii) * dx;
        dy = location(ii) * dy; % for that panel specifically if given
    else                        %   different values for each panel
        dx = location * dx;
        dy = location * dy;
    end
    
    x_locations( ii ) = x_panels( ii ) + dx;
    y_locations( ii ) = y_panels( ii ) + dy;
end
%}
