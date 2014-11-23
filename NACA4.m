% Gerard Boberg and Trevor Buck and Zane Patterson
% 18 Nov 2014
%  Adapted from Gerard's solution to hwk5, p2 for Dr. Marshall's Aero 306
function [ camber, outline_x, outline_y, trailing_edge ] = NACA4( ...
                                                            M, P, TT, n )
%NACA4 Returns the camber and outline of a 2-d NACA 4-digit airfoil
%   M   -- the first , Max Camber
%   P   -- the second digit, Location of Max Camber
%   TT  -- the third AND fourth digits, Thickness
%   n   -- the number of x-coodinates to calculate.
%
%   camber -- the points of the camber line. 
%                (1,:) are the x positions, (2,:) are the y positions
%   outline_x -- a 1xn array of x-locations of the points of the outline 
%                   of the foil.
%   outline_y -- a 1xn array of y-locations of the points of the outline 
%                   of the foil.
%                
%   trailing_edge -- the location of the trailing edge point.
%                (1) is the x-position, (2) is the y-position

m =  M / 100; % max camber
p =  P / 10;  % location of max camber
t = TT / 100; % thickness
x     = linspace( 0, 1, n );  % points to calculate

%% Calculate the camber line
yc    = 1:n;  % camber line
theta = 1:n;  % slope of the camber line
for ii = 1:n                % Camber line equation taken from Wikipedia
    if( x(ii) <= p )        % before the max camber point
        yc(ii)    =    m/p^2  * ( 2*p*x(ii) - x(ii)^2 );   
        theta(ii) = 2*(m/p^2) * ( p - x(ii) );
    else                    % after the max camber point
        yc(ii)    =   m/(1-p)^2 * ( (1-2*p) + 2*p*x(ii) - x(ii)^2 );
        theta(ii) = 2*m/(1-p)^2 * ( p - x(ii) );
    end
end

theta  = atan( theta );
camber        = [ x; yc];             % array we return for the camber line
trailing_edge = [ x(end), yc(end) ];  % the trailing edge point

%% Calculate thickness distribution
% constants a0 -> a4 from wikipedia NACA page
a0 =  0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 =  0.2843;
a4 = -0.1015;
dt = (t/0.2) * (a0 * sqrt(x) + a1 * x + a2 * x.^2 + a3 * x.^3 + a4 * x.^4 );

%% Calculate upper and lower surfaces
xu = x - dt .* sin( theta );
xl = x + dt .* sin( theta );

yu = yc + dt .* cos( theta );
yl = yc - dt .* cos( theta );

% The calculations produced and array that goes left to right on both
%   the upper and lower points.
% If layed end-to-end, it would jump to the front of the airfoil at the
%  junction between lower and upper points.
outline_x = [ xu, ( xl(end:-1:1) ) ];  % As such, reversed the lower points
outline_y = [ yu, ( yl(end:-1:1) ) ];  %  before laying arrays end-to-end


end % End of File
