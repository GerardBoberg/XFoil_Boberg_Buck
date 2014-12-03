%% Aero 306 XFOIL Project
% Gerard Boberg, Trevor Buck, Zane Patterson
%
% 2 Dec 2014
% For Dr. Marshall's Class
%
% We implemented option C where the airfoil is made of purely vortex panels
%

clc;
clear all;
close all;

% Global Variables
n = 20; % Number of points for the chord line. 
        % Number of vortex Panels is 2 * (n - 1)
drop_row_for_kutta = true;

alpha = 4 * ( pi/180 ); % Angle of attack, radians

%% calculate vortex panels

% Airfoil generation
[ camber, x_panels, y_panels, trailing_edge ] = NACA4( 2, 2, 12, n );

% Find Colocation Points
[ x_colocate, y_colocate ] = find_colocation_points( x_panels, y_panels, 1/2 );

% Find the normal vectors
panel_normals  = calc_normal_vectors( x_panels, y_panels );

% Create the big matrix of coefficients on induced velocities.
A = calc_b_matrix_line_vortex( x_colocate, y_colocate, x_panels, y_panels );
A = dot_coefficient_matrix( A, panel_normals );


% free stream b matrix
u_bar = calc_freestream_matrix( alpha, length( panel_normals) );  % sin and cos of angle attack
u_bar = dot_coefficient_matrix( u_bar, panel_normals );  
                                         % dotted with panel normals

% Handle Kutta Condition and over-determined system
if ( drop_row_for_kutta )
    % drop the row for the colocation point for the middle of the bottom
    
    new_row = zeros( 1, 2 * n - 1 );  % the kutta condition states the vortex
    new_row( 1, n:n+1 ) = 1;      % strength for the top and bottom trailing
                                  % edge panels must be the negative of
                                  % each other. This translates to a row
                                  % of zeros with two ones at the trailing
                                  % edge panels.
    
    index = ceil(  n * 3 / 2 ); % length is 2 n. 1:n top, n:2n bottom
    A( index, : ) = new_row;    % middle of the bottom is n * 3/2
    u_bar( index ) = 0;           
    % row replaced, and RHS of equation zero'd out.                                                          
end



% Calculate the vortex strengths
lambdas = A \ u_bar;

Cl = 4 * pi * sum( lambdas )

%% Calculate and plot Cl, Cm c/4, Cp distribution

%% Plot the stream lines
S = 20; % number of points per dimension of streamline plot

free_stream_x = cos( alpha );
free_stream_y = sin( alpha );

U = ones( S, S ) * free_stream_x;
V = ones( S, S ) * free_stream_y;
xp = linspace( -1, 2, S );
yp = linspace( -0.75, 0.75, S );



for ii = 1:S     % For Each X location
    for jj = 1:S % For Each Y location
        for kk = 1:(2*n - 1) % For each vortex Panel
            
            [ dux, duy ]= line_vortex_constant(...
                                lambdas(kk),...
                                x_panels(kk:kk+1), y_panels(kk:kk+1),...
                                xp(ii), yp(jj) );
                            
            U( ii, jj ) = U( ii, jj ) + dux;
            V( ii, jj ) = V( ii, jj ) + duy;
        end
    end
end

 streamline_y = linspace( -0.75, 0.75, 20 );
 streamline_x = linspace( -1, -1, 20 );

%[ streamline_x, streamline_y ] = meshgrid( -1:.2:2, -.75:.15:.75 );


figure();
plot( camber(1,:), camber(2,:), 'g-.' );
hold on;
plot( x_panels, y_panels, 'b--', x_colocate, y_colocate, 'ro' )
streamline( xp, yp, U, V, streamline_x, streamline_y );
axis equal

UT = zeros( size( U ) );
VT = zeros( size( V ) );
for ii = 1:S
    for jj = 1:S
        [ dux, duy ] = line_vortex_constant(...
                                1,...
                                [0, 1], [-0.3, 0.3],...
                                xp(ii), yp(jj) );
        UT( ii, jj ) = UT( ii, jj ) + dux;
        VT( ii, jj ) = VT( ii, jj ) + duy;
    end
end
figure();
quiver( xp, yp, UT, VT );
hold on;
plot( [0, 1], [-0.3, 0.3], 'r' );
