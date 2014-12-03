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
n = 5; % Number of points for the chord line. 
        % Number of vortex Panels is 2 * (n - 1)
drop_row_for_kutta = true;
debug_norms = false;

alpha = 4 * ( pi/180 ); % Angle of attack, radians

%% calculate vortex panels

% Airfoil generation
[ camber, x_panels, y_panels, trailing_edge ] = NACA4( 2, 2, 12, n );

% Find Colocation Points
[ x_colocate, y_colocate ] = find_colocation_points( x_panels, y_panels, 1/2 );

% Find the normal vectors
panel_normals  = calc_normal_vectors( x_panels, y_panels );

% Create the big matrix of coefficients on induced velocities.
AAA = calc_b_matrix_line_vortex( x_colocate, y_colocate, x_panels, y_panels );
A   = dot_coefficient_matrix( AAA, panel_normals );


% free stream b matrix
u_bar = calc_freestream_matrix( alpha, length( panel_normals) );  % sin and cos of angle attack
u_bar = dot_coefficient_matrix( u_bar, panel_normals );  
                                         % dotted with panel normals

% Handle Kutta Condition and over-determined system
if ( drop_row_for_kutta )
    % drop the row for the colocation point for the middle of the bottom
    
    new_row = zeros( 1, 2 * n - 2 );  % the kutta condition states the vortex
    new_row( 1, n-1:n ) = 1;      % strength for the top and bottom trailing
    %new_row( 1, end) = 1;                  % edge panels must be the negative of
                                  % each other. This translates to a row
                                  % of zeros with two ones at the trailing
                                  % edge panels.
    
    index = ceil(  (n) * 3 / 2 ) - 2; % length is 2 n. 1:n top, n:2n bottom
    A( index, : ) = new_row;    % middle of the bottom is n * 3/2
    u_bar( index ) = 0;           
    % row replaced, and RHS of equation zero'd out.                                                          
end



% Calculate the vortex strengths
lambdas = A \ (u_bar);

Cl = 4 * pi * sum( lambdas )

%% Calculate and plot Cl, Cm c/4, Cp distribution

%% Plot the stream lines
S = 100; % number of points per dimension of streamline plot

free_stream_x = cos( alpha );
free_stream_y = sin( alpha );

U = ones( S, S ) * free_stream_x;
V = ones( S, S ) * free_stream_y;
x_span = linspace( -1, 2, S );
y_span = linspace( -0.75, 0.75, S );

[xp, yp] = meshgrid( x_span, y_span );
%[ux, uy] = meshgrid( x_span, y_span );



for ii = 1:S     % For Each X location
    for jj = 1:S % For Each Y location
        for kk = 1:(2*n - 2) % For each vortex Panel
            
            [ dux, duy ]= line_vortex_constant(...
                                lambdas(kk),...
                                x_panels(kk:kk+1), y_panels(kk:kk+1),...
                                xp(ii,jj), yp(ii,jj) );
                            
            U( ii, jj ) = U( ii, jj ) + dux;
            V( ii, jj ) = V( ii, jj ) + duy;
        end
    end
end

 streamline_y = linspace( -0.75, 0.75, 50 );
 streamline_x = linspace( -1, -1, 50 );

%[ streamline_x, streamline_y ] = meshgrid( -1:.2:2, -.75:.15:.75 );

figure();
%plot( camber(1,:), camber(2,:), 'g-.', xxx, yyy, 'g-'  );
hold on;
plot( x_panels, y_panels, 'r--', x_colocate, y_colocate, 'go' )
streamline( xp, yp, U, V, streamline_x, streamline_y );
axis equal


if ( debug_norms)

    xx = 1:2;
    yy = 1:2;
    for ii = 1:length( panel_normals )
        xx = x_colocate( ii );
        yy = y_colocate( ii );

        xx(2) = xx(1) + panel_normals( ii, 1 );
        yy(2) = yy(1) + panel_normals( ii, 2 );
        plot( xx, yy, 'g-' )
    end
end