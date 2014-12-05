
%% xfoil project Aero 306
% Gerard Boberg, Trevor Buck, Zane Patterson
%
% 4 Dec 2014
%
% This script runs the full vortex panel analysis a single time,
%    and then outputs Coefficient of Pressure graph, as well as 
%    streamlines and quiver plots.

% if running this script by itself, make sure to:
% clc, clear all, close all

%% parameters
n_foil              = 161;   % number of points to give to NACA4
                             %    n_panels = 2 * ( n_foil - 1 )
alpha               = 4 * (pi / 180);  % rads = degrees * degrees_to_rads
coloc_percent       = 0.5;   % uses the 50% point for colocation points


M                   = 55; % rendering resolution is M x M
debug_vort_render   = false;
num_to_roll         = 13;


kutta_drop          = false; % true drops a row, false does least-squares
finite_end          = false; % attemtps to model a finite edge. less acc
Cl_offset           = 1;     % Cl, Cm are multiplied by this number
                             %      used for calibration


%% Calculate Airfoil Parameters


% get an airfoil
[ camber, panels_x, panels_y, trailing_edge ] = NACA4( 2, 2, 12, n_foil,...
                                                            finite_end );
n_panels = length( panels_x ) - 1; % n_panels = 2 * n_foil - 2; always even

% pre-allocate
    
[ lambda_t, Cl_t, Cm_le_t, Cm_c4_t, Cp_dist_t ] = vortex_panel_analysis(...
    panels_x, panels_y, alpha, coloc_percent, kutta_drop, finite_end );
            
lambda( 1:n_panels ) = lambda_t;
Cl                    = Cl_t    * Cl_offset;
Cm_le                 = Cm_le_t * Cl_offset;
Cm_c4                 = Cm_c4_t * Cl_offset;
Cp_dist( 1:n_panels ) = Cp_dist_t;
    
disp( [ 'alpha_radians = ', num2str(alpha),...
        '      alpha_degrees = ', num2str( 180/pi*alpha ) ] );
            
%% Rendering
% output basic information
disp( [ 'Coef of Lift  = ', num2str( Cl ) ] );
disp( [ 'Coef of leading edge Moment = ', num2str( Cm_le ) ] );
disp( [ 'Coef of c/4 Moment = ', num2str( Cm_c4 ) ] );
disp( ' ' );
disp( 'Hold on, rendering streamlines...' );

% Display streamline and quiver plots
render_vortex_panels( panels_x, panels_y, lambda(:), M, alpha );
if ( debug_vort_render )
    render_vortex_panels( [0,1], [-0.1, 0.2], 1, 20 );
end

% plot Cp
% Cp_dist is periodic. Run a rolling average to smooth the curve.
Cp_dist_top =  Cp_dist( 1: floor(end/2) );
Cp_dist_bot =  Cp_dist( floor(end/2)+1:end );
Cp_dist_rolling_top = find_rolling_min( Cp_dist_top, num_to_roll );
Cp_dist_rolling_bot = find_rolling_min( Cp_dist_bot, num_to_roll);
                                  
Cp_dist_rolling_top = find_rolling_mean( Cp_dist_rolling_top  , num_to_roll );
Cp_dist_rolling_bot = find_rolling_mean( Cp_dist_rolling_bot, num_to_roll );

% We want to render at the colocation points, not the x points
coloc_x(1:n_panels) = panels_x( 1:end-1) + ...
                coloc_percent * ( panels_x( 2:end ) - panels_x( 1:end-1 ) );
            
figure();
plot(   coloc_x(1:end/2)  , Cp_dist_rolling_top, 'r',...
        coloc_x(end/2+1:end), Cp_dist_rolling_bot, 'b' );
    hold on;
plot(   coloc_x(1:end/2), Cp_dist_top, 'm',...
        coloc_x(end/2+1:end), Cp_dist_bot, 'c' );
title(  'Coefficient of Pressure Distribution' );
xlabel( 'x / chord length' )
ylabel( 'Cp = 1 - ( u / u_inf )^2' )
legend( 'Upper surface Min', 'Lower surface Min',... 
        'Upper surface actual value', 'Lower surface actual value' );
% axis( [0, 1, -5.5, 1] );
axis ij;



% End of File

