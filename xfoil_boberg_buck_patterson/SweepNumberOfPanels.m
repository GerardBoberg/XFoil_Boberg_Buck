
%% xfoil project Aero 306
% Gerard Boberg, Trevor Buck, Zane Patterson
%
% 4 Dec 2014

% This file is a script that will run the full panel analysis, changing the
%    number of panels used to analyize on each iteration.
%
% warning: runs in cubic time. Can easily take sevral minutes to calculate.

clc
clear all
close all

%% parameters
L                   = 35;
n_foil              = ceil( linspace( 8, 160, L ) ); 
alpha               = 0;%(pi / 180) * linspace( 0, 10, L );
coloc_percent       = 0.5;
kutta_drop          = false;
debug_vort_render   = false;
finite_end          = true;  % toggles finite trailing edge
Cl_offset           = 1;


M     = 35; % points to calculate induced velocity at for rendering


%% Calculate Airfoil Parameters


% get an airfoil
[ camber, panels_x, panels_y, trailing_edge ] = NACA4( 2, 2, 12,...
                                                 n_foil(1), finite_end );
n_panels = length( panels_x ) - 1; % n_panels = 2 * n_foil - 2; always even

% pre-allocate
lambda = zeros( L, n_panels(end)) ;
Cl     = 1:L;
Cm_le  = 1:L;
Cm_c4  = 1:L;
Cp_dist= zeros( L, n_panels(end)) ;
for ii = 1:L
    
    % get an airfoil
    [ camber, panels_x, panels_y, trailing_edge ] = NACA4( 2, 2, 12,...
                                                 n_foil(ii), finite_end );
    n_panels = length( panels_x ) - 1;
    
    [ lambda_t, Cl_t, Cm_le_t, Cm_c4_t, Cp_dist_t ] = vortex_panel_analysis(...
                panels_x, panels_y, alpha, coloc_percent, kutta_drop, finite_end );
    lambda(ii, 1:n_panels ) = lambda_t;
    Cl(ii)         = Cl_t * Cl_offset;
    Cm_le(ii)      = Cm_le_t * Cl_offset;
    Cm_c4(ii)      = Cm_c4_t * Cl_offset;
    Cp_dist(ii, 1:n_panels)  = Cp_dist_t;
    
    disp( [ 'n_panels = ', num2str(n_panels) ] )
end
            
%% Rendering
% output basic information
disp( [ 'Coef of Lift  = ', num2str( Cl(1) ) ] );
disp( [ 'Coef of c/4 Moment = ', num2str( Cm_c4(1) ) ] );

% plot Cl vs n panels
figure();
plot( (2.*n_foil-2), Cl )
title( 'n_panels vs Coefficient of Lift' )
xlabel( 'number of panels' )
ylabel( 'Coefficient of Lift' )

% plot Cm vs n panels
figure();
plot( (2.*n_foil-2), Cm_c4 )
title( 'n_panels vs Cm c/4' )
xlabel( 'number of panels' )
ylabel( 'Coefficient of c/4 Moment' )

% Render the streamlines and Quiver
render_vortex_panels( panels_x, panels_y, lambda(end,:), M, alpha(end) );
if ( debug_vort_render )
    render_vortex_panels( [0,1], [-0.1, 0.2], 1, 20 );
end

% Plot Coefficient of Pressure



% End of File

