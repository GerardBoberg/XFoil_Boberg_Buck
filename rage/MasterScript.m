
%% Forward and notes

clc
clear all
close all

%% parameters
L                   = 11;
n_foil              = 19; 
alpha               = (pi / 180) * linspace( 0, 10, L );
coloc_percent       = 0.75;
kutta_drop          = true;
debug_vort_render   = false;

M     = 35; % points to calculate induced velocity at for rendering


%% Calculate Airfoil Parameters

% get an airfoil
[ camber, panels_x, panels_y, trailing_edge ] = NACA4( 2, 2, 12, n_foil );

n_panels = length( panels_x ) - 1; % n_panels = 2 * n_foil - 2; always even

% pre-allocate
lambda = zeros( L, n_panels) ;
Cl     = 1:L;
Cm_le  = 1:L;
Cm_c4  = 1:L;
Cp_dist= zeros( L, n_panels) ;
for ii = 1:L
    [ lambda_t, Cl_t, Cm_le_t, Cm_c4_t, Cp_dist_t ] = vortex_panel_analysis(...
                panels_x, panels_y, alpha(ii), coloc_percent, kutta_drop );
    lambda(ii, : ) = lambda_t;
    Cl(ii)         = Cl_t / 10;
    Cm_le(ii)      = Cm_le_t;
    Cm_c4(ii)      = Cm_c4_t;
    Cp_dist(ii,:)  = Cp_dist_t;
end
            
%% Rendering
% output basic information
disp( [ 'Coef of Lift  = ', num2str( Cl(1) ) ] );
disp( [ 'Coef of c/4 Moment = ', num2str( Cm_c4(1) ) ] );

% plot Cl vs Alpha
figure();
plot( 180/pi * alpha, Cl )
title( 'Angle of Attack vs Coefficient of Lift' )
xlabel( 'Angle of Attack, degrees' )
ylabel( 'Coefficient of Lift' )
axis equal;

% Render the streamlines and Quiver
render_vortex_panels( panels_x, panels_y, lambda(1,:), M, alpha(1) );
if ( debug_vort_render )
    render_vortex_panels( [0,1], [-0.1, 0.2], 1, 20 );
end

% Plot Coefficient of Pressure
figure();
coloc_x = panels_x( 1:end-1) + ...
                coloc_percent * (panels_x( 2:end ) - panels_x( 1:end-1 ));
plot( coloc_x( 1:end/2), Cp_dist(1,1:end/2), 'r',...
    coloc_x(end/2+1:end), Cp_dist(1,end/2+1:end), 'b' );
title(  'Coefficient of Pressure Distribution' );
xlabel( 'x / chord length' )
ylabel( 'Cp = 1 - ( u / u_inf )^2' )
legend( 'Upper surface', 'Lower surface' );
axis ij;



% End of File

