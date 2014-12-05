
%% xfoil project Aero 306
% Gerard Boberg, Trevor Buck, Zane Patterson
%
% 4 Dec 2014

clc
clear all
close all

%% parameters
n_foil              = 9;   % number of points to give to NACA4
                             %    n_panels = 2 * ( n_foil - 1 )
L                   = 31;    % number of times to analyize during sweep
alpha               = (pi / 180) * linspace( -15, 15, L );
coloc_percent       = 0.5;   % where to place the colocation points

M     = 35; % points to calculate induced velocity at for rendering
debug_vort_render   = false;

kutta_drop          = false; % true drops a row, false does a least-squares
finite_end          = false; % creates problems, don't use
Cl_offset           = 1;     % Cl, Cm are multiplied by this. 
                             %    Used for tuning



%% Calculate Airfoil Parameters


% get an airfoil
[ camber, panels_x, panels_y, trailing_edge ] = NACA4( 2, 2, 12, n_foil,...
                                                            finite_end );
n_panels = length( panels_x ) - 1; % n_panels = 2 * n_foil - 2; always even

% pre-allocate
lambda = zeros( L, n_panels(end)) ;
Cl     = 1:L;
Cm_le  = 1:L;
Cm_c4  = 1:L;
Cp_dist= zeros( L, n_panels(end)) ;
for ii = 1:L
    
    [ lambda_t, Cl_t, Cm_le_t, Cm_c4_t, Cp_dist_t ] = vortex_panel_analysis(...
        panels_x, panels_y, alpha(ii), coloc_percent, kutta_drop, finite_end );
            
    lambda(ii, 1:n_panels ) = lambda_t;
    Cl(ii)                  = Cl_t    * Cl_offset;
    Cm_le(ii)               = Cm_le_t * Cl_offset;
    Cm_c4(ii)               = Cm_c4_t * Cl_offset;
    Cp_dist(ii, 1:n_panels) = Cp_dist_t;
    
    disp( [ 'alpha_radians = ', num2str(alpha(ii)),...
            '      alpha_degrees = ', num2str( 180/pi*alpha(ii) ) ] );
end
            
%% Rendering
% output basic information
disp( [ '--- finished sweep ---'] );
disp( [ 'Coef of Lift  = ', num2str( Cl ) ] );
disp( [ 'Coef of c/4 Moment = ', num2str( Cm_c4 ) ] );
slope_Cl = ( Cl(end) - Cl(1) ) / ( pi* ( alpha(end) - alpha(1) ) ) ;
disp( [ 'Slope of Cl = ', num2str( slope_Cl ), ' *pi' ] );

% plot Cl vs Alpha
figure();
plot( (180/pi) .* alpha, Cl )
title( 'Angle of Attack vs Coefficient of Lift' )
xlabel( 'Angle of Attack, degrees' )
ylabel( 'Coefficient of Lift' )

% plot Cm_c4 vs Alpha
figure();
plot( (180/pi) .* alpha, Cm_c4 ); %...
    %(180/pi) .* alpha, Cm_le )
title( 'Angle of Attack vs Coefficient of Moment' )
xlabel( 'Angle of Attack, degrees' )
ylabel( 'Coefficient of Moment' )
legend( 'Cm, quarter chord'); %, 'Cm, leading edge' )


% End of File

