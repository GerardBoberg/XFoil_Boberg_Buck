clc;
clear all;
close all;

use_A        = false;
flip_airfoil = false;
kutta        = false;

n = 8;
n_panels = 2 * n - 2;
[ camber, x_panels, y_panels, trailing_edge ] = NACA4( 2, 2, 12, n );

if ( flip_airfoil ) % reverses le-te-te-le into te-le-le-te
    x_panels = [ x_panels( n_panels/2:end ), x_panels( 2:(n_panels/2) ) ]; 
    y_panels = [ y_panels( n_panels/2:end ), y_panels(2:(n_panels/2) ) ];
end


L = 11;
alpha = ( pi / 180 ) * linspace( 0, 10, L );

Cl        = 1:L;
Cm_le     = 1:L;
Cm_c4     = 1:L;
Cp_dist   = ones( L, n_panels );
lambda    = ones( L, n_panels );

for ii = 1:L
    if ( use_A ) 
    [ Cll, Cm_lee, Cm_c44, Cp_distt, lambdaa  ] =...
            line_vortex_method( x_panels, y_panels, alpha(ii), camber,...
            kutta, flip_airfoil );
    else
    [ Cll, Cm_lee, Cm_c44, Cp_distt, lambdaa  ] =...
            line_vortex_solver( x_panels, y_panels, alpha(ii), camber,...
            kutta, flip_airfoil );
    
    end
    Cl(ii)        = Cll;
    Cm_le(ii)     = Cm_lee(3);
    Cm_c4(ii)     = Cm_c44(3);
    Cp_dist(ii,:) = Cp_distt(:);
    lambda(ii,:)  = lambdaa(:);
end

Coef_lift_at_zero    = Cl(1)
Coef_mom_at_c4       = Cm_c4(1)
figure();
plot( x_panels(1:end-1), Cp_dist, '', 1.25*ones(1,n_panels), max( Cp_dist), 'ro' );
title( 'Cp Distribution at alpha = 4 degrees' );
xlabel( 'x-location' );
ylabel( 'Cp' );

figure();
subplot( 1, 2, 1 ), plot( alpha, Cl );
title( 'angle of attack vs. Cl' )
subplot( 1, 2, 2 ), plot( alpha, Cm_le );
title( 'angle of attack vs. Cm_c4' )

M = 35;
render_streamlines( x_panels, y_panels, lambda(1,:), M, alpha(1) );


%render_streamlines( [0,1.2], [-0.2,0.1], 1, M );


