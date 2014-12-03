clc;
clear all;
close all;

n = 5;
[ camber, x_panels, y_panels, trailing_edge ] = NACA4( 2, 2, 12, n );
alpha = 0 * ( pi / 180 )

[ Cl, Cm_le, Cm_c4, Cp_dist, lambda  ] = line_vortex_method( ...
                                x_panels, y_panels, alpha, camber );
                            
Cl
Cm_le
Cm_c4


figure();
plot( x_panels(1:end-1), Cp_dist );
title( 'Cp Distribution' );
xlabel( 'x-location' );
ylabel( 'Cp' );

M = 50;
render_streamlines( x_panels, y_panels, lambda, M, alpha );


render_streamlines( [-0.5,0.5], [-0.1,0.1], 10, M );


