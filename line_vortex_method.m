function [ Cl, Cm_le, Cm_c4, Cp_dist, lambda  ] = line_vortex_method( ...
                                x_panels, y_panels, alpha, camber )
%LINE_VORTEX_METHOD returns the relevant parameters for a line vortex
%analysis

% pull out some basic information
n_panels     = length( x_panels ) - 1;
chord_length = x_panels( n_panels / 2 ) - x_panels(1);

% extract panel-by-panel data
delta_x = x_panels( 2:end ) - x_panels( 1: end-1 ); 
delta_y = y_panels( 2:end ) - y_panels( 1: end-1 );
panel_lengths = sqrt( delta_x.^2 + delta_y.^2 );

% find the colocation points
x_colocate = x_panels(1:n_panels) + 0.5 .* delta_x;
y_colocate = y_panels(1:n_panels) + 0.5 .* delta_y;

% find the normal vectors
normal_x      = -delta_y ./ panel_lengths;
normal_y      =  delta_x ./ panel_lengths;

% Create the A matrix
A = zeros( n_panels, n_panels );

for ii = 1:n_panels     % the ith colocation point
    for jj = 1:n_panels % the jth vortex's induced velocity on i
        [ du, dv ] = line_vortex_constant ( 1,...
                                x_panels( jj:jj+1 ), y_panels( jj:jj+1 ),...
                                x_colocate( ii ),    y_colocate( ii )   );
        A( ii, jj ) = du * normal_x(ii) + dv * normal_y(ii);
    end
end

% Create the u_bar matrix
u_bar = ones( n_panels, 1 );
u_bar = u_bar .* transpose( cos(alpha) * normal_x + sin(alpha) * normal_y );

% Solve for the lambdas
% 0 = A * lambda + u_bar;
% A^-1 * u_bar = lambda
lambda = A \ ( -u_bar ); % inverse A * u_bar

% and derive the performance results
Cl_perpanel = 2 * transpose(lambda) / chord_length;
Cl          = sum( Cl_perpanel );

Cp_dist = Cl_perpanel ./ panel_lengths; 

Cm_le = 0;
Cm_c4 = 0;

x_le = x_panels( 1 );  % Calculating Cm_le and Cm_c4
y_le = y_panels( 1 );
x_c4 = camber( 1, ceil( end/4 ) );
y_c4 = camber( 2, ceil( end/4 ) );

for ii = 1:n_panels
    r_le  = [x_colocate(ii) - x_le, y_colocate(ii) - y_le, 0];
    r_c4  = [x_colocate(ii) - x_c4, y_colocate(ii) - y_c4, 0];
    
    Cm_le = Cm_le +  Cl_perpanel(ii) *...
                    cross( r_le, [normal_x(ii), normal_y(ii), 0] );
    
    Cm_c4 = Cm_c4 +  Cl_perpanel(ii) *...
                    cross( r_c4, [normal_x(ii), normal_y(ii), 0] );
end

Cm_le = Cm_le ./ chord_length;
Cm_c4 = Cm_c4 ./ chord_length;

end % End of File

