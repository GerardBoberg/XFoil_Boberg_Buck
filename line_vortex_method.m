function [ Cl, Cm_le, Cm_c4, Cp_dist, lambda  ] = line_vortex_method( ...
                                x_panels, y_panels, alpha, camber,...
                                kutta_drop, flip_airfoil, co_percent  )
%LINE_VORTEX_METHOD returns the relevant parameters for a line vortex
%analysis

% version F

% pull out some basic information
n_panels     = length( x_panels ) - 1;
chord_length = x_panels( (n_panels / 2) + 1 ) - x_panels(1);

% extract panel-by-panel data
delta_x = x_panels( 2:end ) - x_panels( 1: end-1 );
delta_y = y_panels( 2:end ) - y_panels( 1: end-1 );
panel_lengths = sqrt( delta_x.^2 + delta_y.^2 );

% find the colocation points
x_colocate = x_panels(1:n_panels) + co_percent .* delta_x;
y_colocate = y_panels(1:n_panels) + co_percent .* delta_y;

% find the normal vectors
normal_x      =   -delta_y ./ panel_lengths;
normal_y      =    delta_x ./ panel_lengths;

% Create the A matrix
A = zeros( n_panels, n_panels );

for ii = 1:n_panels     % the ith colocation point
    for jj = 1:n_panels % the jth vortex's induced velocity on i
        [ du, dv ] = line_vortex_constant ( 1,...
                                x_panels( jj:jj+1 ), y_panels( jj:jj+1 ),...
                                x_colocate( ii ),    y_colocate( ii )   );
        A( ii, jj ) = du * normal_x(ii) + dv * normal_y(ii);
        %disp( [num2str(du), num2str(normal_x(ii)), num2str(dv), num2str(normal_y(ii)) ] );
    end
end


% Create the u_bar matrix
u_bar = ones( n_panels, 1 );
u_bar = u_bar .* transpose( cos(alpha) * normal_x + sin(alpha) * normal_y );

% Handle Kutta Condition
te = n_panels / 2;              % find the trailing edge panel
new_row = zeros( 1, n_panels ); % row of zeros
if ( flip_airfoil )        % except for trailing edge, which are ones
    new_row(1)   = 1;
    new_row(end) = 1;
else
    new_row( te:te+1 ) = 1;         
end

if( kutta_drop )
    % drop a row from both A and u_bar for the kutta condition
    index   = ceil( n_panels * 3/4 );  % middle of the bottom of airfoil
    A( index, : )  = new_row;          % replace the row
    u_bar( index ) = 0;                % and zero out the RHS
else
    % simply apend the row to the bottom, MATLAB will least-squares approx
    A( end+1, : )  = new_row;
    u_bar( end+1 ) = 0;
end

% Solve for the lambdas
% 0 = A * lambda + u_bar;
% A^-1 * u_bar = lambda
lambda = A \ ( -u_bar ); % inverse A * u_bar

% and derive the performance results
Cl_perpanel = 2 * transpose(lambda) / chord_length;
Cl          = sum( Cl_perpanel );

Cp_dist = Cl_perpanel ./ ( panel_lengths ); 

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

