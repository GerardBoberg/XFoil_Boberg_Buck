function [ Cl, Cm_le, Cm_c4, Cp_dist, lambda  ] = line_vortex_solver( ...
                                x_panels, y_panels, alpha, camber,...
                                    kutta_drop, flip_airfoil,co_percent )
%LINE_VORTEX_SOLVER returns the relevant parameters for a line vortex
%analysis

% version G

% pull out some basic information
n_panels      = length( x_panels ) - 1;
chord_length  = x_panels( 1+ n_panels / 2) - x_panels(1);

panel_lengths_x = x_panels( 2:end ) - x_panels( 1:end-1 );
panel_lengths_y = y_panels( 2:end ) - y_panels( 1:end-1 );
panel_lengths   = sqrt( panel_lengths_x .^2 + panel_lengths_y .^2 );

% Find colocation points
[ x_colocate, y_colocate ] = find_colocation_points( x_panels, y_panels, co_percent ); 

% Find the normal vectors
panel_normals = calc_normal_vectors( x_panels, y_panels );

% Create the big matrix of induced velocity coefficients
AAA = calc_b_matrix_line_vortex( x_colocate, y_colocate, x_panels, y_panels );
A   = dot_coefficient_matrix( AAA, panel_normals );

% Create free-stream b matrix
UUU   = calc_freestream_matrix( alpha, n_panels );
u_bar = dot_coefficient_matrix( UUU, panel_normals ); 

% Handle Kutta Condition
te = n_panels / 2;              % find the trailing edge panel
new_row = zeros( 1, n_panels ); % row of zeros
if ( flip_airfoil )        % except for trailing edge, which are ones
    new_row(1,1)   = 1;
    new_row(1,end) = 1;
else
    new_row( te )   = 1;
    new_row( te+1 ) = 1;
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
lambda = A\( -u_bar ); % inverse A * u_bar

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
                    cross( r_le, [ panel_normals(ii,:), 0 ] );
    
    Cm_c4 = Cm_c4 +  Cl_perpanel(ii) *...
                    cross( r_c4, [ panel_normals(ii,:), 0 ] );
end

Cm_le = Cm_le ./ chord_length;
Cm_c4 = Cm_c4 ./ chord_length;

end % End of File

