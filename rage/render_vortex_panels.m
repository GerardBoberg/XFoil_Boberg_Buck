function [ ] = render_vortex_panels( panels_x, panels_y, lambda, M, alpha )
%RENDER_VORTEX_PANELS

n_panels = length( panels_x ) - 1; 

% generate coordinates
x_span = linspace( -1, 2, M );
y_span = linspace( -1, 1, M );
[ xp, yp ] = meshgrid( x_span, y_span );

% setup velocity grid
if ( nargin >= 5 )
    ux = cos(alpha) * ones( size( xp ) ); % free-stream velocity is
    vy = sin(alpha) * ones( size( yp ) ); %   < cos, sin >
else
    alpha = 0;
    ux = zeros( size( xp ) );  % no alpha given ==> no free stream
    vy = zeros( size( yp ) );
end

% calculate induced velocities at each point by each panel
for ii = 1:M      % for each point in the MxM grid
    for jj = 1:M
        for kk = 1:n_panels % for each of the k panels,
            x  = panels_x( kk:kk+1 );   % find the kth vortex panel's 
            y  = panels_y( kk:kk+1 );   % induction on the ith-jth point
            [ dux, dvy ] = line_vortex_constant_2d( lambda(kk) , x, y,...
                                            xp(ii,jj), yp(ii,jj) );
        
            ux(ii, jj) = ux(ii, jj) + dux;
            vy(ii, jj) = vy(ii, jj) + dvy;
        end
    end
end

% Quiver Plot the results
figure();
quiver( xp, yp, ux, vy );
hold on;
plot( panels_x, panels_y, 'r' );
title('Quiver plot of vortex foil in free stream')
xlabel('x / chord length')
ylabel('y / chord length')
legend(['Angle of Attack = ' num2str( alpha * 180/pi ), 'degrees'] );
axis equal;

% Streamline plot the results
stream_x = -1 * ones( 1, ceil( M/2 ) );   % starting locations of the
stream_y = linspace( -1, 1, ceil( M/2 ) );% streamlines

figure();
streamline( xp, yp, ux, vy, stream_x, stream_y );
hold on;
plot( panels_x, panels_y, 'r' );
title('Streamline plot of vortex foil in free stream')
xlabel('x / chord length')
ylabel('y / chord_length')
legend(['Angle of Attack = ' num2str( alpha * 180/pi ), 'degrees'] );
axis equal;


end % End of File

