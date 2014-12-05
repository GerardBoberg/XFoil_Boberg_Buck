function [ ] = render_streamlines( x_panels, y_panels, lambda, M, alpha )
%RENDER_STREAMLINES Simply renders the streamlines


% basic prameters
n_panels = length( x_panels ) - 1;
if ( n_panels <= 1 )
    chord_length = 1;
else
    chord_length = abs( x_panels( (n_panels / 2) + 1 ) - x_panels(1) );
end
    
% setup mesh grid
x_span = linspace(  -0.5, 1.5, M );
y_span = linspace( -1, 1, ceil(M) );
[ x, y ] = meshgrid( x_span, y_span );

U = ones( size(x) ); %pre-allocate
V = ones( size(y) ); %pre-allocate

for ii = 1:M
    for jj = 1:M
        for kk = 1:length( lambda )
            [U(ii,jj), V(ii,jj)] = line_vortex_constant( lambda(kk),...
                                    x_panels(kk:kk+1), y_panels(kk:kk+1),...
                                    x(ii, jj), y(ii,jj) );
                                
            if ( nargin >= 5 )
                U(ii,jj) = U(ii,jj) + cos( alpha );
                V(ii,jj) = V(ii,jj) + sin( alpha );
            end
        end
    end
end

U

stream_x = -0.5*chord_length * ones( 1, ceil(M/2) );
stream_y = linspace( -1, 1, ceil(M/2) );


% Plot
figure();   % streamline plot
streamline( x, y, U, V, stream_x, stream_y );
hold on;
plot( x_panels, y_panels, 'r' );
title( 'streamlines' )
xlabel( 'x / c' )
ylabel( 'y / c' )

figure()    % quiver plot
quiver( x, y, U, V )
hold on;
plot( x_panels, y_panels, 'r' );
title( 'quiver' )
xlabel( 'x / c' )
ylabel( 'y / c' )

end

