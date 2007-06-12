function z = times( x, y, oper )
error( nargchk( 2, 3, nargin ) );
if nargin < 3, oper = '.*'; end

%
% Check sizes
%

sx = size( x );
sy = size( y );
xs = all( sx == 1 );
ys = all( sy == 1 );
if xs,
    sz = sy;
elseif ys,
    sz = sx;
elseif ~isequal( sx, sy ),
    error( 'Matrix dimensions must agree.' );
else
    sz = sx;
end
nn = prod( sz );
zs = all( sz == 1 );

%
% Determine the computation methods
%

persistent remap_m remap_l remap_r
if isempty( remap_r ),
    % Constant result
    temp_0   = cvx_remap( 'zero' );
    temp_1   = cvx_remap( 'constant' );
    temp_2   = cvx_remap( 'valid' );
    remap_1  = ( temp_0' * temp_2 ) | ( temp_1' * temp_1 );
    remap_1  = remap_1 | remap_1';
    remap_1n = ~remap_1;
    
    % Constant * affine, Real * convex/concave/log-convex, positive * log-concave
    temp_3   = cvx_remap( 'affine' );
    temp_4   = cvx_remap( 'real' );
    temp_5   = cvx_remap( 'convex', 'concave', 'log-convex' );
    temp_6   = cvx_remap( 'positive' );
    temp_7   = cvx_remap( 'log-concave' ) .* ~temp_1;
    remap_2  = ( temp_1' * temp_3 ) | ( temp_4' * temp_5 ) | ( temp_6' * temp_7 );
    remap_2  = remap_2 .* remap_1n;
    
    % Real / log-concave, positive / log-convex
    temp_8   = cvx_remap( 'log-convex' ) .* ~temp_1;
    remap_2r = ( temp_4' * temp_7 ) | ( temp_6' * temp_8 );
    remap_2r = remap_2r .* remap_1n;
    
    % Affine * constant, Convex/concave/log-convex * real, log-concave * positive
    remap_3  = remap_2';

    % Affine * affine
    temp_9  = temp_3 .* ~temp_1;
    remap_4 = temp_9' * temp_9;

    % log-concave * log-concave, log-convex * log-convex
    remap_5  = ( temp_7' * temp_7 ) | ( temp_8' * temp_8 );
    
    % log-concave / log-convex, log-convex / log-concave
    temp_9   = temp_7' * temp_8;
    remap_5r = temp_9 | temp_9';

    remap_m = remap_1 + ( 2 * remap_2  + 3 * remap_3 + 4 * remap_4 + 5 * remap_5 ) .* remap_1n;
    remap_r = remap_1 + ( 2 * remap_2r + 3 * remap_3 + 5 * remap_5r ) .* remap_1n;
    remap_r(:,2) = 0;
    remap_l = remap_r';
end
switch oper,
    case '.*',
        remap = remap_m;
        r_recip = 0;
        l_recip = 0;
    case './',
        remap = remap_r;
        r_recip = 1;
        l_recip = 0;
    case '.\',
        remap = remap_l;
        r_recip = 0;
        l_recip = 1;
end
vx = cvx_classify( x );
vy = cvx_classify( y );
vr = remap( vx + size( remap, 1 ) * ( vy - 1 ) );
vu = unique( vr );
nv = length( vu );

%
% Process each computation type separately
%

x   = cvx( x );
y   = cvx( y );
xt  = x;
yt  = y;
xts = xs;
yts = ys;
if nv ~= 1,
    z = cvx( sz, [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    if nv ~= 1,
        t = vr == vu( k );
        if ~xs,
            xt = cvx_subsref( x, t );
            sz = size( xt );
        end
        if ~ys,
            yt = cvx_subsref( y, t );
            sz = size( yt );
        end
    end

    %
    % Apply the appropriate computation
    %

    switch vu( k ),
    case 0,

        % Invalid
        error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation: {%s} %s {%s}', cvx_class( xt, true, true ), oper, cvx_class( yt, true, true ) ) );

    case 1,

        % constant .* constant
        xt = cvx_constant( xt );
        if l_recip, xt = 1.0 ./ xt; end
        yt = cvx_constant( yt );
        if r_recip, yt = 1.0 ./ yt; end
        cvx_optval = xt .* yt;

    case 2,

        % constant .* something
        xb = cvx_constant( xt );
        if l_recip, xb = 1.0 ./ xt; end
        if r_recip, yt = exp( - log( yt ) ); end
        yb = yt.basis_;
        if ~xs,
            nn = prod( size( xb ) );
            if ys,
                xb = cvx_reshape( xb, [ 1, nn ] );
                if issparse( yb ) & ~issparse( xb ), 
                    xb = sparse( xb ); 
                end
            else
                n1 = 1 : nn;
                xb = sparse( n1, n1, xb( : ), nn, nn );
            end
        end
        cvx_optval = cvx( sz, yb * xb );

    case 3,

        % something .* constant
        if l_recip, xt = exp( - log( xt ) ); end
        xb = xt.basis_;
        yb = cvx_constant( yt );
        if r_recip, yb = 1.0 ./ yb; end
        if ~ys,
            nn = prod( size( yb ) );
            if xs,
                yb = cvx_reshape( yb, [ 1, nn ] );
                if issparse( xb ) & ~issparse( yb ),
                    yb = sparse( yb );
                end
            else
                n1 = 1 : nn;
                yb = sparse( n1, n1, yb( : ), nn, nn );
            end
        end
        cvx_optval = cvx( sz, xb * yb );

    case 4,

        % affine .* affine
        nn = prod( sz );
        xA = xt.basis_; yA = yt.basis_;
        if xs & ~ys, xA = xA( :, ones( 1, nn ) ); end
        if ys & ~xs, yA = yA( :, ones( 1, nn ) ); end
        mm = max( size( xA, 1 ), size( yA, 1 ) );
        if size( xA, 1 ) < mm, xA( mm, end ) = 0; end
        if size( yA, 1 ) < mm, yA( mm, end ) = 0; end
        xB = xA( 1, : ); xA( 1, : ) = 0;
        yB = yA( 1, : ); yA( 1, : ) = 0;
        cyA   = conj( yA );
        alpha = sum( real( xA .* yA ), 1 ) ./ max( sum( cyA .* yA, 1 ), realmin );
        adiag = sparse( 1 : nn, 1 : nn, alpha, nn, nn );
        if all( sum( abs( xA - cyA * adiag ), 2 ) <= 2 * eps * sum( abs( xA ), 2 ) ),
            beta  = xB - alpha .* conj( yB );
            alpha = reshape( alpha, sz );
            if isreal( y ),
                cvx_optval = alpha .* square( y ) + reshape( beta, sz ) .* y;
            elseif all( abs( beta ) <= 2 * eps * abs( xB ) ),
                cvx_optval = alpha .* square_abs( y );
            else
                error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): product is not real.\n' ) );
            end
        else
            error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): not a square.\n' ) );
        end

    case 5,

        % posynomial .* posynomial
        xt = log( xt );
        if l_recip, xt = - xt; end
        yt = log( yt );
        if r_recip, yt = - yt; end
        cvx_optval = exp( xt + yt );

    otherwise,

        error( 'Shouldn''t be here.' );

    end

    %
    % Store the results
    %

    if nv == 1,
        z = cvx_optval;
    else
        z = cvx_subsasgn( z, t, cvx_optval );
    end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
