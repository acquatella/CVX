function z = uminus( x )

persistent remap
if isempty( remap ),
    remap = cvx_remap( 'invalid', 'log-concave' ) & ~cvx_remap( 'log-affine' );
end
tt = remap( cvx_classify( x ) );
if nnz( tt ),
    xt = cvx_subsref( x, tt );
    error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: - {%s}', cvx_class( xt ) ) );
end

z = cvx( x.size_, -x.basis_ );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
