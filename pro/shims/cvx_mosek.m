function shim = cvx_mosek( shim )

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    shim.name = 'Mosek';
    shim.dualize = true;
    found = exist( [ 'mosekopt.', mexext ], 'file' );
    if ~found,
        mosek_dir = getenv( 'MOSEKLM_LICENSE_FILE' );
        if ~isempty( mosek_dir ),
            if strncmp( computer, 'PC', 2 ), fs = '\'; ps = ';'; else fs = '/'; ps = ':'; end
            k = strfind( mosek_dir, fs );
            if length( k ) > 1,
                mex_dir = [ mosek_dir(1:k(end-1)), 'toolbox', fs, 'r2009b' ];
                mex_file = [ mex_dir, fs, 'mosekopt.', mexext ];
                if exist( mex_file, 'file' ),
                    shim.path = [ mex_dir, ps ];
                    found = true;
                end
            end
            if ~found,
                shim.error = 'MOSEKLM_LICENSE_FILE is set, but no MOSEK MEX file was found.';
            end
        end
    end
    if found,
        if ~isempty( shim.path ),
            old_dir = pwd;
            cd( shim.path(1:end-1) );
        end
        try
            [rr,res]=mosekopt('minimize echo(0)',struct('c',1,'a',sparse(1,1,1),'blc',0)); %#ok
            if res.rcode ~= 0,
                shim.error = sprintf( 'error %s\n%s', res.rcodestr, res.rmsg );
            end
        catch %#ok
            errmsg = lasterror; %#ok
            shim.error = sprintf( 'unexpected MEX file failure:\n%s\n', errmsg.message );
        end
        if ~isempty( shim.path ),
            cd( old_dir );
        end
    end
    if ~found && isempty( shim.error ),
        shim.error = 'Could not find a MOSEK MEX file.';
    end
end
if isempty( shim.error ) && ~usejava('jvm'),
    shim.error = 'Java support is required to use MOSEK.';
end
if isempty( shim.error ),
    error_msg = 'A CVX Professional license is required.';
    check_license;
    shim.error = error_msg;
end
if isempty( shim.error ),
    shim.check = @check;
    shim.solve = @solve;
end

function check_license
global cvx___
try
    lic = cvx___.license;
    if isempty( lic ) || isempty( lic.signature ), return; end
    expiration = [ 10000, 100, 1 ] * sscanf( lic.expiration, '%d-%d-%d' );
    today = java.util.Date();
    today = today.getDay() + 100 * ( today.getMonth() + 100 * ( today.getYear() + 1900 ) );
    if today > expiration, return; end
    if ~isempty(lic.username)
        username = char(java.lang.System.getProperty('user.name'));
        if ~isequal(username,lic.username), return; end
    end
    if ~isempty(lic.hostid),
        found_hostid = false;
        networks = java.net.NetworkInterface.getNetworkInterfaces();
        while networks.hasMoreElements(),
            ni = networks.nextElement();
            hostid = sprintf('%02x',rem(double(ni.getHardwareAddress())+256,256));
            if isequal( hostid, lic.hostid ),
                found_hostid = true;
                break
            end
        end
        if ~found_hostid,
            return
        end
    end
    message = [ lic.name, '|', lic.organization, '|', lic.email, '|', lic.username, '|', lic.hostid, '|', lic.expiration ];
    dsa = java.security.Signature.getInstance('SHA1withDSA');
    dsa.initVerify(cvx_license('*key*'));
    dsa.update(unicode2native(message));
    if ~dsa.verify(lic.signature), return; end
    assignin('caller','error_msg','');
catch %#ok
end

function found_bad = check( nonls )
found_bad = false;
for k = 1 : length( nonls ),
if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4,
    warning( 'CVX:Mosek:Semidefinite', 'This nonlinearity requires use of semidefinite cones which Mosek does not support.\n%s', ...
        'You will need to use a different solver for this model.' );
    found_bad = true;
end
end

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )

n = numel(c);
m = numel(b);
prob     = [];
prob.c   = c;
prob.a   = At';
prob.blc = b;
prob.buc = b;
prob.blx = -Inf * ones(numel(c),1);
prob.bux = -prob.blx;
prob.ints.sub = [];
if ~isempty( nonls ),
    prob.cones = {};
end
xscale = [];
for k = 1 : length( nonls ),
    nonl = nonls(k);
    tt = nonl.type;
    ti = nonl.indices;
    ni = size( ti, 1 );
    switch tt,
        case 'i_integer',
            prob.ints.sub = [ prob.ints.sub ; ti(:) ];
        case 'i_binary',
            prob.ints.sub = [ prob.ints.sub ; ti(:) ];
            prob.blx( ti ) = 0;
            prob.bux( ti ) = 1;
        case 'i_semicontinuous',
            error( 'CVX:SolverIncompatible', 'MOSEK does not support semicontinous variables.' );
        case 'i_semiinteger',
            error( 'CVX:SolverIncompatible', 'MOSEK does not support semiinteger variables.' );
        case 'nonnegative',
            prob.blx( ti ) = 0;
        case 'lorentz',
            if ni == 1,
                prob.blx( ti ) = 0;
            else
                ti = ti([end,1:end-1],:);
                for qq = 1 : size(ti,2),
                    prob.cones{end+1}.type = 'MSK_CT_QUAD';
                    prob.cones{end}.sub  = ti(:,qq)';
                end
            end
        case 'semidefinite',
            if ni == 1,
                prob.blx( ti ) = 0;
            elseif ni == 3,
                if isempty(xscale), xscale = false(n,1); end
                ti = ti([1,3,2],:);
                xscale(ti(1:2,:)) = true; %#ok
                for qq = 1 : size(ti,2),
                    prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                    prob.cones{end}.sub  = ti(:,qq)';
                end
            else
                error( 'CVX:SolverIncompatible', 'MOSEK does not support semidefinite cones larger than 2x2.' );
            end
        case 'hermitian-semidefinite',
            if ni == 1,
                prob.blx( ti ) = 0;
            elseif ni == 4,
                if isempty(xscale), xscale = false(n,1); end
                ti = ti([1,4,2,3],:);
                xscale(ti(1:2,:)) = true; %#ok
                for qq = 1 : size(ti,2),
                    prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                    prob.cones{end}.sub  = ti(:,qq)';
                end
            else
                error( 'CVX:SolverIncompatible', 'MOSEK does not support semidefinite cones larger than 2x2.' );
            end
        case 'exponential',
            error( 'CVX:SolverIncompatible', 'MOSEK does not support the exponential cone.\nYou must use another solver for this problem.' );
        otherwise,
            error( 'Invalid cone type: %s', tt );
    end
end
if ~isempty( xscale ),
    alpha = sqrt(2.0);
    prob.c(xscale) = prob.c(xscale) * alpha;
    prob.a(:,xscale) = prob.a(:,xscale) * alpha;
end
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = prec(1);
param.MSK_DPAR_INTPNT_TOL_PFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_DFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_INFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = max(1e-14,prec(1));
command = sprintf( 'minimize echo(%d)', 3*~quiet );
[ rr, res ] = cvx_run_solver( @mosekopt, command, prob, param, 'rr', 'res', settings, 3 ); %#ok
if isfield( res.sol, 'int' ),
    sol = res.sol.int;
    has_dual = false;
elseif isfield( res.sol, 'bas' ),
    sol = res.sol.bas;
    has_dual = true;
else
    sol = res.sol.itr;
    has_dual = true;
end
tol = prec(2);
x = sol.xx;
if has_dual,
    y = sol.slc - sol.suc;
    z = sol.slx - sol.sux;
    if isfield( sol, 'snx' ), 
        z = z + sol.snx; 
    end
else
    y = NaN * ones(m,1);
    z = NaN * ones(n,1);
end
status = '';
switch sol.solsta,
    case { 'NEAR_PRIMAL_INFEASIBLE_CER', 'PRIMAL_INFEASIBLE_CER' },
        status = 'Infeasible';
        x(:) = NaN; z = z / abs(b'*y); y = y / abs(b'*y);
    case { 'NEAR_DUAL_INFEASIBLE_CER', 'DUAL_INFEASIBLE_CER' },
        status = 'Unbounded';
        y(:) = NaN; z(:) = NaN; x = x / abs(prob.c'*x);
    case { 'OPTIMAL', 'NEAR_OPTIMAL' },
        status = 'Solved';
        if sol.solsta(1) == 'N', tol = prec(3); end
    case 'UNKNOWN',
        if has_dual,
            uerr = Inf; ferr = Inf;
            if prob.c' * x < 0,
                uerr = norm(prob.c)*norm(sol.xc)/max(norm(b),1)/(-prob.c'*x);
            end
            if b' * y < 0,
                ferr = norm(b)*norm(prob.a'*y+z)/max(norm(c),1);
            end
            perr = norm( sol.xc - prob.blc ) /  ( 1 + norm( prob.blc ) );
            derr = norm( prob.c - prob.a' * y - z ) / ( 1 + norm( prob.c ) );
            gerr = max( 0, prob.c' * x - b' * y ) / max(1,abs(prob.c'*x));
            tol2 = max( [ perr, derr, gerr ] );
            tol  = min( [ tol2, ferr, uerr ] );
            if tol == tol2,
                status = 'Solved';
            elseif tol == ferr,
                status = 'Infeasible';
                x(:) = NaN; z = z / abs(b'*y); y = y / abs(b'*y);
            else
                status = 'Unbounded';
                y(:) = NaN; z(:) = NaN; x = x / abs(prob.c'*x);
            end
        else
            status = 'Failed';
            x(:) = NaN;
        end
end
if isempty(status),
    switch sol.prosta,
        case 'Failed',
            status = 'Failed'; 
            x(:) = NaN; y(:) = NaN; z(:) = NaN; 
        case 'PRIMAL_INFEASIBLE',
            status = 'Infeasible';
            x(:) = NaN; z = z / abs(b'*y); y = y / abs(b'*y);
        case 'DUAL_INFEASIBLE',
            status = 'Unbounded'; 
            y(:) = NaN; z(:) = NaN; x = x / abs(prob.c'*x);
        case 'PRIMAL_FEASIBLE',
            status = 'Solved';
        case 'PRIMAL_AND_DUAL_FEASIBLE',
            status = 'Solved';

        otherwise,
            error( 'Unknown MOSEK status: %s', sol.prosta );
    end
end
if tol > prec(3),
    status = 'Failed';
elseif sol.solsta(1) == 'N',
    tol = prec(3);
    status = [ 'Inaccurate/', status ];
elseif tol > prec(2),
    status = [ 'Inaccurate/', status ];
end
if ~isempty( xscale ),
    x(xscale) = x(xscale) * alpha;
    z(xscale) = z(xscale) / alpha;
end
iters = 0;

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.