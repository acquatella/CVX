function shim = cvx_gurobi( shim )

% GUROBI/CVX_SOLVER_SHIM   Initializes a connection between CVX and Gurobi.
%
% This function verifies that Gurobi is properly installed, licensed, and
% correctly configured for use this version of MATLAB and CVX. If the
% verification succeeds, the function returns a structure containin the
% functions and data necessary to use Gurobi with CVX.
%
% When called with no arguments, CVX_SOLVER_SHIM must construct the entire
% shim structure from scratch. It also performs more complete checking.
% This more complete check is always performed by CVX_SETUP.
%
% However, once the setup has been completed, we would like to save this
% configuration information to a preference file for faster recall. The
% function handles cannot be saved, however. So when CVX_SOLVER_SHIM is
% called with a saved copy of the shim, it restores those handles.

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    shim.name = 'Gurobi';
    shim.dualize = false;
    mexname = [ 'gurobi.', mexext ];
    found = exist( mexname, 'file' );
    if ~found,
        if strncmp( computer, 'PC', 2 ), fs = '\'; ps = ';'; else fs = '/'; ps = ':'; end
        gurobi_dir = getenv( 'GUROBI_HOME' );
        if ~isempty( gurobi_dir ),
            if exist( gurobi_dir, 'dir' ),
                mpath = [ gurobi_dir, fs, 'matlab' ];
                if exist( [ mpath, fs, mexname ], 'file' ),
                    shim.path = [ mpath, ps ];
                    found = true;
                end
            end
            if ~found,
                shim.error = sprintf( 'GUROBI_HOME is set (%s), but no Gurobi MEX file was found. Please correct your installation and try again.', gurobi_dir );
            end
        end
    end
    if ~found && fs(1) == '/' && exist( '/usr/local/bin/gurobi_cl', 'file' ),
        [ status, mpath ] = unix('ls -l /usr/local/bin/gurobi_cl'); %#ok
        temp = strfind( mpath, ' -> ' );
        if ~isempty( temp ),
            mpath = mpath(temp(end)+4:end-1);
            if exist( mpath, 'file' ),
                temp = strfind( mpath, fs );
                if length( temp ) > 2,
                    mpath = [ mpath(1:temp(end-1)), 'matlab' ];
                    if exist( [ mpath, fs, mexname ], 'file' ),
                        shim.path = [ mpath, ps ];
                        found = true;
                    end
                end
            end
        end
    end
    if ~found,
        shim.error = 'Could not find a Gurobi MEX file.';
    end
end
if isempty( shim.error ),
    [ shim.error, shim.warning ] = install_check( shim );
end
if isempty( shim.error ),
    shim.check = @check;
    shim.solve = @solve;
else
    shim.check = [];
    shim.solve = [];
end

% CVX_GUROBI_ERROR
%
% If the attempt to call Gurobi fails, we need to diagnose it as best as
% possible. It might be a licensing error, or it might be an installation
% error. This routine analyzes the structure of the error message and 
% provides a formatted string to display to the user.

function [ error_msg, warning_msg ] = install_check( shim )
% Test the installation with a small sample problem
prob = struct( 'Obj', 1, 'A', sparse(1,1,1), 'Sense', '>', 'RHS', 0 );
params = struct( 'OutputFlag', false );
error_msg = '';
warning_msg = '';
if ~isempty( shim.path ), 
    cur_d = pwd; 
    cd( shim.path(1:end-1) ); 
end
myext = builtin( 'mexext' );
if isempty( regexp( which( 'gurobi' ), [ myext, '$' ] ) ), %#ok
    error_msg = sprintf( 'The original MEX file gurobi.%s cannot be run. Please reinstall Gurobi and re-run cvx_setup. If the situation persists, please contact CVX Research support.', myext );
else
    try
        res = gurobi( prob, params );
    catch errmsg
        error_msg = check_gurobi_error( errmsg );
    end
end
if ~isempty( shim.path ), 
    cd( cur_d );
end
if isempty( error_msg )
    if ~isfield( res, 'versioninfo' ) || res.versioninfo.major < 5,
        error_msg = 'CVX requires Gurobi 5.0 or later.';
    else
        switch res.versioninfo.license,
            case 1, 
                warning_msg = sprintf( 'A trial Gurobi license is detected.\nFull CVX support is enabled, but a CVX Professional license will be required to use CVX with Gurobi after the trial has completed.' );
            case 5,
                warning_msg = sprintf( 'An academic Gurobi license is detected.\nFull CVX support is enabled, but please note that a CVX Professional license is required for any non-academic use.' );
            otherwise,
                if ~usejava('jvm'),
                    error_msg = 'Java support is required to use Gurobi.';
                else
                    error_msg = 'A CVX Professional license is required.';
                    check_license;
                end
        end
    end
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

function error_msg = check_gurobi_error( error_msg )
errtxt = cvx_error( error_msg, 67, false, '    '  );
switch error_msg.identifier,
case { 'MATLAB:invalidMEXFile', 'gurobi:Error' },
    error_msg = sprintf( 'An unexpected error was encountered:\n%sPlease consult the Gurobi documentation for assistance.', errtxt );
otherwise,
    error_msg = sprintf( 'An unexpected error was encountered:\n%sPlease contact CVX support for assistance.', errtxt );
end

% GUROBI_CHECK
%
% We don't actually use this yet. The intention is to provide a gentle
% warning to the user if he selects the Gurobi solver while a model is
% being built; AND if the nonlinearities already present in the model are
% incompatible. We can also use it to provide warnings as constraints are
% entered. For now, however, we simply wait until the GUROBI_SOLVE sees 
% the problem and exits with a fatal error.

function found_bad = check( nonls )
found_bad = false;
for k = 1 : length( nonls ),
    if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4,
        warning( 'CVX:SolverIncompatible', ...
            [ 'Gurobi does not support semidefinite cones larger than 2x2.\n', ...
              'You will need to use a different solver for this model.' ] );
        found_bad = true;
        break;
    end
end

% GUROBI_SOLVE
%
% This routine accepts the problem to solve in internal CVX form and
% performs the conversions necessary for Gurobi to solve it.

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )
need_duals = nargout > 4;

n = numel(c);
m = numel(b);
prob       = [];
prob.Obj   = full(c);
prob.A     = At';
prob.RHS   = full(b);
prob.Sense = '=';
prob.LB    = -Inf * ones(n,1);
prob.UB    = -prob.LB;
prob.Cones = struct( 'Index', {} );
prob.vtype = 'C';
prob.vtype = prob.vtype(1,ones(1,n));
xscale = zeros(2,0);
for k = 1 : length( nonls ),
    nonl = nonls(k);
    tt = nonl.type;
    ti = nonl.indices;
    ni = size( ti, 1 );
    switch tt,
        case 'i_integer',
            prob.vtype( ti ) = 'I';
        case 'i_binary',
            prob.vtype( ti ) = 'B';
            prob.LB( ti ) = 0;
            prob.UB( ti ) = 1;
        case 'i_semicontinuous',
            prob.vtype( ti ) = 'S';
        case 'i_semiinteger',
            prob.vtype( ti ) = 'N';
        case 'nonnegative',
            prob.LB( ti ) = 0;
        case 'lorentz',
            if ni == 1,
                prob.LB( ti ) = 0;
            else
                ti = ti([end,1:end-1],:);
                for qq = 1 : size(ti,2),
                    prob.Cones(end+1).Index = ti(:,qq)';
                end
            end
        case 'semidefinite',
            if ni == 1,
                prob.LB( ti ) = 0;
            elseif ni == 3,
                xscale = [ xscale, ti([1,3],:) ]; %#ok
                for qq = 1 : size(ti,2),
                    prob.Cones(end+1).Index = ti(:,qq)';
                end
            else
                error( 'CVX:SolverIncompatible', 'Gurobi does not support semidefinite cones larger than 2x2.\nYou must use another solver for this problem.' );
            end
        case 'hermitian-semidefinite',
            if ni == 1,
                prob.LB( ti ) = 0;
            elseif ni == 4,
                xscale = [ xscale, ti([1,4],:) ]; %#ok
                for qq = 1 : size(ti,2),
                    prob.Cones(end+1).Index = ti(:,qq)';
                end
            else
                error( 'CVX:SolverIncompatible', 'Gurobi does not support semidefinite cones larger than 2x2.\nYou must use another solver for this problem.' );
            end
        case 'exponential',
            error( 'CVX:SolverIncompatible', 'Gurobi does not support the exponential cone.\nYou must use another solver for this problem.' );
        otherwise,
            error( 'Invalid cone type: %s', tt );
    end
end
if isempty( prob.Cones ),
    prob = rmfield( prob, 'Cones' );
end
if ~isempty( xscale ),
    nlmis = size(xscale,2);
    tempc = 1:2*nlmis; tempc = reshape([tempc;tempc],4,nlmis);
    tempr = reshape(1:2*nlmis,2,nlmis); tempr = [tempr;tempr];
    tempv = [1;1;1;-1]; tempv = tempv(:,ones(1,nlmis));
    tmat = sparse(tempr,tempc,tempv);
    xscale = xscale(:);
    prob.Obj(xscale) = tmat*prob.Obj(xscale);
    prob.A(:,xscale) = prob.A(:,xscale)*tmat;
end
params.OutputFlag = double(~quiet);
params.InfUnbdInfo = 1;
params.QCPDual = +need_duals;
params.BarConvTol = max(prec(2)/10,prec(1));
params.BarQCPConvTol = max(prec(2)/10,prec(1));
params.FeasibilityTol = max([1e-9,prec(2)/10,prec(1)]);
params.OptimalityTol = max([1e-9,prec(2)/10,prec(1)]);
try
    res = cvx_run_solver( @gurobi, prob, params, 'res', settings, 2 );
catch errmsg
    error( 'CVX:SolverError', check_gurobi_error( errmsg ) );
end
tol = prec(2);
while true,
    switch res.status,
        case { 'NUMERIC', 'INF_OR_UNBD' },
            x = NaN*ones(n,1);
            if need_duals,
                y = NaN * ones(m,1);
                z = NaN * ones(n,1);
            end
            tol = Inf;
        case 'INFEASIBLE',
            status = 'Infeasible';
            x = NaN*ones(n,1); 
            if need_duals,
                y = - res.farkasdual / abs( b' * res.farkasdual );
                z = prob.A' * y;
            end
        case 'UNBOUNDED',
            status = 'Unbounded'; 
            x = res.unbdray / abs( prob.Obj' * res.unbdray );
            if need_duals,
                y = NaN * ones(m,1);
                z = NaN * ones(n,1);
            end
        case { 'OPTIMAL', 'SUBOPTIMAL' },
            status = 'Solved';
            x = res.x;
            if need_duals,
                y = res.pi;
                z = prob.Obj - prob.A' * y;
            end
            if res.status(1) == 'S',
                if prec(1) == prec(2),
                    tol = prec(3); 
                else
                end
            end
        otherwise,
            error( 'Unknown Gurobi status: %s', res.status );
    end
    break;
end
if tol > prec(3),
    status = 'Failed';
elseif tol > prec(2),
    status = [ 'Inaccurate/', status ];
end
if ~isempty( xscale ),
    x(xscale) = tmat*x(xscale);
    if need_duals,
        z(xscale) = 0.5*tmat*z(xscale);
    end
end
iters = 0;

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
