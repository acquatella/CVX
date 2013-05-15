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

global cvx___
if ~isempty( shim.solve ),
    return
end
shim.check = [];
shim.solve = [];
try_internal = false;
if ~usejava('jvm'),
    shim.error = 'Java support is required.';
elseif isempty( cvx___.license ),
    try_internal = true;
else
    try
        cvx___.license = full_verify( cvx___.license );
        if cvx___.license.days_left >= 0, 
            switch cvx___.license.license_type,
            case { 'academic', 'trial' },
                try_internal = true;
            otherwise,
                try_internal = any( strfind( cvx___.license.license_type, '+gurobi' ) );
            end
            if try_internal,
                hostid = cvx___.license.hostid;
                if isempty( hostid ),
                    try_internal = false;
                elseif any( strcmp( hostid, '**' ) ),
                    try_internal = any( cellfun( @(x)any(strcmp(x,hostid)), get_hostid ) );
                end
            end
        elseif cvx___.license.days_left == -Inf,
            shim.error = 'An error occurred verifying the CVX Professional license.';
        else
            shim.error = 'This CVX Professional license has expired.';
        end
    catch exc
        shim.error = my_get_report( exc );
    end
end
if ~isempty( shim.error ),
    shim.name = 'Gurobi';
    return
end
[ fs, ps, int_path, mext ] = cvx_version;
fname = [ 'gurobi.', mext ];
int_plen = length( int_path );
if isempty( shim.name ),
    shim.name = 'Gurobi';
    shim.dualize = false;
    flen = length(fname);
    fpaths = { [ int_path, fs, 'gurobi', fs, mext(4:end), fs, fname ] };
    fpaths = [ fpaths ; which( fname, '-all' ) ];
    no_native = length( fpaths ) == 1;
    switch mext,
        case 'mexmaci64', d1 = '/Library/gurobi*'; d2 = '/Library/'; d3 = 'mac64';
        case 'mexmaci',   d1 = '/Library/gurobi*'; d2 = '/Library/'; d3 = 'mac32';
        case 'mexa64',    d1 = '/opt/gurobi*'; d2 = '/opt/'; d3 = 'linx64';
        case 'mexglx',    d1 = '/opt/gurobi*'; d2 = '/opt/'; d3 = 'linux32';
        case 'mexw32',    d1 = 'C:\gurobi*'; d2 = 'C:\'; d3 = 'win32';
        case 'mexw64';    d1 = 'C:\gurobi*'; d2 = 'C:\'; d3 = 'win64';
        otherwise,        d1 = [];
    end
    if ~isempty( d1 ),
        dd = dir( d1 );
        dd = { dd([dd.isdir]).name };
        dd = strcat( strcat( d2, dd ), [ fs, d3, fs, 'matlab', fs, fname ] );
        fpaths = [ fpaths ; dd(:) ];
    end
    gurobi_dir = getenv( 'GUROBI_HOME' );
    if ~isempty( gurobi_dir ),
        fpaths{end+1} = [ gurobi_dir, fs, 'matlab', fs, fname ];
    end
    if exist( '/usr/local/bin/gurobi_cl', 'file' ),
        [ status, mpath ] = unix('ls -l /usr/local/bin/gurobi_cl'); %#ok
        mpath = regexprep( mpath, '^.*-> ', '' );
        if ~isempty( mpath ),
            mpath = regexprep( mpath, '^/+', '/' );
            mpath = regexprep( mpath, '/[^/]+/[^/]+$', '' );
            mpath = [ mpath, fs, 'matlab', fs, fname ];
            fpaths{end+1} = mpath;
        end
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    prob = struct( 'Obj', 1, 'A', sparse(1,1,1), 'Sense', '>', 'RHS', 0 );
    params = struct( 'OutputFlag', false );
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        new_dir = fpath(1:end-flen-1);
        is_internal = strncmp( new_dir, int_path, int_plen );
        if is_internal,
            tshim.location = [ '{cvx}', new_dir(int_plen+1:end-length(mext)+2) ];
        else
            tt = strfind( new_dir, fs );
            tshim.location = new_dir(1:tt(end-1)-1);
        end
        if ~exist( fpath, 'file' ),
            if exist( strrep( fpath, mext(4:end), 'a64' ), 'file' ),
                switch mext,
                    case 'mexmaci', pform = '32-bit OSX';
                    case 'glx',     pform = '32-bit Linux';
                    case 'w32',     pform = '32-bit Windows';
                    otherwise,      pform = [ 'the "', mext, '" platform' ];
                end
                tshim.error = sprintf( 'This version of Gurobi does not support %s.', pform );
                if pform(1) ~= '"',
                    tshim.error = [ tshim.error, ' If possible, please switch to the 64-bit version of MATLAB.' ];
                end
                shim = [ shim, tshim ]; %#ok
            end
            continue
        elseif any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        is_internal = strncmp( new_dir, int_path, int_plen );
        if is_internal,
            tshim.location = [ '{cvx}', new_dir(int_plen+1:end-length(mext)+2) ];
            if fs == '\', fsre = '\\'; else fsre = fs; end
            lfile = regexprep( prefdir, [ fsre, 'R\d\d\d\d\w$' ], [ fs, 'cvx_gurobi.lic' ] );
            if ~exist( lfile, 'file' ), lfile = []; end
        else
            tt = strfind( new_dir, fs );
            tshim.location = new_dir(1:tt(end-1)-1);
            lfile = [];
        end
        try
            res = gurobi5( lfile, prob, params );
        catch errmsg
            tshim.error = errmsg.message;
        end
        try
            version = res.versioninfo.major * 100 + 10 * res.versioninfo.minor + res.versioninfo.technical;
            tshim.version = sprintf( '%.2f', version / 100 );
        catch %#ok
            version = 0;
        end
        try
            ltype = res.versioninfo.license;
        catch %#ok
            ltype = 0;
        end
        tshim.params = struct( 'license', lfile, 'ltype', ltype );
            if isempty( cvx___.license ),
                switch ltype,
                case 1,
                    tshim.warning = sprintf( 'A trial Gurobi license is detected.\nFull CVX support is enabled, but a CVX Professional license will be required to use CVX with Gurobi after the trial has completed.' );
                case 5,
                    tshim.warning = sprintf( 'An academic Gurobi license is detected.\nFull CVX support is enabled, but please note that a paid CVX Professional license is required for any non-academic use.' );
                otherwise,
                    tshim.error = 'A CVX Professional license is required.';
                end
        elseif is_internal && ~try_internal,
                tshim.error = 'This CVX Professional license does not include the internal Gurobi solver.';
            end
        if isempty( tshim.error ) && version < 500,
            tshim.error = 'CVX requires Gurobi 5.0 or later.';
        end
        if isempty( tshim.error ),
            grb = @(x,y)gurobi5(lfile,x,y);
            tshim.check = @check;
            tshim.solve = @(varargin)solve(grb,varargin{:});
            if no_native || k ~= 2,
                tshim.path = [ new_dir, ps ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'Could not find a Gurobi MEX file.';
    end
else
    shim.check = [];
    shim.solve = [];
    try
        fpath = shim.fullpath;
        lfile = shim.params.license;
        ltype = shim.params.ltype;
    catch %#ok
        shim.error = 'The CVX/Gurobi interface has been updated. Please re-run CVX_SETUP.';
        return
    end
    if ~strcmp( mext, fpath(end-length(mext)+1:end) ),
        opath = regexp( shim.fullpath, 'mex\w+$', 'match' );
        if ~isempty(opath),
            opath = opath{1}(4:end);
            npath = mext(4:end);
            shim.path = strrep( shim.path, [fs,opath,ps], [fs,npath,ps] );
            shim.fullpath = [ strrep( shim.fullpath(1:end-length(opath)+1), [fs,opath,fs], [fs,npath,fs] ), npath ];
        end
    end
    if ~exist( fpath, 'file' ),
        shim.error = 'The Gurobi MEX file has been moved. Please re-run CVX_SETUP.';
    elseif isempty( shim.path ) && ~strcmp( shim.fullpath, which( fname ) ),
        shim.error = 'A new Gurobi MEX file has been installed. Please re-run CVX_SETUP.';
    elseif isempty( cvx___.license ) && ltype~= 1 && ltype ~= 5,
        shim.error = 'A CVX Professional license is required.';
    else
        grb = @(x,y)gurobi5(lfile,x,y);
        shim.check = @check;
        shim.solve = @(varargin)solve(grb,varargin{:});
    end
end

function res = gurobi5( lfile, prob, params )
if ~isempty(lfile),
    oenv = getenv('GRB_LICENSE_FILE');
    setenv('GRB_LICENSE_FILE',lfile);
    params.isvname = 'CVX';
    params.appname = 'CVX';
    params.isv_key = '';
end
try
    res = gurobi( prob, params );
catch exc
    if ~isempty( lfile ),
        setenv( 'GRB_LICENSE_FILE', oenv );
    else
        oenv = getenv( 'GRB_LICENSE_FILE' );
    end
    if strfind( exc.message, 'No Gurobi license found.' ),
        gfile = [];
        if ~isempty( lfile ),
            if exist( lfile, 'file' ),
                gfile = lfile;
                exc = sprintf([
                    'The Gurobi/CVX license found at\n',...
                    '    %s\n', ...
                    'is invalid or has expired.'], lfile );                
            else
                exc = sprintf([
                    'The Gurobi/CVX license once found at\n',...
                    '    %s\n', ...
                    'is now missing.'], lfile );
            end
        elseif isempty( oenv ),
            exc = sprintf([...
                'No valid Gurobi license was found. Possible causes include:\n',...
                '- A license has not yet been installed.\n',...
                '- The license has expired.\n',...
                '- The license was issued for a different user or host ID.\n',...
                '- It was installed in a non-standard location, and the\n',...
                '  environment variable GRB_LICENSE_FILE needs to be set.\n',...
                'Please consult the Gurobi documentation for assistance.'] );
        elseif exist( oenv, 'file'  ),
            gfile = oenv;
            exc = sprintf([...
                'The Gurobi license found at\n',...
                '    %s\n', ...
                'is invalid or has expired.'], oenv );
        else
            exc = sprintf([
                'The environment variable GRB_LICENSE_FILE is set to\n', ...
                '    %s\n', ...
                'but no Gurobi license file is found there.'], oenv );
        end
        if ~isempty( gfile ),
            fid = fopen( gfile, 'r' );
            if fid ~= 0,
                fstr = fread( fid, Inf, 'uint8=>char' )';
                fclose( fid );
                match = regexp( fstr, 'HOSTID=([0-9a-f]+)', 'once' );
                if match && ~any( cellfun( @(x)strcmp(x(5:end),fstr(match+7:match+14)), get_hostid ) ),
                    exc = sprintf( '%s\n    %s', exc, fstr(match:match+14) );
                end
                match = regexp( fstr, 'EXPIRATION=\d\d\d\d-\d\d-\d\d', 'once' );
                if match && floor(datenum(fstr(match+11:match+20),'yyyy-mm-dd'))<floor(datenum(clock))
                    exc = sprintf( '%s\n    %s', exc, fstr(match:match+20) );
                end
            end
        end
        throw( MException( 'CVX:GurobiLicense', exc ) );
    else
        throw( MException( 'CVX:GurobiError', my_get_report( exc ) ) );
    end
end
if ~isempty(lfile),
    setenv('GRB_LICENSE_FILE',oenv);
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

function [ x, status, tol, iters, y, z ] = solve( mfunc, At, b, c, nonls, quiet, prec, settings )
need_y = nargout > 4;
need_z = nargout > 5;

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
prec(1) = prec(2);
params.OutputFlag = double(~quiet);
params.InfUnbdInfo = 1;
params.QCPDual = need_y;
params.BarConvTol = prec(1);
params.BarQCPConvTol = prec(1);
params.FeasibilityTol = max([1e-9,prec(1)]);
params.OptimalityTol = max([1e-9,prec(1)]);
res = cvx_run_solver( mfunc, prob, params, 'res', settings, 2 );
tol = prec(2);
x = []; y = []; z = [];
lbound = [];
switch res.status,
    case { 'NUMERIC', 'INF_OR_UNBD' },
        tol = Inf;
    case 'INFEASIBLE',
        status = 'Infeasible';
        lbound = -Inf;
        if isfield( res, 'farkasdual' ),
            y = - res.farkasdual / abs( prob.RHS' * res.farkasdual );
            if need_z, z = prob.A' * y; end
        elseif need_y,
            tol = Inf;
        end
    case 'UNBOUNDED',
        status = 'Unbounded';
        lbound = Inf;
        if isfield( res, 'unbdray' ),
            x = res.unbdray / abs( prob.Obj' * res.unbdray );
        else
            tol = Inf;
        end
    case { 'OPTIMAL', 'SUBOPTIMAL', 'INTERRUPTED' },
        status = 'Solved';
        if isfield( res, 'x' ),
            x = res.x;
        else
            tol = Inf;
        end
        if isfield( res, 'pi' ),
            y = res.pi;
            if need_z, z = prob.Obj - prob.A' * y; end
            lbound = prob.RHS' * y;
        elseif need_y,
            tol = Inf;
        end
        if isfield( res, 'objbound' ),
            lbound = res.objbound;
        end
        if tol == Inf,
            status = 'Failed';
        elseif res.status(1) == 'S',
            status = 'Inaccurate/Solved';
        elseif res.status(1) == 'I',
            status = 'Suboptimal';
        end
    otherwise,
        tol = Inf;
        warning( 'CVX:SolverWarning', 'Gurobi returned an unknown status "%s". Please contact CVX Research Support.', res.status );
end
if isempty(x), 
    x = NaN * ones(n,1); 
end
if need_y && isempty(y),
    y = NaN * ones(m,1); 
end
if need_z && isempty(z),
    z = NaN * ones(n,1); 
end
if tol == Inf,
    status = 'Failed';
elseif tol > prec(2),
    status = [ 'Inaccurate/', status ];
end
if ~isempty( lbound ),
    tol(2) = lbound;
end
if ~isempty( xscale ),
    x(xscale) = tmat * x(xscale);
    if need_duals,
        z(xscale) = 0.5 * tmat * z(xscale);
    end
end
iters = 0;

%%%%%%%%%%%%%%%%%%%%%
% BEGIN SHIM_COMMON %
%%%%%%%%%%%%%%%%%%%%%

function public_key = get_public_key
try
    public_key = cvx_license( '*key*' );
catch %#ok
    public_key = int8(0);
end

%%%%%%%%%%%%%%%%
% BEGIN COMMON %
%%%%%%%%%%%%%%%%

function lic = full_verify( lic )
try
    signature = lic.signature;
    lic.signature = [];
    parser = java.text.SimpleDateFormat('yyyy-MM-dd');
    expire = parser.parse(lic.expiration);
    today = java.util.Date;
    days_left = ceil( ( double(expire.getTime()) - double(today.getTime) ) / 86400000 );
    if ~isfield( lic, 'prefix' ),
        lic.prefix = '';
    end
    if isempty( lic.username ),
        t_username = '';
    elseif ischar( lic.username ),
        t_username = lic.username;
    else
        t_username = sprintf( '%s,', lic.username{:} );
        t_username(end) = [];
    end
    if isempty( lic.hostid ),
        t_hostid = '';
    elseif ischar( lic.hostid ),
        t_hostid = lic.hostid;
    else
        t_hostid = sprintf( '%s,', lic.hostid{:} );
        t_hostid(end) = [];
    end
    message = sprintf( '%s|', lic.prefix, lic.name, lic.organization, lic.email, lic.license_type, t_username, t_hostid, lic.expiration, lic.prefix(end:-1:1) );
    dsa = java.security.Signature.getInstance('SHA1withDSA');
    dsa.initVerify(get_public_key);
    dsa.update(unicode2native(message,'UTF-8'));
    lic.status = {};
    if ~dsa.verify(int8(signature)),
        lic.status{end+1} = 'SIGNATURE';
    end
    if ~isempty( lic.hostid ) && ~any( cellfun( @(x)any(strcmp(x,lic.hostid)), get_hostid ) ) && ~any( strncmp( lic.hostid, '*', 1 ) ),
        lic.status{end+1} = 'HOSTID';
    end
    if ~isempty( lic.username ) && ~any( strcmpi( get_username, lic.username ) ),
        lic.status{end+1} = 'USER';
    end
    if days_left < 0,
        lic.status{end+1} = 'EXPIRED';
    end
    if isempty( lic.status ),
        lic.signature = signature;
        lic.status = 'VERIFIED';
    else
        lic.status = sprintf( '%s,', lic.status{:} );
        lic.status = sprintf( 'INVALID:%s', lic.status(1:end-1) );
    end
    lic.days_left = days_left;
catch exc 
    if ~isstruct( lic ) || numel( lic ) ~= 1, lic = []; end
    lic.status = my_get_report(exc);
    lic.days_left = -Inf;
end

function username = get_username
persistent p_username
if isempty( p_username )
    p_username = char(java.lang.System.getProperty('user.name'));
end
username = p_username;

function [ hostid_addr, hostid_name ] = get_hostid
persistent p_hostid_name p_hostid_addr
if isempty( p_hostid_addr )
    hostid_name = {}; 
    hostid_addr = {};
    networks = java.net.NetworkInterface.getNetworkInterfaces();
    while networks.hasMoreElements(),
        ni = networks.nextElement();
        try
            hostid = ni.getHardwareAddress();
        catch %#ok
            hostid = [];
        end
        if ~isempty(hostid),
            hostid_name{end+1} = char(ni.getName); %#ok
            hostid_addr{end+1} = sprintf('%02x',rem(double(hostid)+256,256)); %#ok
        end
    end
    if isempty( hostid_addr ),
        if strncmp( computer, 'MAC', 3 ),
            [status,str] = system('/sbin/ifconfig'); %#ok
            str = regexp( str, '^en\d:(\s+.*\n)*', 'match', 'lineanchors', 'dotexceptnewline' );
            for k = 1 : length(str),
                str2 = regexp( str{k}, '([\w\d]+):.*\sether\s([0-9a-f:]+)',  'tokens' );
                if ~isempty( str2 ),
                    hostid_name{end+1} = str2{1}{1};
                    hostid_addr{end+1} = strrep( str2{1}{2}, ':', '' );
                end
            end
        end
    end
    if ~isempty( hostid_name )
        if strncmp( computer, 'MAC', 3 ), 
            master = 'en'; 
        else
            master = 'eth'; 
        end
        ndxs = find( strncmp( hostid_name, master, length(master) ) );
        if ~isempty( ndxs ),
            hostid_name = hostid_name(ndxs); 
            hostid_addr = hostid_addr(ndxs);
        end
        [ hostid_name, ndxs2 ] = sort( hostid_name );
        hostid_addr = hostid_addr(ndxs2);
        if isempty( ndxs )
            % If this computer does not have any 'en*' or 'eth*' ports, we
            % are only going to trust the first hostID we find.
            hostid_name = hostid_name(1);
            hostid_addr = hostid_addr(1);
        end
    end
    p_hostid_name = hostid_name;
    p_hostid_addr = hostid_addr;
else
    hostid_name = p_hostid_name;
    hostid_addr = p_hostid_addr;
end

function lines = my_get_report( exc, debug )
try
    errmsg = getReport( exc, 'extended', 'hyperlinks', 'off' );
    errmsg = regexprep( errmsg,'</?a[^>]*>', '' );
catch %#ok
    errmsg = sprintf( '%s\n    Line %d: %s\n', exc.message, exc.stack(1).line, exc.stack(1).file );
end
width = 64;
lines = { 'UNEXPECTED ERROR: ---------------------------------------------' };
rndx = [ 0, regexp( errmsg, '\n' ), length(errmsg) + 1 ];
for k = 1 : length(rndx) - 1,
    line = errmsg( rndx(k)+1 : rndx(k+1) - 1 );
    if ~isempty( line ),
        emax     = length( line );
        n_indent = 0;
        if emax > width,
            f_indent = sum( regexp( line, '[^ ]', 'once' ) - 1 );
            sndxs = find( line == ' ' );
        end
        while true,
            if emax + n_indent <= width || isempty( sndxs ),
                lines{end+1} = [ 32 * ones(1,n_indent), line ]; %#ok
                break;
            end
            sndx = sndxs( sndxs <= width - n_indent + 1 );
            if isempty( sndx ), sndx = sndxs(1); end
            chunk = line(1:sndx(end)-1);
            lines{end+1} = [ 32*ones(1,n_indent), chunk ]; %#ok
            line(1:sndx(end)) = [];
            sndxs = sndxs(length(sndx)+1:end) - sndx(end);
            emax = emax - sndx(end);
            n_indent = f_indent + 4;
        end
    end
end
if nargin < 2 || ~debug,
    lines{end+1} = 'Please report this error to CVX Support by visiting';
    lines{end+1} = '    http://support.cvxr.com/support/tickets/new';
    lines{end+1} = 'or by sending an email to cvx@cvxr.com. Please include the full';
    lines{end+1} = 'output of this function in your report. Thank you!';
end
lines{end+1} = '---------------------------------------------------------------';
if nargin < 2 || ~debug,
    lines = sprintf( '%s\n', lines{:} );
end

%%%%%%%%%%%%%%
% END COMMON %
%%%%%%%%%%%%%%

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
