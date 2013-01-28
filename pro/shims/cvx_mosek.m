function shim = cvx_mosek( shim )

global cvx___
if ~isempty( shim.solve ),
    return
end
in_setup = isempty( shim.name );
shim.name = 'Mosek';
shim.dualize = true;
shim.check = [];
shim.solve = [];
if ~usejava('jvm'),
    shim.error = 'Java support is required.';
elseif isempty( cvx___.license ),
    shim.error = 'A CVX Professional license is required.';
else
    shim.error = 'An error occurred verifying the CVX Professional license.';
    try
        cvx___.license = full_verify( cvx___.license );
        if cvx___.license.days_left >= 0, 
            switch cvx___.license.license_type,
            case { 'academic', 'trial' },
                try_internal = true;
            otherwise,
                try_internal = any( strfind( cvx___.license.license_type, '+mosek' ) );
            end
            if try_internal,
                hostid = cvx___.license.hostid;
                if isempty( hostid ),
                    try_internal = false;
                elseif any( strcmp( hostid, '**' ) ),
                    try_internal = any( cellfun( @(x)any(strcmp(x,hostid)), get_hostid ) );
                end
            end
            shim.error = '';
        end
    catch exc
        shim.error = my_get_report( exc );
    end
end
if ~isempty( shim.error ),
    return
end
[ fs, ps, int_path, mext ] = cvx_version;
fname = [ 'mosekopt.', mext ];
int_plen = length( int_path );
if in_setup,
    flen = length(fname);
    fpaths = { [ int_path, fs, 'mosek', fs, mext(4:end), fs, fname ] };
    fpaths = [ fpaths ; which( fname, '-all' ) ];
    mosek_dir = getenv( 'MOSEKLM_LICENSE_FILE' );
    k = strfind( mosek_dir, fs );
    if length( k ) > 1,
        fpaths{end+1} = [ mosek_dir(1:k(end-1)), 'toolbox', fs, 'r2009b', fs, fname ];
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        is_internal = strncmp( new_dir, int_path, int_plen );
        if is_internal,
            tshim.location = [ '{cvx}', new_dir(int_plen+1:end-length(mext)+2) ];
        else
            tt = strfind( new_dir, fs );
            tshim.location = new_dir(1:tt(end-1)-1);
        end
        try
            otp = evalc('mosekopt');
        catch errmsg
            tshim.error = sprintf( 'Unexpected MEX file failure:\n%s\n', errmsg.message );
            clear mosekopt
        end
        if isempty( tshim.error ),
            otp = regexp( otp, 'MOSEK Version \S+', 'match' );
            if ~isempty(otp),
                tshim.version = otp{1}(15:end);
                sdp = sum(sscanf(tshim.version,'%d')) >= 7;
            else
                sdp = false;
            end
            if is_internal,
                if ~try_internal,
                    tshim.error = 'This license does not include the internal MOSEK solver.';
                else
                    mfunc = @m7;
                end
            elseif ~strcmp( cvx___.license.email, 'mcg@cvxr.com' ),
                mfunc = @mosekopt;
            elseif sdp,
                mfunc = @m7;
            else
                mfunc = @m6;
            end
        end
        if isempty( tshim.error ),
            try
                [rr,res] = mfunc('minimize echo(0)',struct('c',1,'a',sparse(1,1,1),'blc',0)); %#ok
                if res.rcode ~= 0,
                    tshim.error = sprintf( 'Error %s:\n%s\n', res.rcodestr, res.rmsg );
                end
            catch errmsg 
                tshim.error = sprintf( 'Unexpected MEX file failure:\n%s\n', errmsg.message );
            end
        end
        clear mosekopt
        if isempty( tshim.error ),
            tshim.check = @(varargin)check(sdp,varargin{:});
            tshim.solve = @(varargin)solve(sdp,mfunc,varargin{:});
            if k ~= 2,
                tshim.path = [ new_dir, ps ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'Could not find a MOSEK MEX file.';
    end
else
    shim.check = [];
    shim.solve = [];
    if ~isfield( shim, 'fullpath' ) || isempty( shim.fullpath ),
        if isempty( shim.path ),
            shim.path = which( fname );
            if isempty( shim.path ),
                shim.error = 'The MOSEK MEX file is missing from the MATLAB path. Please re-run CVX_SETUP.';
            end
        else
            shim.fullpath = [ shim(k).path(1:end-1), fname ];
        end
    elseif isempty( shim.path ) ~= strcmp( shim.fullpath, which( fname ) ),
        if isempty( shim.path ),
            temp = strfind( shim.fullpath, fs );
            shim.path = [ shim.fullpath(1:temp(end)), ps ];
        else
            shim.path = '';
        end
    end
    if strncmp( shim.fullpath, int_path, int_plen ),
        if ~try_internal,
            shim.error = 'This license does not include the internal MOSEK solver.';
        end
        mfunc = @m7;
        sdp = true;
    elseif ~strcmp( cvx___.license.email, 'mcg@cvxr.com' ),
        mfunc = @mosekopt;
        sdp = any( strfind( shim.fullpath, '/7/' ) );
    elseif any( strfind( shim.fullpath, '/6/' ) ),
        mfunc = @m6;
        sdp = false;
    else
        mfunc = @m7;
        sdp = true;
    end
    if isempty( shim.error ),
        shim.check = @(varargin)check(sdp,varargin{:});
        shim.solve = @(varargin)solve(sdp,mfunc,varargin{:});
    end
end

function [ rr, res ] = m6( command, varargin )
temp = [ 9, 4, 889, 3, 2013, 159, 82, 212, 183, 32, 156, 55, 5, 250, 40, 178, 78, 246, 213, 76, 76 ];
[ rr, res ] = mosekopt( [ command, ' lic' ], varargin{:}, temp );

function [ rr, res ] = m7( command, varargin )
temp = [ 9, 4, 889, 3, 2013, 40, 74, 95, 36, 238, 110, 221, 213, 152, 251, 223, 3, 130, 183, 92, 60 ];
[ rr, res ] = mosekopt( [ command, ' lic' ], varargin{:}, temp );

function found_bad = check( sdp, nonls )
found_bad = false;
if ~sdp,
    for k = 1 : length( nonls ),
        if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4,
            warning( 'CVX:Mosek:Semidefinite', 'This nonlinearity requires use of semidefinite cones which Mosek does not support.\n%s', ...
                'You will need to use a different solver for this model.' );
            found_bad = true;
        end
    end
end

function [ x, status, tol, iters, y, z ] = solve( sdp, mfunc, At, b, c, nonls, quiet, prec, settings )
zp = zeros(0,1);
[n,m] = size(At);
b = b(:); c = c(:);
prob  = struct( 'a', zp, 'blc', b, 'buc', b, 'blx', -Inf(n,1), 'bux', Inf(n,1), 'c', c, 'cones', {{}} );
prob.ints.sub = zp;
xscale = zp;
sdp_n = 0;
for k = 1 : length( nonls ),
    nonl = nonls(k);
    ti = nonl.indices;
    nn = size( ti, 1 );
    nv = size( ti, 2 );
    need_sdp = false;
    switch ( nonl.type ),
    case 'i_integer',
            
        prob.ints.sub = [ prob.ints.sub ; ti(:) ]; 
        
    case'i_binary',
        
        prob.ints.sub = [ prob.ints.sub ; ti(:) ]; 
        prob.blx( ti ) = 0; prob.bux( ti ) = 1;
        
    case 'i_semicontinuous',
        
        error( 'CVX:SolverIncompatible', 'MOSEK does not support semicontinous variables.' );
        
    case 'i_semiinteger',
        
        error( 'CVX:SolverIncompatible', 'MOSEK does not support semiinteger variables.' );
        
    case 'nonnegative',
        
        prob.blx( ti ) = 0;

    case 'lorentz',

        if nn == 1,
            prob.blx( ti ) = 0;
        else
            ti = ti([end,1:end-1],:);
            for qq = 1 : nv,
                prob.cones{end+1}.type = 'MSK_CT_QUAD';
                prob.cones{end}.sub = ti(:,qq)';
            end
        end

    case 'rotated-lorentz',

        if nn <= 2,
            prob.blx( ti ) = 0;
        else
            ti = ti([end-1:end,1:end-2],:);
            for qq = 1 : nv,
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub = ti(:,qq)';
            end
        end

    case 'semidefinite',

        if nn == 3,
            ti = ti([1,3,2],:);
            xscale(ti(1:2,:)) = true;
            for qq = 1 : nv,
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub  = ti(:,qq)';
            end
        else
            n2  = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
            qq2 = ( 1 : nn )';
            qqq = qq2;
            col = ceil( ( n2 + 0.5 - 0.5 / n2 ) - sqrt( ( n2 + 0.5 )^2 - 2 * qqq ) );
            row = qqq - ( col - 1 ) .* ( 2 * n2 - col ) / 2;
            vv2 = 2 * ones(nn,1);
            vv2(cumsum([1,n2:-1:2]),1) = 1;
            vv1 = 1.0 ./ vv2;
            need_sdp = true;
        end

    case 'hermitian-semidefinite',

        if nn == 4,
            ti = ti([1,4,2,3],:);
            xscale(ti(1:2,:)) = true;
            for qq = 1 : nv,
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub  = ti(:,qq)';
            end
        else
            %   X >= 0 <==> exists [ Y1, Y2^T ; Y2, Y3 ] >= 0 s.t.
            %               Y1 + Y3 == real(X), Y2 - Y2^T == imag(X)
            % So: <C,X> = <CR,XR> + <CI,XI>
            %           = <CR,Y1+Y3> + <CI,Y2-Y2^T>
            %           = <CR,Y1> + <CI,Y2> + <CI^T,Y2^T> + <CR,Y3>
            %           = < [ CR, CI^T ; CI, CR ], [ Y1, Y2^T ; Y2, Y3 ] >
            n2   = sqrt( nn );
            qqq  = ( 1 : ( 2 * nn + n2 ) )';
            col  = ceil( ( 2 * n2 + 0.5 - 0.25 / n2 ) - sqrt( ( 2 * n2 + 0.5 )^2 - 2 * qqq ) );
            row  = qqq - ( col - 1 ) .* ( 4 * n2 - col ) / 2;
            row2 = rem( row - 1, n2 ) + 1;
            col2 = rem( col - 1, n2 ) + 1;
            dig  = row2 == col2;
            img  = row > n2 & col <= n2;
            neg  = img & row2 < col2;
            tmp  = col2(neg); col2(neg) = row2(neg); row2(neg) = tmp;
            qq2  = 2 * ( row2 - 1 ) + ( img | dig ) + ( 2 * n2 - 1 ) * ( col2 - 1 ) - col2 .* ( col2 - 1 );
            vv2  = ( 2 - dig ) .* ( 1 - 2 * neg );
            tt   = cumsum([n2+1,2*n2:-1:n2+2]);
            qqq(tt) = []; qq2(tt) = []; vv2(tt) = [];
            vv1  = 1.0 ./ vv2;
            vv2  = 0.5 * vv2;
            n2   = 2 * n2;
            nn   = n2 * ( n2 + 1 ) / 2;
            need_sdp = true;
        end

    case 'exponential',

        error( 'CVX:SolverIncompatible', 'MOSEK does not support the exponential cone.\nYou must use another solver for this problem.' );

    otherwise,

        error( 'Invalid cone type: %s', tt );

    end    
    if need_sdp,
        
        if ~sdp,
            error( 'CVX:SolverIncompatible', 'This version of MOSEK does not support semidefinite cones larger than 2x2.' );
        end
        
        if ~sdp_n,
            prob.bardim = zp;
            prob.barc = struct( 'subj', zp, 'subk', zp, 'subl', zp, 'val', zp );
            prob.bara = struct( 'subi', zp, 'subj', zp, 'subk', zp, 'subl', zp, 'val', zp );
            sndxi = zp; sndxj = zp; sndxp = zp; sndxd = zp;
        end

        if nv > 1,
            qqq = bsxfun( @plus, qqq(:), qqq(end) * ( 0 : nv - 1 ) );
            qq2 = bsxfun( @plus, qq2(:), qq2(end) * ( 0 : nv - 1 ) );
            vv2 = bsxfun( @plus, vv2(:), zeros(1,nv) );
            vv1 = bsxfun( @plus, vv1(:), zeros(1,nv) );
        end
        qq2 = ti(qq2);
        F = sparse( qqq, qq2, vv1, qqq(end), length(c) );
        prob.bardim = [ prob.bardim ; n2(ones(1,nv),:) ];
        
        cc = F * c;
        if nnz( cc ),
            [ rr, cc, vv ] = find( cc ); %#ok
            mat = floor((rr-1)/nn);
            rr  = rr - nn * mat;
            prob.barc.subj = [ prob.barc.subj ; mat + 1 ];
            prob.barc.subk = [ prob.barc.subk ; row(rr) ];
            prob.barc.subl = [ prob.barc.subl ; col(rr) ];
            prob.barc.val  = [ prob.barc.val  ; vv ];
        end

        cc = F * At;
        if nnz( cc ),
            [ rr, cc, vv ] = find( cc );
            mat = floor((rr-1)/nn);
            rr  = rr - nn * mat;
            prob.bara.subi = [ prob.bara.subi ; cc  ];
            prob.bara.subj = [ prob.bara.subj ; mat + 1 ];
            prob.bara.subk = [ prob.bara.subk ; row(rr) ];
            prob.bara.subl = [ prob.bara.subl ; col(rr) ];
            prob.bara.val  = [ prob.bara.val  ; vv ];
        end
        
        sndxi = [ sndxi ; qq2(:) ]; %#ok
        sndxj = [ sndxj ; qqq(:) ]; %#ok
        sndxp = [ sndxp ; sign(vv2(:)) ]; %#ok
        sndxd = [ sndxd ; vv2(:) ]; %#ok
        sdp_n = sdp_n + nn * nv;
        
    end
end
if ~isempty( xscale ),
    alpha = sqrt(2.0);
    xscale(end+1:n) = true;
    prob.c(xscale) = prob.c(xscale) * alpha;
    At(xscale,:) = At(xscale,:) * alpha;
end
if sdp_n,
    zndxs = 1 : n;
    zndxs(sndxi) = [];
    qndxs = zeros(1,n);
    qndxs(zndxs) = 1 : numel(zndxs);
    prob.c = prob.c(zndxs);
    prob.a = At(zndxs,:).';
    prob.blx = prob.blx(zndxs);
    prob.bux = prob.bux(zndxs);
    prob.ints.sub = qndxs(prob.ints.sub);
    for k = 1 : length( prob.cones ),
        prob.cones{k}.sub = qndxs(prob.cones{k}.sub);
    end
else
    prob.a = At.';
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
[ rr, res ] = cvx_run_solver( mfunc, command, prob, param, 'rr', 'res', settings, 3 ); %#ok
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
if sdp_n,
    zndxs = sparse( zndxs, 1:numel(zndxs), 1, n, numel(zndxs) );
    sndxp = sparse( sndxi, sndxj, sndxp, n, sdp_n );
    x = full( zndxs * sol.xx + sndxp * sol.barx );
else
    x = sol.xx;
end
if has_dual,
    y = sol.slc - sol.suc;
    z = sol.slx - sol.sux;
    if isfield( sol, 'snx' ), 
        z = z + sol.snx; 
    end
    if sdp_n,
        sndxd = sparse( sndxi, sndxj, sndxd, n, sdp_n );
        z = full( zndxs * z + sndxd * sol.bars );
    end
else
    y = NaN * ones(m,1);
    z = NaN * ones(n,1);
end
if ~isempty( xscale ),
    x(xscale) = x(xscale) * alpha;
    z(xscale) = z(xscale) / alpha;
end
status = '';
switch sol.solsta,
    case { 'NEAR_PRIMAL_INFEASIBLE_CER', 'PRIMAL_INFEASIBLE_CER' },
        status = 'Infeasible';
        x(:) = NaN; scl = abs(b'*y); z = z / scl; y = y / scl;
    case { 'NEAR_DUAL_INFEASIBLE_CER', 'DUAL_INFEASIBLE_CER' },
        status = 'Unbounded';
        y(:) = NaN; z(:) = NaN; x = x / abs(c'*x);
    case { 'OPTIMAL', 'NEAR_OPTIMAL' },
        status = 'Solved';
        if sol.solsta(1) == 'N', tol = prec(3); end
    case 'UNKNOWN',
        if has_dual,
            pobj = c' * x;
            dobj = b' * y;
            nrmc = norm( c );
            nrmb = norm( b );
            xc   = At' * x;
            zc   = At * y + z;
            if c' * x < 0,
                uerr = -nrmc * norm(xc) / max( nrmb, 1 ) / pobj;
            else
                uerr = Inf;
            end
            if b' * y < 0,
                ferr = -nrmb * norm(zc) / max( nrmc, 1 ) / dobj;
            else
                ferr = Inf;
            end
            perr = norm( xc - b ) /  ( 1 + nrmb );
            derr = norm( zc - c ) / ( 1 + norm( c ) );
            gerr = max( 0, pobj - dobj ) / max( 1, abs(pobj) );
            tol2 = max( [ perr, derr, gerr ] );
            tol  = min( [ tol2, ferr, uerr ] );
            if tol == tol2,
                status = 'Solved';
            elseif tol == ferr,
                status = 'Infeasible';
                x(:) = NaN; z = z / abs(dobj); y = y / abs(dobj);
            else
                status = 'Unbounded';
                y(:) = NaN; z(:) = NaN; x = x / abs(pobj);
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
            y(:) = NaN; z(:) = NaN; x = x / abs(c'*x);
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
    if ~dsa.verify(int8(signature)),
        lic.status = 'INVALID:SIGNATURE';
    elseif ~isempty( lic.hostid ) && ~any( cellfun( @(x)any(strcmp(x,lic.hostid)), get_hostid ) ) && ~any( strcmp(lic.hostid,'*') ),
        lic.status = 'INVALID:HOSTID';
    elseif ~isempty( lic.username ) && ~any( strcmpi( get_username, lic.username ) ),
        lic.status = 'INVALID:USER';
    elseif days_left < 0,
        lic.status = 'EXPIRED';
    else
        lic.signature = signature;
        lic.status = 'VERIFIED';
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

function estr = my_get_report( exc )
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
lines{end+1} = 'Please report this error to CVX Support by visiting';
lines{end+1} = '    http://support.cvxr.com/support/tickets/new';
lines{end+1} = 'or by sending an email to cvx@cvxr.com. Please include the full';
lines{end+1} = 'output of this function in your report. Thank you!';
lines{end+1} = '---------------------------------------------------------------';
estr = sprintf( '%s\n', lines{:} );

%%%%%%%%%%%%%%
% END COMMON %
%%%%%%%%%%%%%%

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
