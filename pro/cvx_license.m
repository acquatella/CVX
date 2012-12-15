function [ license, ltext ] = cvx_license( varargin )

% CVX_LICENSE   License processing for CVX Professional.
%    This file performs various functions needed to perform license
%    management for the professional features of CVX.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick exit for cvx_global %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~usejava( 'jvm' ),
    error( 'CVX:Licensing', 'The CVX licensing mechanism requires Java.' );
elseif nargin == 1 && isstruct( varargin{1} ) && numel( varargin{1} ) == 1,
    license = full_verify( varargin{1} );
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full processing for cvx_setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent blank fs
if isempty( fs ),
    blank = struct( 'name', '', 'organization', '', 'email', '', ...
        'license_type', '', 'username', {{}}, 'hostid', {{}}, ...
        'expiration', '0000-00-00', 'signature', int8([]), 'prefix', '', ...
        'status', 'NOTFOUND', 'days_left', -Inf, 'filename', '' );
    if strncmp( computer, 'PC', 2 ), 
        fs = '\'; 
    else
        fs = '/'; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input arguments %
%%%%%%%%%%%%%%%%%%%%%%%%%

clear_all = false;
lnames = { '' };
verbose = true;
for k = 1 : nargin,
    arg = varargin{k};
    if isempty( arg ),
        continue;
    elseif ischar( arg ),
        switch arg,
            case '*key*',
                license = get_public_key;
                return
            case '-quiet',
                verbose = false;
            case '-clear',
                clear_all = true;
            otherwise,
                if exist( arg, 'file' ),
                    lnames{end+1} = arg; %#ok
                elseif exist( [ arg, '.mat' ], 'file' ),
                    lnames{end+1} = [ arg, '.mat' ], %#ok
                else
                    error( 'CVX:Licensing', 'License file "%s" does not exist.', arg );
                end
        end
    else
        error( 'CVX:Licensing', 'Invalid use of the CVX licensing system.' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for all available license files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build search order:
% --- CVX saved preferences
% --- Matlab path
% --- Home directory
% --- Desktop directory
if length( lnames ) <= 1,
    dname = 'cvx_license.dat';
    lnames = [ lnames ; which( dname, '-all' ) ];
    homedir = char(java.lang.System.getProperty('user.home'));
    lnames{end+1} = [ homedir, fs, dname ];
    lnames{end+1} = [ homedir, fs, 'Desktop', fs, dname ];
    lnames{end+1} = [ homedir, fs, 'Downloads', fs, dname ];
    [ dummy, ndxs ] = sort( lnames ); dummy = dummy(:);
    tt = [ true ; ~strcmp( dummy(2:end), dummy(1:end-1) ) ];
    lnames = lnames( sort(ndxs(tt)) );
end

%%%%%%%%%%%%%%%%%%%%%
% Test each license %
%%%%%%%%%%%%%%%%%%%%%

best_days = -1;
best_ndx = 0;
licenses = blank;
found_saved = false;
for k = 1 : length(lnames),
    lic = load_and_verify( lnames{k}, blank, fs );
    if isequal( lic.status, 'NOTFOUND' ), continue; end
    if k == 1, 
        found_saved = true; 
    end
    found = false;
    for kk = 2 : length(licenses),
        if ( isequal( licenses(kk).expiration, lic.expiration ) && ...
             isequal( licenses(kk).signature, lic.signature ) ),
            licenses(kk).filename{end+1} = lic.filename; %#ok
            found = true;
            break
        end
    end
    if ~found,
        lic.filename = { lic.filename };
        try
            licenses(end+1) = lic; %#ok
        catch %#ok
            lic2 = blank;
            for k = fieldnames(lic2)',
                try lic2.(k{1}) = lic.(k{1}); catch end %#ok
            end
            licenses(end+1) = lic2; %#ok
        end
        if isequal( lic.status, 'VERIFIED' ) && lic.days_left > best_days,
            best_days = lic.days_left;
            best_ndx = length(licenses);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print current host info and licenses found %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose,
    ltext = {};
    if nargout == 0,
        ltext{end+1} = '';
    end
    ltext{end+1} = 'License host:';
    ltext{end+1} = sprintf( '    Username: %s', get_username );
    [ hostid_addr, hostid_name ] = get_hostid;
    if isempty( hostid_addr )
        ltext{end+1} = '    Host IDs: none';
    else
        ltext{end+1} = sprintf( '    Host IDs: %s (%s)', hostid_addr{1}, hostid_name{1} );
    end
    ndxs = false(1,length(licenses));
    ndxs(1) = true;
    ltext{end+1} = 'Installed license:';
    if found_saved,
        ltext = print_license(licenses(2),'    ',ltext);
        ndxs(2) = true;
    else
        ltext{end+1} = '    No license installed.';
    end
    if best_days >= 0 && ~ndxs(best_ndx),
        if found_saved,
            ltext{end+1} = 'Replacement license:';
        elseif nargout > 0,
            ltext{end+1} = 'Installing license:';
        else
            ltext{end+1} = 'Valid license found (run "cvx_setup" to install):';
        end
        ndxs(best_ndx) = true;
        ltext = print_license(licenses(best_ndx),'    ',ltext);
    end
    if any(~ndxs),
        if any(ndxs),
            prefix = '    ';
            if nnz(ndxs) > 1,
                ltext{end+1} = 'Other licenses found:';
                prefix2 = '        ';
            else
                ltext{end+1} = 'Other license found:';
                prefix2 = '    ';
            end
        else
            prefix = '';
            prefix2 = '    ';
        end
        for k = 1 : length(ndxs),
            if ~ndxs(k),
                ltext = print_license(licenses(k),prefix,ltext,prefix2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print all licenses found %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if best_days < 0,
    if verbose,
        ltext{end+1} = 'No valid licenses found.';
    end
    license = [];
elseif clear_all,
    if verbose,
        ltext{end+1} = 'License clear requested.';
    end
    license = [];
else
    license = licenses( best_ndx );
    license.filename = license.filename{1};
end
if nargout == 0,
    clear license
    if verbose,
        ltext{end+1} = '';
    end
end
if verbose,
    fprintf( '%s\n', ltext{:} );
end

function ltext = print_license( lic, prefix, ltext, prefix2 )
if isfield( lic, 'filename' ) && ~isempty( lic.filename ),
    fprefix = '%sFile: %s';
    fprefix2 = '%sAlso in file: %s';
    fnames = lic.filename;
    if ~iscell( fnames ),
        fnames = { fnames };
    end
    if ~isempty( fnames ) && ~isempty( fnames{1} ) && fnames{1}(1) == '(',
        fnames(1) = [];
        fprefix = fprefix2;
    end
    for k = 1 : length(fnames),
        ltext{end+1} = sprintf( fprefix, prefix, fnames{k} ); %#ok
        if nargin > 3, prefix = prefix2; end
        fprefix = fprefix2;
    end
end
if any( strcmp( lic.status, { 'EMPTY', 'CORRUPT' } ) )
    ltext{end+1} = sprintf( '%sStatus: %s', prefix, lic.status );
    return
end
if ~isempty( lic.organization ),
    ltext{end+1} = sprintf( '%sOrganization: %s', prefix, lic.organization );
end
if ~isempty( lic.name ),
    if isempty( lic.email ),
        ltext{end+1} = sprintf( '%sContact: %s', prefix, lic.name );
    else
        ltext{end+1} = sprintf( '%sContact: %s (%s)', prefix, lic.name, lic.email );
    end
elseif ~isempty( lic.email ),
    ltext{end+1} = sprintf( '%sContact: %s', prefix, lic.email );
end
if ~isempty( lic.license_type ),
    ltext{end+1} = sprintf( '%sLicense type: %s', prefix, lic.license_type );
end
if ~isempty( lic.username ),
    if any( strcmpi( get_username, lic.username ) ),
        status = '';
    else
        status = ' (MISMATCH)';
    end
    if ischar( lic.username ),
        l_username = lic.username;
    else
        l_username = sprintf( '%s, ', lic.username{:} );
        l_username = l_username(1:end-2);
    end
    ltext{end+1} = sprintf( '%sNamed user: %s%s', prefix, l_username, status );
end
if ~isempty( lic.hostid ),
    [ hostid_addr, hostid_name ] = get_hostid;
    ndxs = cellfun( @(x)any(strcmp(x,lic.hostid)), hostid_addr );
    if any( ndxs ),
        hostid_name = hostid_name(ndxs);
        hostid_name = hostid_name(~cellfun('isempty',hostid_name));
        if isempty( hostid_name ),
            status = '';
        else
            status = sprintf( '%s,', hostid_name{ndxs} );
            status = sprintf( ' (%s)', status(1:end-1) );
        end
    elseif isempty( hostid_addr ),
        status = ' (MISMATCH: no host id)';
    else
        status = ' (MISMATCH)';
    end
    if ischar( lic.hostid ),
        l_hostid = lic.hostid;
    else
        l_hostid = sprintf( '%s, ', lic.hostid{:} );
        l_hostid = l_hostid(1:end-2);
    end
    ltext{end+1} = sprintf( '%sHost ID: %s%s', prefix, l_hostid, status );
end
parser = java.text.SimpleDateFormat('yyyy-MM-dd');
try
    expire = parser.parse(lic.expiration);
catch %#ok
    lic.expiration = '0000-01-01';
    expire = parser.parse(lic.expiration);
end
if lic.days_left < -365
    ltext{end+1} = sprintf( '%sEXPIRED: %s', prefix, lic.expiration );
elseif lic.days_left < 0
    if lic.days_left == -1, plural = ''; else plural = 's'; end
    ltext{end+1} = sprintf( '%sEXPIRED: %s (%d day%s ago)', prefix, lic.expiration, -lic.days_left, plural );
else
    lic.days_left = abs(lic.days_left);
    if lic.days_left == 1, plural = ''; else plural = 's'; end
    ltext{end+1} = sprintf( '%sExpiration: %s (%g day%s remaining)', prefix, lic.expiration, lic.days_left, plural );
end
if isequal( lic.status, 'INVALID:SIGNATURE' ) || isequal( lic.status, 'CORRUPT:SIGNATURE' ),
    ltext{end+1} = sprintf( '%sSignature: INVALID', prefix );
else
    ltext{end+1} = sprintf( '%sSignature: valid', prefix );
end
if strncmp( lic.status, 'UNEXPECTED ERROR:', 17 ),
    rndx = [ 0, regexp( lic.status, '\n' ), length(lic.status) + 1 ];
    for k = 1 : length(rndx)-1;
        ltext{end+1} = sprintf('%s%s',prefix,lic.status(rndx(k)+1:rndx(k+1)-1));
        spacer = '    ';
    end
end

function public_key = get_public_key
persistent p_public_key
if isempty( p_public_key ),
    keyfac = java.security.KeyFactory.getInstance('DSA');
    p_public_key = '308201b73082012c06072a8648ce3804013082011f02818100fd7f53811d75122952df4a9c2eece4e7f611b7523cef4400c31e3f80b6512669455d402251fb593d8d58fabfc5f5ba30f6cb9b556cd7813b801d346ff26660b76b9950a5a49f9fe8047b1022c24fbba9d7feb7c61bf83b57e7c6a8a6150f04fb83f6d3c51ec3023554135a169132f675f3ae2b61d72aeff22203199dd14801c70215009760508f15230bccb292b982a2eb840bf0581cf502818100f7e1a085d69b3ddecbbcab5c36b857b97994afbbfa3aea82f9574c0b3d0782675159578ebad4594fe67107108180b449167123e84c281613b7cf09328cc8a6e13c167a8b547c8d28e0a3ae1e2bb3a675916ea37f0bfa213562f1fb627a01243bcca4f1bea8519089a883dfe15ae59f06928b665e807b552564014c3bfecf492a038184000281804b5b19f08b0852bfd6527750cd897cf44909e5f875accee48fc43a1482c23d8203c63d63dfb353d4be3a3152ed2b15ed7cdb6ea3bfc5278ca3975d7f32f920bc6317417c0145afc8513b4920c75f0cec5eddbe325a8a311a2a5b4d84e5ba9675ba95c169e449637abc54baa1c95d87ae2cbd1d797bbbbfc63e00d2674f86d317';
    p_public_key = uint8(hex2dec(reshape(p_public_key,2,length(p_public_key)/2)'));
    p_public_key = java.security.spec.X509EncodedKeySpec(p_public_key);
    p_public_key = keyfac.generatePublic(p_public_key);
end
public_key = p_public_key;

function lic = load_and_verify( fname, blank, fs )
persistent base64
if isempty( base64 ),
    base64 = uint8(zeros(1,256)+65);
    base64(uint8(['A':'Z', 'a':'z', '0':'9', '+/=']))= 0:64;
    base64(uint8('-_'))= 62:63;
end
try
    lic = [];
    if isempty( fname ) && exist( [ prefdir, fs, 'cvx_prefs.mat' ], 'file' ),
        lic = load( [ prefdir, fs, 'cvx_prefs.mat' ], 'license' );
        lic = lic.license;
        fname = '(from saved preferences)';
    elseif exist( fname, 'file' ),
        fid = fopen( fname );
        lic = textscan(fid,'%s%s','delimiter','=');
        fclose( fid );
        lic{2}(1:end-1) = cellfun(@(x)native2unicode(x,'UTF-8'),lic{2}(1:end-1),'UniformOutput',false);
        lic = cell2struct( lic{2}, lic{1} );
        if isempty(lic.username),
            lic.username = {};
        else
            tmp = [0,strfind(lic.username,','),length(lic.username)+1];
            usernames = cell(1,length(tmp)-1);
            for k = 1 : length(usernames),
                usernames{k} = lic.username(tmp(k)+1:tmp(k+1)-1);
            end
            lic.username = usernames;
        end
        if isempty(lic.hostid),
            lic.hostid = {};
        else
            tmp = [0,strfind(lic.hostid,','),length(lic.hostid)+1];
            hostids = cell(1,length(tmp)-1);
            for k = 1 : length(hostids),
                hostids{k} = lic.hostid(tmp(k)+1:tmp(k+1)-1);
            end
            lic.hostid = hostids;
        end
        x = base64( lic.signature );
        nbytes = length(x);
        nchunks = ceil( nbytes / 4 );
        x( end + 1 : 4 * nchunks ) = 0;
        x = reshape( x, 4, nchunks );
        y = [ bitshift( x(1,:), 2 ) ; bitshift( x(2,:), 4 ) ; bitshift( x(3,:), 6 ) ]; 
        y = bitor( y, [ bitshift( x(2,:), -4 ) ; bitshift( x(3,:), -2 ) ; x(4,:) ] );
        y = y(:)';
        ndxs = strfind(y,'|');
        lic.prefix = char(y(1:ndxs(1)-1));
        y = int16(y(ndxs(1)+1:ndxs(end)-1));
        y(y>127) = y(y>127) - 256;
        lic.signature = int8(y);
    end
    if isempty( lic ), 
        lic = blank;
    else
        lic = full_verify( lic ); 
    end
    lic.filename = fname;
catch exc
    if ~isstruct( lic ) || numel( lic ) ~= 1, lic = blank; end
    lic.status = my_get_report(exc);
    lic.signature = [];
    lic.days_left = -Inf;
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
    elseif ~isempty( lic.hostid ) && ~any( cellfun( @(x)any(strcmp(x,lic.hostid)), get_hostid ) )
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
                lines{end+1} = [ 32 * ones(1,n_indent), line ];
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
% The command 'cvx_where' will show where this file is located.





