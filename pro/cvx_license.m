function [ license, ltext ] = cvx_license( seed_name )

% CVX_LICENSE   License processing for CVX Professional.
%    This file performs various functions needed to perform license
%    management for the professional features of CVX.

persistent public_key

if ~usejava( 'jvm' ),
    license = [];
    return
end

if isempty( public_key ),
    keyfac = java.security.KeyFactory.getInstance('DSA');
    public_key = '308201b73082012c06072a8648ce3804013082011f02818100fd7f53811d75122952df4a9c2eece4e7f611b7523cef4400c31e3f80b6512669455d402251fb593d8d58fabfc5f5ba30f6cb9b556cd7813b801d346ff26660b76b9950a5a49f9fe8047b1022c24fbba9d7feb7c61bf83b57e7c6a8a6150f04fb83f6d3c51ec3023554135a169132f675f3ae2b61d72aeff22203199dd14801c70215009760508f15230bccb292b982a2eb840bf0581cf502818100f7e1a085d69b3ddecbbcab5c36b857b97994afbbfa3aea82f9574c0b3d0782675159578ebad4594fe67107108180b449167123e84c281613b7cf09328cc8a6e13c167a8b547c8d28e0a3ae1e2bb3a675916ea37f0bfa213562f1fb627a01243bcca4f1bea8519089a883dfe15ae59f06928b665e807b552564014c3bfecf492a038184000281804b5b19f08b0852bfd6527750cd897cf44909e5f875accee48fc43a1482c23d8203c63d63dfb353d4be3a3152ed2b15ed7cdb6ea3bfc5278ca3975d7f32f920bc6317417c0145afc8513b4920c75f0cec5eddbe325a8a311a2a5b4d84e5ba9675ba95c169e449637abc54baa1c95d87ae2cbd1d797bbbbfc63e00d2674f86d317';
    public_key = uint8(hex2dec(reshape(public_key,2,length(public_key)/2)'));
    public_key = java.security.spec.X509EncodedKeySpec(public_key);
    public_key = keyfac.generatePublic(public_key);
    clear keyfac
end

if nargin == 1,
    if isstruct( seed_name ),
        license = quick_verify( seed_name, public_key );
        return
    elseif ~isempty( seed_name )
        if ~ischar( seed_name ) && size( seed_name, 1 ) > 1,
            error( 'CVX:Licensing', 'License filename must be a string.' );
        elseif strcmpi( seed_name, '*key*' )
            license = public_key;
            return
        elseif ~exist( seed_name, 'file' ) && ~exist( [ seed_name, '.mat' ], 'file' ),
            error( 'CVX:Licensing', 'License file %s does not exist.', seed_name );
        end
    end
end
if strncmp( computer, 'PC', 2 ), 
    fs = '\'; 
else
    fs = '/'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Check current license %
%%%%%%%%%%%%%%%%%%%%%%%%%

try
    cur_lic = load( [ prefdir, fs, 'cvx_prefs.mat' ], 'license' );
    cur_lic = cur_lic.license;
    if isempty( cur_lic.signature ),
        cur_lic = [];
    else
        cur_lic = full_verify( cur_lic, public_key );
        cur_lic.filename = { '(from saved preferences)' };
    end
catch %#ok
    cur_lic = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for all available license files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build search order:
% --- CVX saved preferences
% --- Matlab path
% --- Home directory
% --- Desktop directory
lnames = {};
if nargin == 1,
    if ~isempty( seed_name ),
        if ~exist( seed_name, 'file' ), 
            seed_name = [ seed_name, '.mat' ]; 
        end
        lnames = { seed_name };
    end
elseif nargin == 0,
    lname = 'cvx_license.mat';
    lnames = which( lname, '-all' );
    homedir = char(java.lang.System.getProperty('user.home'));
    lnames{end+1} = [ homedir, fs, lname ];
    lnames{end+1} = [ homedir, fs, 'Desktop', fs, lname ];
    [ dummy, ndxs ] = sort( lnames );
    tt = [ true, ~strcmp( dummy(2:end), dummy(1:end-1) ) ];
    lnames = lnames( sort(ndxs(tt)) );
else
    lnames = {};
end

%%%%%%%%%%%%%%%%%%%%%
% Test each license %
%%%%%%%%%%%%%%%%%%%%%

ndx = 0;
if ~isempty( cur_lic ),
    ndx = ndx + 1;
    licenses = cur_lic;
else
    licenses = [];
end
for k = 1 : length(lnames),
    fname = lnames{k};
    if ~exist( fname, 'file' ),
        continue;
    end
    lic = load( fname );
    if ( ~isempty( cur_lic ) && ...
         isequal( cur_lic.expiration, lic.expiration ) && ...
         isequal( cur_lic.signature, lic.signature ) ),
        cur_lic.filename{end+1} = fname;
        continue
    end
    for kk = 1 : ndx,
        if ( isequal( licenses(kk).expiration, lic.expiration ) && ...
             isequal( licenses(kk).signature, lic.signature ) ),
            licenses(kk).filename{end+1} = fname;
            continue
        end
    end
    lic.filename = { fname };
    lic.status = 'UNVERIFIED';
    lic = full_verify( lic , public_key );
    ndx = ndx + 1;
    if ndx == 1,
        licenses = lic;
    else
        try
            licenses(ndx) = lic;
        catch %#ok
            for kk = fieldnames(lic)',
                licenses(ndx).(kk{1}) = lic.(kk{1});
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print all licenses found %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ltext = {};
if nargout == 0 && nargin == 0,
    ltext{end+1} = '';
end
best_days = -Inf;
best_ndx = 0;
if ~isempty( cur_lic ),
    ltext{end+1} ='Current license:';
    ltext = print_license( cur_lic, '    ', ltext );
    if length(cur_lic.filename) > 1,
        ltext{end+1} = sprintf( '    Also in file: %s', cur_lic.filename{2:end} );
    end
    if cur_lic.days_left > best_days,
        best_days = cur_lic.days_left;
        best_ndx = 1;
    end
    nfound = 1;
    if nfound < length(licenses),
        ltext{end+1} = 'Other license files found:';
        for k = 2 : length(licenses),
            ltext{end+1} = sprintf( '    File: %s', licenses(k).filename{:} ); %#ok
            ltext = print_license( licenses(k), '        ', ltext );
            if licenses(k).days_left > best_days,
                best_days = licenses(k).days_left;
                best_ndx = k;
            end
        end
    end
elseif ~isempty( licenses ),
    if length(licenses) == 1,
        ltext{end+1} = 'License file found:';
    else
        ltext{end+1} = 'License files found:';
    end
    for k = 1 : length(licenses),
        ltext{end+1} = sprintf( '    File: %s', licenses(k).filename{:} ); %#ok
        ltext = print_license(licenses(k),'    ',ltext);
        if licenses(k).days_left > best_days,
            best_days = licenses(k).days_left;
            best_ndx = k;
        end
    end
end
if best_days < 0,
    ltext{end+1} = 'No valid licenses found.';
    if nargout == 0 || nargin == 0 || isempty( seed_name ),
        if nargout == 0 && nargin == 0,
            ltext{end+1} = '';
            ltext{end+1} = 'Licensing info:';
        end
        username = char(java.lang.System.getProperty('user.name'));
        ltext{end+1} = sprintf('    Named user: %s', username );
        networks = java.net.NetworkInterface.getNetworkInterfaces();
        if strncmp(computer,'MAC',3), prefix = 'en'; else prefix = 'eth'; end
        addrs = {}; names = {};
        while networks.hasMoreElements(),
            ni = networks.nextElement();
            addr = sprintf('%02x',rem(double(ni.getHardwareAddress())+256,256));
            if ~isempty( addr )
                name = char(ni.getName());
                if strncmp( prefix, name, length(prefix) ),
                    names{end+1} = name; %#ok
                    addrs{end+1} = addr; %#ok
                end
            end
        end
        if isempty(names),
            ltext{end+1} = '    Could not determine host ID.';
        else
            [ names, ndxs ] = sort( names );
            ltext{end+1} = sprintf('    Host ID: %s (%s)', addrs{ndxs(1)}, names{1} );
        end
    end
    license = [];
elseif nargin == 0,
    license = licenses(best_ndx);
elseif isempty( seed_name ),
    license = cur_lic;
elseif length( licenses ) < 1 + ~isempty(cur_lic),
    license = [];
else
    license = licenses(1+~isempty(cur_lic));
end
if ~isempty( license ) && iscell( license.filename ),
    license.filename = license.filename{1};
end
if nargout == 0,
    if nargin == 0,
        ltext{end+1} = '';
    end
    clear license
end
if nargout < 2,
    fprintf( '%s\n', ltext{:} );
    clear ltext
end

function ltext = print_license( lic, prefix, ltext )
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
if ~isempty( lic.username ),
    username = char(java.lang.System.getProperty('user.name'));
    if ~isequal( username, lic.username ), 
        ltext{end+1} = sprintf( '%sNamed user: %s (MISMATCH: %s)', prefix, lic.username, username );
    else
        ltext{end+1} = sprintf( '%sNamed user: %s', prefix, lic.username );
    end
end
if ~isempty( lic.hostid ),
    networks = java.net.NetworkInterface.getNetworkInterfaces();
    names = {}; addrs = {};
    master = '';
    while networks.hasMoreElements(),
        ni = networks.nextElement();
        hostid = sprintf('%02x',rem(double(ni.getHardwareAddress())+256,256));
        if isempty(hostid),
            continue;
        elseif isequal( hostid, lic.hostid ),
            master = char(ni.getName);
            break;
        else
            names{end+1} = char(ni.getName); %#ok
            addrs{end+1} = hostid; %#ok
        end
    end 
    if ~isempty( master ),
        ltext{end+1} = sprintf( '%sHost ID: %s (%s)', prefix, lic.hostid, master );
    elseif  isempty(addrs),
        ltext{end+1} = sprintf( '%sHost ID: %s (MISMATCH: no host id)', prefix, lic.hostid );
    else
        if strncmp( computer, 'MAC', 3 ),
            master = 'en0'; 
        else
            master = 'eth0'; 
        end
        ndxs = find(strcmp(names,master));
        if isempty(ndxs),
            [dum,ndxs] = sort(names); %#ok
            ndxs = ndxs(1);
        end
        ltext{end+1} = sprintf( '%sHost ID: %s (MISMATCH: %s (%s))', prefix, lic.hostid, addrs{ndxs}, names{ndxs} );
    end
end
parser = java.text.SimpleDateFormat('yyyy-MM-dd');
try
    expire = parser.parse(lic.expiration);
catch %#ok
    lic.expiration = '0000-01-01';
    expire = parser.parse(lic.expiration);
end
today = java.util.Date;
remain_days = ceil( ( double(expire.getTime()) - double(today.getTime) ) / 86400000 );
if remain_days < -365
    ltext{end+1} = sprintf( '%sEXPIRED: %s', prefix, lic.expiration );
elseif remain_days < 0
    ltext{end+1} = sprintf( '%sEXPIRED: %s (%d days ago)', prefix, lic.expiration, -remain_days );
else
    ltext{end+1} = sprintf( '%sExpiration: %s (%g days remaining)', prefix, lic.expiration, remain_days );
end
if isequal( lic.status, 'VERIFIED'  ),
    ltext{end+1} = sprintf( '%sSignature: valid', prefix );
elseif isequal( lic.status, 'INVALID:SIGNATURE' ),
    ltext{end+1} = sprintf( '%sSignature: INVALID', prefix );
else
    ltext{end+1} = sprintf( '%sSignature: unverified', prefix );
end

function valid = quick_verify( lic, public_key )
valid = false;
try
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
    dsa.initVerify(public_key);
    dsa.update(unicode2native(message));
    if ~dsa.verify(lic.signature), return; end
    valid = true;
catch %#ok
end

function lic = full_verify( lic, public_key )
try 
    signature = lic.signature; 
    lic.signature = []; 
    lic.days_left = '-Inf';
catch %#ok
    signature = [];
end
if isempty( lic ),
    lic = struct( 'signature', [], 'status', 'EMPTY', 'days_left', '-Inf' );
    return
elseif ~isstruct( lic ) || numel( lic ) ~= 1,
    lic = struct( 'signature', [], 'status', 'CORRUPT', 'days_left', '-Inf' );
    return
elseif ~isfield( lic, 'name' ) || ~ischar( lic.name ) || isempty( lic.name ) || ...
       ~isfield( lic, 'organization' ) || ~ischar( lic.organization ) || isempty( lic.organization ) || ...
       ~isfield( lic, 'email' ) || ~ischar( lic.email ) || isempty( lic.email ) || ...
       ~isfield( lic, 'username' ) || ~ischar( lic.username ) || ...
       ~isfield( lic, 'hostid' ) || ~ischar( lic.hostid ) || ...
       ~isfield( lic, 'expiration' ) || ~ischar( lic.expiration ) || isempty( lic.expiration ) || ...
       ~isfield( lic, 'signature' ) || isempty( signature ) || ~isa( signature, 'int8' ),
    lic.status = 'CORRUPT';
    lic.days_left = '-Inf';
    return
end
if ~isempty( lic.username ),
    username = char(java.lang.System.getProperty('user.name'));
    if ~isequal( username, lic.username ),
        lic.status = 'INVALID:USER';
        return
    end
end
if ~isempty( lic.hostid ),
    lic.hostid_name = '';
    networks = java.net.NetworkInterface.getNetworkInterfaces();
    while networks.hasMoreElements(),
        ni = networks.nextElement();
        hostid = sprintf('%02x',rem(double(ni.getHardwareAddress())+256,256));
        if isequal( hostid, lic.hostid ),
            lic.hostid_name = char(ni.getName);
            break
        end
    end
    if isempty( lic.hostid_name ),
        lic.status = 'INVALID:HOSTID';
        return
    end
end
parser = java.text.SimpleDateFormat('yyyy-MM-dd');
try
    expire = parser.parse(lic.expiration);
catch %#ok
    lic.status = 'CORRUPT:EXPIRATION';
    return
end
today = java.util.Date;
lic.days_left = ceil( ( double(expire.getTime()) - double(today.getTime) ) / 86400000 );
if lic.days_left < 0,
    lic.status = 'EXPIRED';
    return
end
try
    message = sprintf( '%s|%s|%s|%s|%s|%s', lic.name, lic.organization, lic.email, lic.username, lic.hostid, lic.expiration );
catch %#ok
    lic.status = 'CORRUPT:DSIGNATURE';
    lic.days_left = -Inf;
    return
end
dsa = java.security.Signature.getInstance('SHA1withDSA');
dsa.initVerify(public_key);
dsa.update(unicode2native(message));
if ~dsa.verify(signature),
    lic.status = 'INVALID:SIGNATURE';
    lic.days_left = -Inf;
else
    lic.signature = signature;
    lic.status = 'VERIFIED';
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.





