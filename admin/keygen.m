lic = struct;
lic.name         = input( 'Name: ', 's' );
lic.organization = input( 'Organization: ', 's' );
lic.email        = input( 'Email address: ', 's' );
lic.license_type = input( 'License type: ', 's' );
lic.username = {};
while true,
    username = input( 'Username: ', 's' );
    if isempty( username ), break; end
    lic.username{end+1} = username;
end
lic.hostid = {};
while true,
    hostid = input( 'Host ID: ', 's' );
    if isempty( hostid ), break; end
    lic.hostid{end+1} = hostid;
end
if isempty( lic.username ),
    m_username = '';
elseif iscell( lic.username ),
    m_username = sprintf( '%s,', lic.username{:} );
    m_username(end) = [];
else
    m_username = lic.username;
end
if isempty( lic.hostid ),
    m_hostid = '';
elseif iscell( lic.hostid ),
    m_hostid = sprintf( '%s,', lic.hostid{:} );
    m_hostid(end) = [];
else
    m_hostid = lic.hostid;
end
lic.expiration = input( 'Expiration date (YYYY-MM-DD): ', 's' );
message = sprintf( '%s|', lic.name, lic.organization, lic.email, lic.license_type, m_username, m_hostid, lic.expiration );
keys = load('keys');
privkey = keys.keypair.getPrivate();
dsa = java.security.Signature.getInstance('SHA1withDSA');
dsa.initSign(privkey);
dsa.update(unicode2native(message,'UTF-8'));
lic.signature = dsa.sign()';
