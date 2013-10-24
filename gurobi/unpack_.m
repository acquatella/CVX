dd = dir('~/Downloads/gurobi*_mac64.pkg' );
if ~isempty(dd),
    system( 'rm -rf /tmp/gurobi' );
    system( 'mkdir /tmp/gurobi' );
    system( [ 'xar -C /tmp/gurobi -xf ~/Downloads/', dd(1).name ] );
    system( 'tar xfz /tmp/gurobi/*/Payload -C /tmp/gurobi' );
    system( 'tar xfz /tmp/gurobi/*.tar.gz -C /tmp/gurobi' );
    system( 'mv -f /tmp/gurobi/*/mac64/bin/grbgetkey ~/Projects/CVX/trunk/gurobi/maci64/' );
    system( 'svn rm ~/Projects/CVX/trunk/gurobi/maci64/libgurobi*.so' );
    system( 'mv -f /tmp/gurobi/*/mac64/lib/libgurobi*.so ~/Projects/CVX/trunk/gurobi/maci64/' );
    system( 'svn add ~/Projects/CVX/trunk/gurobi/maci64/libgurobi*.so' );
    system( 'mv -f /tmp/gurobi/*/mac64/matlab/gurobi.mexmaci64 ~/Projects/CVX/trunk/gurobi/maci64/' );
    system( 'rm -rf /tmp/mosek' );
end
dd = dir('~/Downloads/gurobi*_linux64.tar.gz');
if ~isempty(dd),
    system( 'rm -rf /tmp/gurobi' );
    system( 'mkdir /tmp/gurobi' );
    system( [ 'tar xfz ~/Downloads/', dd(1).name, ' -C /tmp/gurobi' ] );
    system( 'mv -f /tmp/gurobi/*/linux64/bin/grbgetkey ~/Projects/CVX/trunk/gurobi/a64/' );
    system( 'svn rm ~/Projects/CVX/trunk/gurobi/a64/libgurobi*.so' );
    system( 'cp -fH /tmp/gurobi/*/linux64/lib/libgurobi*.so ~/Projects/CVX/trunk/gurobi/a64/' );
    system( 'svn add ~/Projects/CVX/trunk/gurobi/a64/libgurobi*.so' );
    system( 'mv -f /tmp/gurobi/*/linux64/matlab/gurobi.mexa64 ~/Projects/CVX/trunk/gurobi/a64/' );
    system( 'rm -rf /tmp/mosek' );
end
mpath = mfilename( 'fullpath' );
temp = strfind( mpath, filesep );
mpath = mpath( 1 : temp(end) - 1 );
odir = pwd;
switch computer,
    case 'GLNXA64',
        cd( [ mpath, filesep, 'a64' ]  );
        files = dir('*.mexa64');
        files = { files.name };
        for k = files(:)',
            system( sprintf( 'chrpath -r "." %s', k{1} ) );
            system( sprintf( 'chrpath -l %s', k{1} ) );
        end
end
cd( odir )
