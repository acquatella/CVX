if exist( '~/Downloads/mosektoolsosx64x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolsosx64x86.tar.bz2 -C /tmp "*/libmosek64.7.0.dylib" "*/libiomp5.dylib" "*/mosekopt.mexmaci64"' );
    system( 'find /tmp/mosek -type f -exec mv -f {} ~/cvx/trunk/mosek/maci64 \;' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolslinux64x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolslinux64x86.tar.bz2 -C /tmp "*/libmosek64.so.7.0" "*/libiomp5.so" "*/mosekopt.mexa64"' );
    system( 'find /tmp/mosek -type f -exec mv -f {} ~/cvx/trunk/mosek/a64 \;' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolslinux32x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolslinux32x86.tar.bz2 -C /tmp "*/libmosek.so.7.0" "*/libiomp5.so" "*/mosekopt.mexglx"' );
    system( 'find /tmp/mosek -type f -exec mv -f {} ~/cvx/trunk/mosek/glx \;' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolswin32x86.zip', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'unzip -d /tmp ~/Downloads/mosektoolswin32x86.zip "*/libiomp5md.dll" "*/mosek7_0.dll" "*/mosekopt.mexw32"' );
    system( 'find /tmp/mosek -type f -exec mv -f {} ~/cvx/trunk/mosek/w32 \;' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolswin64x86.zip', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'unzip -d /tmp ~/Downloads/mosektoolswin64x86.zip "*/libiomp5md.dll" "*/mosek64_7_0.dll" "*/mosekopt.mexw64"' );
    system( 'find /tmp/mosek -type f -exec mv -f {} ~/cvx/trunk/mosek/w64 \;' );
    system( 'rm -rf /tmp/mosek' );
end
mpath = mfilename( 'fullpath' );
temp = strfind( mpath, filesep );
mpath = mpath( 1 : temp(end) - 1 );
odir = pwd;
switch computer,
    case 'MACI64',
        libs = { 'libmosek64.7.0.dylib', 'libiomp5.dylib' };
        cd( [ mpath, filesep, 'maci64' ] );
        files = [ dir('*.dylib' ); dir('*.mexmaci64') ];
        files = { files.name };
        for k = files(:)',
            k = k{1};
            for l = libs(:)',
                l = l{1};
                if strcmp( l, k ),
                    system( sprintf( 'install_name_tool -id %s %s', k, k ) );
                end
                system( sprintf( 'install_name_tool -change "%s" "@loader_path/%s" %s', l, l, k ) );
                system( sprintf( 'install_name_tool -change "@rpath/%s" "@loader_path/%s" %s', l, l, k ) );
            end
        end
        system( [ 'otool -L', sprintf( ' %s', files{:} ) ] );
    case 'GLNXA64',
        cd( [ mpath, filesep, 'a64' ]  );
        system( 'chrpath -r "." mosekopt.mexa64' );
        system( 'chrpath -l mosekopt.mexa64' );
        cd( [ mpath, filesep, 'glx' ] );
        system( 'chrpath32 -r "." mosekopt.mexglx' );
        system( 'chrpath32 -l mosekopt.mexglx' );
end
cd( odir )
