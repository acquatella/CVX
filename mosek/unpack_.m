if exist( '~/Downloads/mosektoolsosx64x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolsosx64x86.tar.bz2 -C /tmp "*/libmosek64.7.0.dylib" "*/libiomp5.dylib" "*/mosekopt.mexmaci64"' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012a/mosekopt.mexmaci64 ~/Projects/CVX/trunk/mosek/maci64/mosekopt.mexmaci64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2013a/mosekopt.mexmaci64 ~/Projects/CVX/trunk/mosek/maci64/mosekopt801.mexmaci64' );
    system( 'mv -f /tmp/mosek/7/tools/platform/osx64x86/bin/* ~/Projects/CVX/trunk/mosek/maci64/' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolslinux64x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolslinux64x86.tar.bz2 -C /tmp "*/libmosek64.so.7.0" "*/libiomp5.so" "*/mosekopt.mexa64"' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2009b/mosekopt.mexa64 ~/Projects/CVX/trunk/mosek/a64/mosekopt.mexa64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012a/mosekopt.mexa64 ~/Projects/CVX/trunk/mosek/a64/mosekopt714.mexa64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012b/mosekopt.mexa64 ~/Projects/CVX/trunk/mosek/a64/mosekopt800.mexa64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2013a/mosekopt.mexa64 ~/Projects/CVX/trunk/mosek/a64/mosekopt801.mexa64' );
    system( 'mv -f /tmp/mosek/7/tools/platform/linux64x86/bin/* ~/Projects/CVX/trunk/mosek/a64/' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolslinux32x86.tar.bz2', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'tar jvfx ~/Downloads/mosektoolslinux32x86.tar.bz2 -C /tmp "*/libmosek.so.7.0" "*/libiomp5.so" "*/mosekopt.mexglx"' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2009b/mosekopt.mexglx ~/Projects/CVX/trunk/mosek/glx/mosekopt.mexglx' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012a/mosekopt.mexglx ~/Projects/CVX/trunk/mosek/glx/mosekopt714.mexglx' );
    system( 'mv -f /tmp/mosek/7/tools/platform/linux32x86/bin/* ~/Projects/CVX/trunk/mosek/glx/' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolswin32x86.zip', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'unzip -d /tmp ~/Downloads/mosektoolswin32x86.zip "*/libiomp5md.dll" "*/mosek7_0.dll" "*/mosekopt.mexw32"' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2009b/mosekopt.mexw32 ~/Projects/CVX/trunk/mosek/w32/mosekopt.mexw32' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012a/mosekopt.mexw32 ~/Projects/CVX/trunk/mosek/w32/mosekopt714.mexw32' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2013a/mosekopt.mexw32 ~/Projects/CVX/trunk/mosek/w32/mosekopt801.mexw32' );
    system( 'mv -f /tmp/mosek/7/tools/platform/win32x86/bin/* ~/Projects/CVX/trunk/mosek/w32/' );
    system( 'rm -rf /tmp/mosek' );
end
if exist( '~/Downloads/mosektoolswin64x86.zip', 'file' ),
    system( 'rm -rf /tmp/mosek' );
    system( 'unzip -d /tmp ~/Downloads/mosektoolswin64x86.zip "*/libiomp5md.dll" "*/mosek64_7_0.dll" "*/mosekopt.mexw64"' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2009b/mosekopt.mexw64 ~/Projects/CVX/trunk/mosek/w64/mosekopt.mexw64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2012a/mosekopt.mexw64 ~/Projects/CVX/trunk/mosek/w64/mosekopt714.mexw64' );
    system( 'mv -f /tmp/mosek/7/toolbox/r2013a/mosekopt.mexw64 ~/Projects/CVX/trunk/mosek/w64/mosekopt801.mexw64' );
    system( 'mv -f /tmp/mosek/7/tools/platform/win64x86/bin/* ~/Projects/CVX/trunk/mosek/w64/' );
    system( 'rm -rf /tmp/mosek' );
end
mpath = mfilename( 'fullpath' );
temp = strfind( mpath, filesep );
mpath = mpath( 1 : temp(end) - 1 );
odir = pwd;
switch computer,
    case 'MACI64',
        cd( [ mpath, filesep, 'maci64' ] );
        libs = dir('*.dylib');
        files = dir('*.mexmaci64');
        files = { libs.name, files.name };
        libs = { libs.name };
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
        files = dir('*.mexa64');
        files = { files.name };
        for k = files(:)',
            system( sprintf( 'chrpath -r ''$ORIGIN'' %s', k{1} ) );
            system( sprintf( 'chrpath -l %s', k{1} ) );
        end
        cd( [ mpath, filesep, 'glx' ] );
        files = dir('*.mexglx');
        files = { files.name };
        for k = files(:)',
            system( sprintf( 'chrpath32 -r ''$ORIGIN'' %s', k{1} ) );
            system( sprintf( 'chrpath32 -l %s', k{1} ) );
        end
end
cd( odir )
