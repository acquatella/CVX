function cvx_cd( subdir )

%CVX_CD   Change current working directory to a CVX subdirectory. 
%   CVX_CD executes CD <cvx_dir>, where <cvx_dir> is the root directory of the
%   cvx package.
%
%   CVX_CD <subdir> executes CD <cvx_dir>/<subdir>. Useful subdirectories of the
%   cvx package include:
%        doc/          documentation
%        examples/     sample cvx models
%   So for example, 
%       CVX_CD doc 
%   sets the current directory to the documentation directory of cvx.

oldd = pwd;
cd(cvx_where);
if nargin ~= 0,
	if exist(subdir,'dir'),
	    cd(subdir);
	else
	    cd(oldd);
	    error( sprintf( 'Cannot CD to %s%s (Name is nonexistent or not a directory).', cvx_where, subdir ) );
	end
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
