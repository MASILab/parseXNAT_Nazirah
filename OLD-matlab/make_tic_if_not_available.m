function make_tic_if_not_available
%MAKE_TIC_IF_NOT_AVAILABLE - Creates a tic if one isn't available
% Syntax:  make_tic_if_not_available
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: tprintf
%
% Author:  plassaaj
% Date:    09-Mar-2015 12:55:55
% Version: 1.0
% Changelog:
%
% 09-Mar-2015 12:55:55 - initial creation
%
%------------- BEGIN CODE --------------

try
    et = toc;
catch 
    tic;
    return;
end

%------------- END OF CODE --------------