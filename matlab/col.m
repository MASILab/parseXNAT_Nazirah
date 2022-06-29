function y = col(x)
%COL - columizes a matrix x as in y=x(:)
%
% Syntax:  y = col(x)
%
% Inputs:
%    x - matrix to be columised
%
% Outputs:
%    y - column vector of all elements in x
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Author:  harrigr
% Date:    12-Aug-2015 13:28:08
% Version: 1.0
% Changelog:
%
% 12-Aug-2015 13:28:08 - initial creation
%
%------------- BEGIN CODE --------------

y = x(:);

%------------- END OF CODE --------------