function code = getfourcc(fourcc)

code = [char(bitshift(fourcc,0,8)) char(bitshift(fourcc,-8,8)) char(bitshift(fourcc,-16,8)) char(bitshift(fourcc,-24,8))];


%   $Revision: 1.1.6.2 $  $Date: 2004/03/30 13:06:56 $
%   Copyright 1984-2003 The MathWorks, Inc.

