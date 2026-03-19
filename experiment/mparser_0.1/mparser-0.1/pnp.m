% function pnp( x )
%
% "Parse-and-print"
%
% A simple script that wraps the mparser MEX file.  It parses a
% file, checks for parse errors, and then prints out the parsed
% file to the screen.
%
%
% Copyright(c) 2011 David Wingate
% 
% This file is part of the mparser package.
% See the file COPYING for license details.
%

function pnp( x )
  [r,a] = mparser( x, 1 );
  if ( r == 0 )
    matlab_ast_print( a );
  else
    % parse error
    a
  end;
  
end
