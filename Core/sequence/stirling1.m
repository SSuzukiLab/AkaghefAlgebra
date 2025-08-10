function S = stirling1(m, n, c)
%STIRLING2 returns Stirling number of the 1st kind
%
%BACKGROUND Stirling numbers of the 1st kind show how many ways in which a
%set of m elements can be clustered to form n groups. This algorithm uses a
%recursive formula that is exact up to intmax('uint64'). Any Stirling
%number under this will be represented exactly and any number larger than
%this but less than realmax will be approximated using double precision
%(16-digits). This will solve Stirling numbers exactly for S(26,n), for any
%n and S(>26,n) for select n. Solutions are available for S(<220,n) for all
%n, and S(>=220,n) for select n. An error will indicate overflow. If the S
%is uint64 it is exact. If S is double, it may be an approximation good to
%16 digits.
%
%USAGE
%   S = stirling1( m ) 
%   returns S, a vector containing the number of ways in which a set of m
%   elements can be clustered to form 0,1,2,...,m groups
%
%   S = stirling1( m, n )
%   returns scalar S, the number of ways in which a set of m
%   elements can be clustered to form m groups
%
%   S = stirling1( m, n, 'double'). 
%   This convention is used internally, if it is empirically determined
%   that S > intmax('unit64'). You can use it directly to force double
%   precision.


arguments
    m 
    n = 0:m
    c {mustBeMember(c,{'uint64','double','sym'})}= 'double'
end
mMax=max(m);
S = diag(ones(mMax+1,1,c) );       % table of solutions S(m,n)
     
for r = 2:mMax+1
    for k = 2:r-1
        S(r,k) = S(r-1, k-1) + (r-2)*S(r-1,k);% 1-indexed
    end
end

S = S(m+1,n+1);

if isa(S,'double') && any(S > realmax('double'),"all")
    error('overflow');
elseif isa(S,'uint64') && any(S >= intmax('uint64'),"all")
    S = stirling2(m,n,'double');
end

