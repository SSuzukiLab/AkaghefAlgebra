function f = polyfitExact(y, dmax, x)
% POLYFITEXACT  Exact-fitting polynomial (minimal degree) as symfun
%   f = exactPolyFit(y, dmax, x)
%     y    : vector of values y_i (length N)
%     dmax : maximum allowed degree (default N-1)
%     x    : abscissae (default 0:N-1)
%
%   Returns symfun f(n) of the smallest degree ≤ dmax that satisfies
%   f(x_i) == y_i for all i, using exact (rational) arithmetic.
%
%   Example:
%     y = [0;0;-1;-13;-58;-170;-395;-791];
%     f = exactPolyFit(y);          % minimal degree, default x=0:7
%     simplify(subs(f, 0:7).')      % verify exact fit

    arguments
        y sym {mustBeVector}
        dmax (1,1) double {mustBeNonnegative} = numel(y)-1
        x {mustBeVector} = 0:length(y)-1
    end

    y = y(:);
    N = numel(y);
    if numel(x) ~= N, error('x and y must have the same length.'); end

    dmax = min(dmax, N-1);

    % Use exact rationals
    % x = sym(x, 'r');
    % y = sym(y, 'r');

    for k = 0:dmax
        % Vandermonde (ascending powers 0..k)
        A = sym(zeros(N, k+1));
        for j = 0:k
            A(:, j+1) = x.^j;
        end

        % Symbolic solve (works even if overdetermined); check exactness
        % c = A \ y;
        
        % Temporarily suppress warnings for inconsistent systems
        warning('off', 'symbolic:mldivide:InconsistentSystem');
        warning('off','symbolic:sym:isAlways:TruthUnknown')
        c=mldivide(A,y);
        if all(isAlways(simplify(A*c - y) == 0))
            syms n
            poly = poly2sym(flip(c.'), n);  % flip to descending powers
            f = symfun(poly, n);
            warning('on', 'symbolic:mldivide:InconsistentSystem');
            warning('on','symbolic:sym:isAlways:TruthUnknown')
            if k==N-1
                warning('polynomial may be overfitting')
            end
            return
        end
    end

    error('No exact polynomial of degree ≤ %d fits all points. Try increasing dmax (≤ %d).', dmax, N-1);
end