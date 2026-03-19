classdef PolyCalc
    % For usage, refer PolyCalcTest.m
    methods(Static)
        function c = times(a, b)
            % TIMES Multiply two polynomials given by coefficient vectors of same length.
            %   a, b : vectors of coefficients (lowest degree first)
            %          must have the same length.
            %   c    : result coefficients, truncated to the same length as a and b.

            % Ensure row vectors
            a = a(:).';
            b = b(:).';

            if numel(a) ~= numel(b)
                error('times:SizeMismatch', ...
                    'Input polynomials must have the same length.');
            end

            n = numel(a);
            D = n - 1;  % maximum degree inferred from input length

            % Type-aware zero for numeric/symbolic/mixed
            z0 = a(1)*0 + b(1)*0;
            c  = repmat(z0, 1, n);  % length D+1 = n

            % Manual truncated convolution (only terms up to degree D)
            % p(x) = sum_i a(i) x^(i-1), q(x) = sum_j b(j) x^(j-1)
            % c(k+1) = sum_{i+j-2 = k} a(i)b(j), 0 <= k <= D
            for i = 1:n
                ai = a(i);
                if isequal(ai, 0)
                    continue;   % skip zero coefficient (numeric case)
                end
                % j index such that degree k = i+j-2 <= D  =>  j <= D - i + 2
                jmax = min(n, D - i + 2);
                if jmax < 1
                    continue;
                end
                j = 1:jmax;
                k = i + j - 2;              % degrees for these pairs
                c(k+1) = c(k+1) + ai * b(j); % vectorized update
            end
        end

        function c = power(a, n)
            % POWER Raise polynomial to integer power with truncation.
            %   a : coeffs of p(x)
            %   n : nonnegative integer
            %   c : coeffs of p(x)^n, truncated to same length as a.

            % Ensure row
            a = a(:).';

            if ~isscalar(n) || n < 0 || n ~= floor(n)
                error('power:InvalidExponent', ...
                      'Exponent must be a nonnegative integer scalar.');
            end

            m = numel(a);
            D = m - 1; %#ok<NASGU> % kept for symmetry / future use

            % Type-aware zero and one
            z0  = a(1)*0;               % 0 in the coefficient field
            one = [z0+1, repmat(z0, 1, m-1)];  % polynomial "1"

            if n == 0
                c = one;
                return;
            end

            base = a;
            c    = one;
            k    = n;

            % Exponentiation by squaring, using truncated multiplication
            while k > 0
                if bitand(k, 1)
                    c = PolyCalc.times(c, base);
                end
                k = bitshift(k, -1);
                if k > 0
                    base = PolyCalc.times(base, base);
                end
            end
        end

        function r = subs(a, b)
            % SUBS Polynomial composition r(x) = p(q(x)) with fixed max degree.
            %   a : coeffs of p(x)
            %   b : coeffs of q(x)
            %   Both vectors must have the same length; degree inferred from length-1.
            %   r : coeffs of p(q(x)), truncated to the same length.

            % Ensure row vectors
            a = a(:).';
            b = b(:).';

            if numel(a) ~= numel(b)
                error('subs:SizeMismatch', ...
                    'Input polynomials must have the same length.');
            end

            m = numel(a);
            D = m - 1;

            % Type-aware zero
            z0 = a(1)*0 + b(1)*0;
            r  = repmat(z0, 1, m);

            % q^0 = 1
            pwr = [z0+1, repmat(z0, 1, m-1)];

            % p(q(x)) = sum_{i=0}^D a_i * q(x)^i
            for i = 0:D
                ai = a(i+1);
                if ~isequal(ai, 0)
                    r = r + ai * pwr;
                end
                if i < D
                    % q^{i+1} = q^i * q, truncated
                    pwr = PolyCalc.times(pwr, b);
                end
            end
        end
    end
end
