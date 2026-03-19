
% Symbolic tests for PolyCalc (times, power, subs)
syms s x
% 例: p(x) = 1 + sx + 3x^2, q(x) = x + x^2
a = sym([1 s 3]);      % coeffs of p(x)
b = sym([0 1 1]);      % coeffs of q(x)
D = numel(a) - 1;

p = a(1) + a(2)*x + a(3)*x^2;
q = b(1) + b(2)*x + b(3)*x^2;
%% Test 1: times (symbolic)
fprintf('Test 1: times (symbolic) ... ');


expr = expand(p * q);           % 真の積
% 低次からの係数列（不足分は 0 を含む）
coeff_ref = coeffsFixed(expr,D,x).';

coeff_calc = PolyCalc.times(a, b);   % 我々の実装

ok1 = all(simplify(coeff_calc - coeff_ref) == 0);
assert(ok1, 'Test 1 failed: times does not match symbolic reference.');
fprintf('OK\n');

%% Test 2: power (symbolic)
fprintf('Test 2: power (symbolic) ... ');

n = 3;                               % 任意の非負整数
expr = expand(p^n);
coeff_ref = coeffsFixed(expr,D,x).';   % truncate

coeff_calc = PolyCalc.power(a, n);

ok2 = all(simplify(coeff_calc - coeff_ref) == 0);
assert(ok2, 'Test 2 failed: power does not match symbolic reference.');
fprintf('OK\n');

%% Test 3: subs (symbolic composition)
fprintf('Test 3: subs (symbolic composition) ... ');

expr = expand(subs(p, x, q));   % p(q(x))
coeff_ref = coeffsFixed(expr,D,x).';     % truncate

coeff_calc = PolyCalc.subs(a, b);

ok3 = all(simplify(coeff_calc - coeff_ref) == 0);
assert(ok3, 'Test 3 failed: subs does not match symbolic reference.');
fprintf('OK\n');

%% 追加テスト: ゼロ・単位多項式での安全性（軽く）
fprintf('Test 4: edge cases (zero and one) ... ');

a_zero = sym([0 0 0]);
a_one  = sym([1 0 0]);

% 0 * q = 0
c1 = PolyCalc.times(a_zero, b);
assert(all(simplify(c1 - a_zero) == 0), '0*q != 0');

% 1 * q = q (truncate 同じ長さ)
c2 = PolyCalc.times(a_one, b);
assert(all(simplify(c2 - b) == 0), '1*q != q');

% q^0 = 1
c3 = PolyCalc.power(b, 0);
assert(all(simplify(c3 - a_one) == 0), 'q^0 != 1');

fprintf('OK\n');

fprintf('===== All symbolic tests passed. =====\n');