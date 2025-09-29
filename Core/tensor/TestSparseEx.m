%{
   Concept of SparseEx

  SparseEx is a lightweight sparse array class in MATLAB, designed to represent and manipulate arrays (especially 2D matrices, extendable to n-dimensional tensors) using an explicit coordinate (COO) format.

   Core Idea

 - Store only the non-zero entries of an array as:
   - key: integer indices of each entry (rows of multi-indices)
   - val: the corresponding values
 - The full array is thus the sum of entries at those coordinates.

   Motivation
  MATLAB’s built-in sparse works only for numeric 2D matrices. SparseEx generalizes the concept to:
  - Arbitrary element types (val could be numbers, symbolic, custom objects)
  - Arbitrary rank tensors (not limited to 2D)
  - More transparent interface for algebraic experiments

   Design Principles
  - Minimal core: key, val, size
  - Operations defined algebraically:
    - Addition merges entries with same key
    - Multiplication supports scalars and 2D matrix product
    - Simplification removes zero entries
  - Extendibility: users can define their own element type in val if it supports +, *, etc.
  - Transparency: provide disp for inspection, and a future full() for conversion back to dense MATLAB arrays
  - Immutability by default: operations return new objects, not modifying operands

   Position in Workflow
  SparseEx serves as a prototype algebraic tensor/matrix backend:
  - For small dense inputs → compresses them to COO
  - For symbolic / custom objects → allows sparse manipulation where MATLAB sparse fails
  - For research code → provides clarity, debuggability, and customization

   1) Construction & Basic Invariants
  - Numeric input (2D): correct size, key, val, Nelem (zeros excluded)
  - Empty constructor: fields defaulting correctly (size unset or [], Nelem==0)
  - All‑zero numeric input: Nelem==0 after construction
  - Single nonzero at each corner / along diagonal: indices correct
  - Repeated values with same linear index (if added later): confirm de-dup on simplify
  - Large but sparse input: memory footprint vs dense
  - Non-2D numeric input: current behavior (should error or be TODO) is consistent

   2) simplify
  - Removes zeros in val and synchronizes key rows
  - Idempotence: calling twice doesn’t change result
  - Stability of order (if you care): index order preserved or defined

   3) Display (disp)
  - Prints expected header and entries for small examples
  - Handles Nelem==0 gracefully

   4) Nelem (Dependent prop)
  - Matches numel(val) after add/remove/simplify operations
  - Stays correct after plus, scalar mtimes

   5) plus
  - Happy path: add disjoint supports → union of key, untouched values
  - Overlapping supports → values summed at equal keys
  - Commutativity: A+B == B+A
  - Associativity: (A+B)+C == A+(B+C)
  - Identity: A + 0 == A (zero array of same size)
  - Size mismatch: throws with clear message
  - Simplification: zero sums vanish after simplify

   6) mtimes (current spec)
  - Scalar × SparseEx (left): scales val only
  - SparseEx × scalar (right): scales val only
  - 2D matmul: A*B equals dense baseline for random sparse A,B
  - Shape rules: A.size(2) == B.size(1); otherwise assert throws
  - Mixed zero/sparse: 0*A == 0, A*0 == 0
  - Distributivity (where defined): A*(B+C) == A*B + A*C
  - Associativity with scalar: (α*A)*β == α*(A*β)
  - NOTE: mtimes uses full(arg): add/verify existence of a full helper (see §10)

   7) Error Handling & Messages
  - Size mismatch on plus
  - Non-2D mtimes blocked with clear error
  - Bad inputs (non-numeric scalar, mismatched size types) give meaningful ids/messages

   8) Numerical Correctness vs MATLAB Baselines
  - Small random 2D matrices: compare SparseEx(A) * SparseEx(B) to A*B
  - Additive and multiplicative identities: compare to zeros/eye when representable
  - Symmetry cases (A symmetric): A*A' vs baseline

   9) Ordering & Uniqueness of Keys
  - Different row orders of identical (key,val) represent the same tensor after simplify
  - Duplicate keys upon construction (if supported later) collapse to single with summed val

   10) Missing/Auxiliary Methods (define & then test)
  - Your mtimes calls full(arg) but the class has no full method yet. Add:
    - full(obj): returns numeric dense array for numeric val (2D), else error
  - Tests:
    - full(SparseEx(A)) equals A for numeric A
    - Round-trip: SparseEx(full(obj)) equals obj (after simplify), for numeric cases

   11) Performance & Scalability (smoke/benchmark)
  - Construction time vs sparse for large random 2D inputs
  - plus and scalar mtimes time grow ~linearly with Nelem
  - A*B (2D) cost roughly proportional to baseline dense multiply for comparable density

   12) Edge & Corner Cases
  - 1×1 arrays, 1×N and N×1
  - Extremely rectangular shapes
  - Very large indices near dimension bounds
  - val containing NaN/Inf (numeric): preserved through plus, scalar mtimes
  - After simplify, class remains consistent when Nelem==0 (e.g., key becomes 0xD)

   13) Property-based / Randomized Checks
  - Random sparse A,B (varying density): algebraic laws (commutativity of plus, associativity where valid, distributivity)
  - Random sparse A and random scalars α,β: scalar laws

   14) API Consistency & Immutability Expectations
  - Methods don’t mutate inputs unintentionally (unless documented); results are new objects
  - After operations, size preserved/updated correctly

   15) Future-proofing Hooks (optional, if you’ll add later)
  - nD support: constructor, full, plus defined for D>2
  - Type-general val (non-numeric): ensure plus/scalar mtimes behavior is guarded/defined
  - Conversion: to/from MATLAB sparse for 2D numeric

   Test Implementation for SparseEx

  This script tests the SparseEx class according to the specification above
   scalar
a=SparseEx(3);
b=a+a;
disp(b)
b=a-a;
disp(b)

   Test 1: Construction & Basic Invariants
fprintf('Testing Construction & Basic Invariants...\n');


  Test numeric 2D input
A = [1 0 3; 0 2 0; 0 0 4];
sparse_A = SparseEx(A);
assert(isequal(sparse_A.size, [3, 3]), 'Size should be [3, 3]');
assert(sparse_A.Nelem == 4, 'Should have 4 non-zero elements');

  Test empty constructor
empty_sparse = SparseEx();
assert(empty_sparse.Nelem == 0, 'Empty constructor should have Nelem==0');

  Test all-zero input
zeros_A = zeros(3, 3);
sparse_zeros = SparseEx(zeros_A);
assert(sparse_zeros.Nelem == 0, 'All-zero matrix should have Nelem==0');

fprintf('✓ Construction tests passed\n\n');

   Test 2: simplify method
fprintf('Testing simplify method...\n');

  Create object with zeros that need removal
sparse_with_zeros = sparse_A;
sparse_with_zeros.val(end+1) = 0;   Add a zero
sparse_with_zeros.key(end+1, :) = [1, 1];   At position (1,1)

simplified = sparse_with_zeros.simplify();
assert(simplified.Nelem == sparse_A.Nelem, 'Simplify should remove zeros');

  Test idempotence
double_simplified = simplified.simplify();
assert(isequal(simplified.key, double_simplified.key), 'Simplify should be idempotent');
assert(isequal(simplified.val, double_simplified.val), 'Simplify should be idempotent');

fprintf('✓ Simplify tests passed\n\n');

   Test 3: Display (disp)
fprintf('Testing display...\n');
disp('Display of sparse_A:');
disp(sparse_A);
disp('Display of empty sparse:');
disp(empty_sparse);
fprintf('✓ Display tests completed\n\n');

   Test 4: Nelem property
fprintf('Testing Nelem property...\n');
assert(sparse_A.Nelem == length(sparse_A.val), 'Nelem should match length of val');
fprintf('✓ Nelem tests passed\n\n');

   Test 5: plus operation
fprintf('Testing plus operation...\n');

B = [0 1 0; 1 0 0; 0 0 1];
sparse_B = SparseEx(B);
%}
%{
>> s = struct(a=[-3,1;2,4 ],b=[1 2; 3 4] );
>> s.a*s.b
ans =
     0    -2
    14    20
>> A=calcTensorExpression('s.a{1,2}s.b{2,3}',[1,3])
A.disp0
A = 
SparseEx [2,2] with 3 non-zero entries:
  (1,2)  -2
  (2,1)  14
  (2,2)  20
  SparseEx のプロパティ:

     zero: 0
      key: [3×2 double]
      val: [3×1 double]
     size: [2 2]
    Nelem: 3
     rank: 2
>> A.disp
SparseEx [2,2] with 3 non-zero entries:
  (1,2)  -2
  (2,1)  14
  (2,2)  20
%}



%% test (2,2)
s = struct(a=[-3,1;2,4 ],b=[1 2; 3 4] );
expected=s.a*s.b;
A=calcTensorExpression('s.a{1,2}s.b{2,3}',[1,3]);
assert(isequal(A.toMatrix(), expected), 'Result does not match expected');
assert(isequal(A.size, [2, 2])&&isequal(A.Nelem, 3)&&isequal(A.rank, 2)&& ...
    isequal(A.val, [-2; 14; 20])&&isequal(A.key, [1, 2; 2, 1; 2, 2]));


%%


% Test addition
C = sparse_A + sparse_B;
expected_C = A + B;
assert(isequal(full(C), expected_C), 'Addition should match dense baseline');

% Test commutativity
C2 = sparse_B + sparse_A;
assert(isequal(full(C), full(C2)), 'Addition should be commutative');

% Test identity (adding zero)
zero_sparse = SparseEx(zeros(3, 3));
A_plus_zero = sparse_A + zero_sparse;
assert(isequal(full(A_plus_zero), A), 'A + 0 should equal A');

fprintf('✓ Plus operation tests passed\n\n');

%% Test 6: mtimes operation
fprintf('Testing mtimes operation...\n');

% Test scalar multiplication (left)
alpha = 2.5;
scaled_A = alpha * sparse_A;
expected_scaled = alpha * A;
assert(isequal(full(scaled_A), expected_scaled), 'Scalar multiplication should work');

% Test scalar multiplication (right)
scaled_A2 = sparse_A * alpha;
assert(isequal(full(scaled_A), full(scaled_A2)), 'Scalar multiplication should be commutative');

% Test matrix multiplication
D = sparse_A * sparse_B;
expected_D = A * B;
assert(isequal(full(D), expected_D), 'Matrix multiplication should match dense baseline');

fprintf('✓ Mtimes operation tests passed\n\n');

%% Test 7: Error Handling
fprintf('Testing error handling...\n');

% Test size mismatch in addition
try
    wrong_size = SparseEx([1 2; 3 4]); % 2x2
    result = sparse_A + wrong_size; % Should error (3x3 + 2x2)
    error('Should have thrown size mismatch error');
catch ME
    assert(contains(ME.message, 'size') || contains(ME.message, 'dimension'), ...
        'Should throw size mismatch error');
end

fprintf('✓ Error handling tests passed\n\n');

%% Test 8: Numerical Correctness
fprintf('Testing numerical correctness...\n');

% Random test matrices
rng(42); % For reproducibility
A_rand = rand(5, 4);
B_rand = rand(4, 3);
A_rand(A_rand < 0.7) = 0; % Make sparse
B_rand(B_rand < 0.7) = 0;

sparse_A_rand = SparseEx(A_rand);
sparse_B_rand = SparseEx(B_rand);

% Test multiplication
C_sparse = sparse_A_rand * sparse_B_rand;
C_dense = A_rand * B_rand;
assert(max(abs(full(C_sparse) - C_dense), [], 'all') < 1e-12, ...
    'Random matrix multiplication should be accurate');

fprintf('✓ Numerical correctness tests passed\n\n');

%% Test 9: full method
fprintf('Testing full method...\n');

% Test conversion back to dense
A_recovered = full(sparse_A);
assert(isequal(A_recovered, A), 'full() should recover original matrix');

% Test round-trip
sparse_recovered = SparseEx(A_recovered);
assert(isequal(full(sparse_recovered), A), 'Round-trip should preserve matrix');

fprintf('✓ Full method tests passed\n\n');

%% Test 10: Edge Cases
fprintf('Testing edge cases...\n');

% Test 1x1 matrix
single_val = SparseEx([5]);
assert(full(single_val) == 5, '1x1 matrix should work');

% Test 1xN and Nx1 vectors
row_vec = SparseEx([1 0 3 0 5]);
col_vec = SparseEx([1; 0; 3; 0; 5]);
assert(isequal(size(full(row_vec)), [1, 5]), 'Row vector size should be correct');
assert(isequal(size(full(col_vec)), [5, 1]), 'Column vector size should be correct');

fprintf('✓ Edge case tests passed\n\n');

%% Test 11: Algebraic Properties
fprintf('Testing algebraic properties...\n');

% Test associativity of addition
A1 = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));
A2 = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));
A3 = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));

left_assoc = (A1 + A2) + A3;
right_assoc = A1 + (A2 + A3);
assert(max(abs(full(left_assoc) - full(right_assoc)), [], 'all') < 1e-12, ...
    'Addition should be associative');

% Test distributivity
A_test = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));
B_test = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));
C_test = SparseEx(rand(3, 3) .* (rand(3, 3) > 0.5));

left_dist = A_test * (B_test + C_test);
right_dist = A_test * B_test + A_test * C_test;
assert(max(abs(full(left_dist) - full(right_dist)), [], 'all') < 1e-12, ...
    'Multiplication should distribute over addition');

fprintf('✓ Algebraic property tests passed\n\n');

%% Test 12: Immutability
fprintf('Testing immutability...\n');

original_A = sparse_A;
original_val = sparse_A.val;
original_key = sparse_A.key;

% Perform operations
result = sparse_A + sparse_B;
scaled = 3 * sparse_A;

% Check original unchanged
assert(isequal(sparse_A.val, original_val), 'Original val should be unchanged');
assert(isequal(sparse_A.key, original_key), 'Original key should be unchanged');

fprintf('✓ Immutability tests passed\n\n');

fprintf('All SparseEx tests completed successfully! ✓\n');
