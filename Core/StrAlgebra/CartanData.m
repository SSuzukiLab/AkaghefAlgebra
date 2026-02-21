classdef CartanData
    % CartanData - A class to represent simple Lie Algebras (Classical Types)
    %   Properties:
    %       TypeCharacter : Char ('A', 'B', 'C', 'D') - Classification symbol
    %       N_DynkinNodes : Integer (n, the number of simple roots/nodes)
    %       Type          : String (e.g., 'A3') - Dependent property
    %       SimpleRoots   : Matrix (N_DynkinNodes x Dimension)
    %       CartanMatrix  : Matrix (N_DynkinNodes x N_DynkinNodes)
    %       Symmetrizer   : Vector (1 x N_DynkinNodes)
    
    properties (SetAccess = private)
        TypeCharacter % Stores the classification char (A, B, C, D)
        N_DynkinNodes   % Stores the number of simple roots (n)
        SimpleRoots   % Alpha
        CartanMatrix  % A
        Symmetrizer   % d
    end
    
    properties (Dependent)
        Type % Dynamically generated string (e.g., "A3")
    end
    
    methods
        function obj = CartanData(typeString)
            % Constructor: CartanData('A3')
            
            % Parse the input string
            tokens = regexp(typeString, '([A-Z])(\d+)', 'tokens', 'once');
            
            if isempty(tokens)
                error('Invalid type string format. Expected format: A1, B3, C4, etc.');
            end
            
            typeChar = tokens{1};
            rankInt = str2double(tokens{2});
            
            if isnan(rankInt) || rankInt < 1
                 error('Invalid size/rank detected. It must be a positive integer.');
            end
            
            % Assign to stored properties
            obj.TypeCharacter = upper(typeChar);
            obj.N_DynkinNodes = rankInt; % Updated property name
            
            % 1. Construct Simple Roots using TypeCharacter
            obj.SimpleRoots = obj.constructRoots(obj.TypeCharacter, obj.N_DynkinNodes);
            
            % 2. Compute Cartan Matrix
            obj.CartanMatrix = obj.computeCartanFromRoots(obj.SimpleRoots);
            
            % 3. Compute Symmetrizer
            obj.Symmetrizer = obj.computeSymmetrizer(obj.SimpleRoots);
        end
        
        function val = get.Type(obj)
            % Getter for dependent property Type
            % Combines TypeCharacter and N_DynkinNodes into a string
            val = sprintf('%s%d', obj.TypeCharacter, obj.N_DynkinNodes);
        end
        
        function res = verifySymmetrization(obj)
            % Verifies if diag(d) * A is a symmetric matrix
            D = diag(obj.Symmetrizer);
            A = obj.CartanMatrix;
            S = D * A;
            
            res = isequal(S, S');
            if ~res
                disp('Symmetry check failed.');
                disp(S);
            end
        end
        
        function disp(obj)
            % Custom display method uses the dynamic Type property
            fprintf('Lie Algebra Type: %s\n', obj.Type);
            fprintf('   (Character: %s, N_DynkinNodes: %d)\n', obj.TypeCharacter, obj.N_DynkinNodes);
            fprintf('------------------------\n');
            disp('Simple Roots (Alpha):');
            disp(obj.SimpleRoots);
            disp('Cartan Matrix (A):');
            disp(obj.CartanMatrix);
            disp('Symmetrizer (d):');
            disp(obj.Symmetrizer);
        end
    end
    
    methods (Access = private)
        function roots = constructRoots(~, type, n)
            % Generates simple roots in standard ambient space
            
            switch type
                case 'A'
                    dim = n + 1;
                    I = eye(dim);
                    roots = zeros(n, dim);
                    for i = 1:n
                        roots(i, :) = I(i, :) - I(i+1, :);
                    end
                    
                case 'B'
                    dim = n;
                    I = eye(dim);
                    roots = zeros(n, dim);
                    for i = 1:n-1
                        roots(i, :) = I(i, :) - I(i+1, :);
                    end
                    roots(n, :) = I(n, :); 
                    
                case 'C'
                    dim = n;
                    I = eye(dim);
                    roots = zeros(n, dim);
                    for i = 1:n-1
                        roots(i, :) = I(i, :) - I(i+1, :);
                    end
                    roots(n, :) = 2 * I(n, :); 
                    
                case 'D'
                    dim = n;
                    I = eye(dim);
                    roots = zeros(n, dim);
                    for i = 1:n-1
                        roots(i, :) = I(i, :) - I(i+1, :);
                    end
                    roots(n, :) = I(n-1, :) + I(n, :);
                    
                otherwise
                    error('Only classical types A, B, C, D are implemented.');
            end
        end
        
        function A = computeCartanFromRoots(~, roots)
            [n, ~] = size(roots);
            A = zeros(n, n);
            
            for i = 1:n
                for j = 1:n
                    alpha_i = roots(i, :);
                    alpha_j = roots(j, :);
                    
                    dot_prod = dot(alpha_i, alpha_j);
                    norm_sq  = dot(alpha_j, alpha_j);
                    
                    val = 2 * dot_prod / norm_sq;
                    A(i, j) = round(val);
                end
            end
        end
        
        function d = computeSymmetrizer(~, roots)
            [n, ~] = size(roots);
            sq_lengths = zeros(1, n);
            
            for i = 1:n
                sq_lengths(i) = dot(roots(i, :), roots(i, :));
            end
            
            min_len = min(sq_lengths);
            d = sq_lengths / min_len;
        end
    end
end