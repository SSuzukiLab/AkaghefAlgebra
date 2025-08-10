classdef SparseEx
    % SparseEx: A lightweight sparse array class using key-value representation
    %   Stores non-zero elements and their indices explicitly.

    properties
        key double       % NxD array: indices of non-zero elements
        val (:,1)        % Nx1 array: corresponding non-zero values
        size 
    end

    properties(Dependent)
        Nelem            % Number of non-zero elements
    end

    methods
        function obj = SparseEx(arg)
            % Constructor: Converts numeric array into SparseEx format
            sz = size(arg);
            if isnumeric(arg)
                obj.val = arg(:);
                tmp = arrayfun(@(x){1:x}, sz);
                obj.key = fliplr(table2array(combinations(tmp{:}))); % Generate linear indices
            end
        end

        function ret = get.Nelem(obj)
            % Return the number of stored non-zero elements
            ret = numel(obj.val);
        end
        function obj=simplify(obj)
            idx=find(obj.val==0);
            obj.val(idx)=[];
            obj.key(idx,:)=[];
        end

        function outputArg = method1(obj, inputArg)
            % Placeholder method (not used)
            outputArg = obj.Property1 + inputArg;
        end

        function ret = plus(arg1, arg2)
            % Addition of two SparseEx objects
            arguments
                arg1 SparseEx
                arg2 SparseEx
            end
            N = arg1.Nelem;
            ret = arg1;
            [ret.key, ~, iC] = unique([arg1.key; arg2.key], 'rows');
            ret.val(iC(N+1:end)) = arg2.val;
            ret.val(iC(1:N)) = ret.val(iC(1:N)) + arg1.val;
            ret=ret.simplify;
        end

        function disp(obj)
            % disp method: Display non-zero entries in sparse-like format
            n = obj.Nelem;
            fprintf('SparseEx with %d non-zero entries:\n', n);
            for i = 1:n
                k = obj.key(i, :);
                kstr = join(string(k), ',');
                v = obj.val(i);
                fprintf('  (%s)  %g\n', kstr, v);
            end
        end
    end
end
