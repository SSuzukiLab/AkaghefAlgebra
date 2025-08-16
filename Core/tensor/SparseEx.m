classdef(InferiorClasses=?sym) SparseEx
    % SparseEx: A lightweight sparse array class using key-value representation
    %   Stores non-zero elements and their indices explicitly.

    properties
        zero (1,1) =0; % Default value for zero elements
        key double       % NxR array: indices of non-zero elements
        val (:,1)        % Nx1 array: corresponding non-zero values
        size 
    end

    properties(Dependent)
        Nelem            % Number of non-zero elements
        rank             % rank of tensor
    end

    methods
        function obj = SparseEx(arg)
            % Constructor: Converts numeric array into SparseEx format
            if ~isa(arg, 'SparseEx')
                obj.zero = zeros(1,1,'like',arg);
                obj.size = size(arg);
                if AlgebraConfig.H.SP_elim_zero 
                    idx= find(arg ~= obj.zero);
                else
                    idx = 1:numel(arg);
                end
                if isempty(arg)
                    obj.key = [];
                    obj.val = [];
                elseif length(obj.size) ==2&& obj.size(2) == 1
                    if obj.size(1) == 1
                        obj.size = []; % scalar case
                    else
                        obj.size = obj.size(1); % vector case
                    end
                    obj.key = idx(:); 
                    obj.val = arg(idx(:));
                else
                    subs=cell(1,length(obj.size));
                    [subs{:}] = ind2sub(obj.size, idx); 
                    obj.key = horzcat(subs{:});
                    obj.val = arg(idx); 
                end
            end
        end

        function ret = get.Nelem(obj)
            % Return the number of stored non-zero elements
            ret = numel(obj.val);
        end
        function ret = get.rank(obj)
            % Return the rank of the tensor (number of dimensions)
            ret = size(obj.key, 2);
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
                v = string(obj.val(i));
                fprintf('  (%s)  %s\n', kstr, v);
            end
        end
    end
end
