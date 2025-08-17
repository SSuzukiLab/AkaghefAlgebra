classdef(InferiorClasses=?sym) SparseEx<IAdditive
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
    methods(Static)
        function obj=convert(arg,Elem2NonzeroIndices)
            if isa(arg, 'SparseEx')
                obj = arg; % Already a SparseEx object
                return;
            end
            obj=SparseEx;
            if nargin<2
                Elem2NonzeroIndices = @find;
            elseif isa(Elem2NonzeroIndices,'double')
                Elem2NonzeroIndices = @(x) find(~eqD(x,0,Elem2NonzeroIndices));
            elseif isequal(Elem2NonzeroIndices,"full")
                Elem2NonzeroIndices =@(x) 1:numel(x);
            end
            % obj.zero = zeros(1,1,'like',arg);
            obj.size = size(arg);
            idx= Elem2NonzeroIndices(arg);
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
    methods
        function obj = SparseEx(arg)
            % Constructor: Converts numeric array into SparseEx format
            if nargin>0&&~isa(arg, 'SparseEx')
                % issue:sparseExへの変換メソッドを持っていた場合はどうする？
                obj=SparseEx.convert(arg);
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
