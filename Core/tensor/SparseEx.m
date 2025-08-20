classdef(InferiorClasses=?sym) SparseEx<IAdditive&ICompare
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
        function obj=set_vk(key,val,size)
            obj = SparseEx();
            obj.key = key;
            obj.val = val;
            obj.size = size;
        end
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
                    obj.key=zeros(1,0);
                    obj.val=arg;
                else
                    obj.size = obj.size(1); % vector case
                    obj.key = idx(:);
                    obj.val = arg(idx(:));
                end
            else
                subs=cell(1,length(obj.size));
                [subs{:}] = ind2sub(obj.size, idx(:));
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
            ret = length(obj.size);
        end
        function ret=toMatrix(obj)
            if isempty(obj.size)
                ret = obj.val; % scalar case
                return;
            end
            ret = zeros([obj.size,1],'like',obj.val);
            subscript=mat2cell(obj.key,obj.Nelem,ones(1,obj.rank));
            Indices = sub2ind([obj.size,1], subscript{:});
            ret(Indices) = obj.val;
        end
        % 簡約化=同次項括り＋零係数項削除+式の簡約化
        function obj=C(obj,arg)
            arguments
                obj SparseEx
                arg.isZero function_handle = @(x)logical(x==0)
                arg.simplify_func = AlgebraConfig.H.simplify_func_Sparse;
                arg.level  {mustBeMember(arg.level, {'','low', 'medium','high'})} =''  % 最適化の度合い
            end
            obj=combineTerm(obj);
            if ~isempty(arg.level)&&strcmp(arg.level,'low')
                arg.simplify_func=@(x)x;
            end
            obj.val=arg.simplify_func(obj.val);
            obj=obj.removeZero(arg.isZero);
        end
        function obj=removeZero(obj,isZero)
            % Remove zero elements from the SparseEx object
            idx=find(isZero(obj.val));
            obj.val(idx)=[];
            obj.key(idx,:)=[];
            if isempty(obj.val)
                obj.val=0;
                obj.key=ones(1,obj.rank);
            end
        end
        function obj=combineTerm(obj)
            % Combine terms with the same key
            if isempty(obj.key)
                obj.val = sum(obj.val);
                return;
            end
            [subscripts, sortIdx] = sortrows(obj.key); 
            sortedValues = obj.val(sortIdx);
            cumulativeSum = cumsum(sortedValues); % cumulative sum in that order
            GrpEnds=[find(~all(subscripts(1:end-1,:)==subscripts(2:end,:),2))];  % group ends
            obj.val = cumulativeSum([GrpEnds;end]) - [0; cumulativeSum(GrpEnds)];                     % sums per group
            obj.key =subscripts([GrpEnds;end],:);
        end
        function ret= ge(arg1,arg2)
            % Greater than or equal comparison
            arguments
                arg1 SparseEx
                arg2 SparseEx
            end
            ret=arg1-arg2;
            ret.val=ret.val>=0;
        end
        function arg= not(arg)
            % Greater than or equal comparison
            arguments
                arg SparseEx
            end
            arg.val = ~arg.val; % Negate the values for inequality
        end
        function arg= abs(arg)
            % Greater than or equal comparison
            arguments
                arg SparseEx
            end
            arg.val = abs(arg.val); 
        end
        function ret= max(arg)
            % Greater than or equal comparison
            arguments
                arg SparseEx
            end
            ret = max(arg.val); 
        end

        function ret= eq(arg1,arg2)
            % Greater than or equal comparison
            arguments
                arg1 SparseEx
                arg2 SparseEx
            end
            ret=arg1-arg2;
            ret.val=ret.val==0;
        end
        function ret=eqD(arg1,arg2,tol)
            A=isapprox(arg1,arg2,tol);
            ret=all(A.val);
        end
        function ret= isapprox(arg1,arg2,tol)
            % Greater than or equal comparison
            arguments
                arg1 SparseEx
                arg2 SparseEx
                tol double
            end
            ret=arg1-arg2;
            ret.val=abs(ret.val)<=tol;
        end
        function ret = plus(arg1, arg2)
            % Addition of two SparseEx objects
            arguments
                arg1 SparseEx
                arg2 SparseEx
            end
            N = arg1.Nelem;
            ret = arg1;
            % [ret.key, ~, iC] = unique([arg1.key; arg2.key], 'rows');
            % if isempty(ret.key)
            %     ret.val=arg1.val+arg2.val;
            % else
            % ret.val(iC(N+1:end)) = arg2.val;
            % ret.val(iC(1:N)) = ret.val(iC(1:N)) + arg1.val;
            % end
            ret.val=[arg1.val;arg2.val];
            ret.key=[arg1.key;arg2.key];
            ret=ret.C;
        end

        function arg=uminus(arg)
            % Negation of SparseEx object
            arg.val = -arg.val;
        end
        function disp0(obj)
            builtin('disp', obj)
        end
        function disp(obj)
            % disp method: Display non-zero entries in sparse-like format
            n = obj.Nelem;
            if isempty(obj.size)
                fprintf('SparseEx scalar: %s\n', string(obj.val));
            else
                fprintf('SparseEx [%s] with %d non-zero entries:\n', ...
                    join(string(obj.size), ','), n);
                for i = 1:n
                    k = obj.key(i, :);
                    kstr = join(string(k), ',');
                    v = string(obj.val(i));
                    fprintf('  (%s)  %s\n', kstr, v);
                end
            end
        end
    end
end
