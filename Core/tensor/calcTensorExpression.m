function val=calcTensorExpression(str,ord,Elem2NonzeroIndices)
    arguments
        str (1,:) char
        ord (1,:) double=[]
        % Elem2NonzeroIndices: Function to find nonzero indices of elements
        % see SparseEx.convert for details
        Elem2NonzeroIndices =@find; 
    end

     % Enumerates all matching term index tuples [t1 ... tN] from tensors in table T,
    % but only for nonzero elements in T.value. Contracts axes based on matching
    % dummy variable labels in legList.
    %
    % T: table with variables:
    %   - labels (string)         : tensor expression
    %   - legList (cell array)    : axis dummy variables
    %   - value (cell array)      : SparseEx object
    %   - rank (scalar)           : tensor rank
    %   - indexTable (cell array)       : table of subscripts and coefficients of terms of tensor
    %   - ID (uint)               : numbering of input tensors
    %
    T=parseFormulaArgument(str);
    % Refactor to avoid cross-frame eval. Treat variables as function inputs/outputs. This is fastest, safest, and debuggable.
    % Use nested functions when you need true shared workspace (lexical scoping) instead of dynamic scoping.
    for i=1:height(T)
        try % Try evaluating the expression in the caller workspace
            T.value{i} = evalin('caller', T.labels(i));
        catch
            disp(T)
            error('Expression "%s" cannot be evaluated in the base workspace.', T.labels(i));
        end
    end
    T=evaluateParsedTable(T,Elem2NonzeroIndices);
    T=setupIndexTable(T);
    T=calcContractOrder(T);
    [tbl,T]=calcIndexTable(T,ord);
    val=multiplyCoeff(tbl,T,ord);
end
function T = parseFormulaArgument(arg)
    % Match chemical symbols with corresponding indices
    pattern = '(?<label>[A-Za-z\(\)\^\-\d]+)\{(?<indices>[0-9,]+)\}';
    pattern = '(?<label>[A-Za-z\(\)\^\-\d]+(?:\.[A-Za-z\(\)\^\-\d]+)?)\{(?<indices>[0-9,]+)\}';
    pattern = '(?<label>[A-Za-z\(\)\^\-\d\.]+)\{(?<indices>[0-9,]+)\}';
    tokens = regexp(arg, pattern, 'names');

    labels = string({tokens.label})';
    legList={tokens.indices};
    legList=cellfun(@(str){str2double(strsplit(str, ','))},legList).';
    ID=(1:length(labels)).';

    T = table(ID,labels, legList);
    % T.rank=cellfun(@length,T.legList);
end
function T=evaluateParsedTable(T,Elem2NonzeroIndices)
    % Validate the parsed table where each label is an expression
    % (not just a variable name), which must be evaluable in the caller workspace,
    % and its rank (i.e., ndims) must equal rank(i)
    idx=horzcat(T.legList{:});
    [vals, ~, grp] = unique(idx);           % Get unique values and group indices
    counts = accumarray(grp, 1)>2;          % Count occurrences
    if any(counts)
        error("contraction indices:%s must appear just 2 times",mat2str(vals(counts)))
    end
    values=repmat(SparseEx(),height(T),1);
    T.size=cell(height(T),1);
    for i = 1:height(T)
        expr = T.labels(i);
        value=T.value{i};
        value = SparseEx.convert(value, Elem2NonzeroIndices);
        values(i) = value;

        expectedRank = numel(value.size);
        actualRank = numel(T.legList{i});
        if expectedRank ~= actualRank
            warning('Expression "%s" has rank %d, expected %d.', expr, actualRank, expectedRank);
        end
    end
    T.value = values;
end
function T=setupIndexTable(T)
    T.indexTable=cell(height(T),1);
    % 各テンソル
    for i = 1:height(T)
        % find index of tensor to contract inside a tensor
        repeatedIdx =T.legList{i};
        [~, uniqueIdx] = ismember(repeatedIdx, repeatedIdx);
        uniqueIdx = unique(uniqueIdx);
        repeatedIdx(uniqueIdx) = [];
        
        % calc inside contraction
        subscripts = T.value(i).key;
        validIdx = true(T.value(i).Nelem,1);
        for idx = repeatedIdx
            pos = find(T.legList{i} == idx,2);
            validIdx = validIdx & (subscripts(:,pos(1)) == subscripts(:,pos(2)));
        end
        subscripts = subscripts(validIdx,uniqueIdx);
        nonzeroIdx = find(validIdx);
        % Create a table with indices and subscripts
        T.indexTable{i} = array2table([nonzeroIdx, subscripts], ...
            VariableNames=["i" + i, "v" + T.legList{i}(uniqueIdx)]);
    end
end
function [tbl,T]=calcIndexTable(T,ord)
    allVars=horzcat(ord,T.legList{:});
    % assert all element in allVars appear just 2 times
    uniqueVars = unique(allVars(:));
    mat=sum(uniqueVars==allVars ,2);
    if ~all(mat==2)
        error('All legs of tensor must appear exactly 2 times: %s', ...
        mat2str(uniqueVars(mat~=2)'));
    end
    tbl=T.indexTable{1};
    for i=2:height(T)
        try
            tbl=innerjoin(tbl,T.indexTable{i});
        catch %if no contraction btw T.indexTable{i} and tbl, calc kronecker product 
            indexTable=T.indexTable{i};
            C=combinations(1:height(tbl),1:height(indexTable));
            C=rowfun(@(x,y)[tbl(x,:),indexTable(y,:)],C);
            tbl=C{:,1};
        end
    end

end
function T=calcContractOrder(T)
    % Calculate the contraction order for the tensor network
    if isempty(T), return; end
    rem=1:height(T);
    acc=[];
    order=[];
    Ncontract=[];
    while ~isempty(rem)
        % Find the next tensor to contract
        [Ncontract(end+1), next] = max(cellfun(@(L)sum(ismember(acc,L)),T.legList(rem)));
        % Update the contraction order
        order(end+1) = rem(next);
        acc=[acc,T.legList{rem(next)}];
        % Remove the contracted tensor from the remaining list
        rem(next) = [];
    end
    T=T(order,:);
    if any(Ncontract(2:end)==0)
        % disp("no contraction product included")
        % issue: evaluate using kronecker product after multiplying
    end
end
function val=multiplyCoeff(tbl,T,ord)
    % Multiply coefficients in tbl with corresponding values in T.value
    % and return the result as a string expression.
    uniqueVars = unique(cat(2, T.legList{:}));
    % tbl=removevars(tbl, "v"+setdiff(uniqueVars,ord));
    % Calculate the number of terms and tensors
    numTerms = height(tbl);
    numTensors = height(T);

    % Extract indices and values from the table
    indices = tbl{:,"i"+T.ID};
    subscripts = num2cell(tbl{:,"v"+ord});
    
    % If there are no tensors, return 1 as the value
    if numTensors == 0
        val = 1;
        return
    end
    
    % Initialize dims based on the order of indices
    dims = zeros(1, length(ord));
    for idx = 1:length(ord)
        tensorIdx = cellfun(@(x)ismember(ord(idx),x),T.legList);
        sz = T.value(tensorIdx).size;
        dims(idx) = sz(find(ord(idx) == T.legList{tensorIdx}, 1));
    end
    
    % Initialize the value with the first tensor's values
    V = T.value(1).val(indices(:,1));
    for idx = 2:numTensors
        % Multiply the values of the tensors together
        V = V .* T.value(idx).val(indices(:,idx));
        % disp("times") % Display a message indicating multiplication
    end
    val=SparseEx;
    val.val=V;
    val.key=tbl{:,"v"+ord};
    val.size=dims;
    val=val.C(level="low");
    % 
    % % Initialize the result array with the appropriate class
    % val = zeros([dims, 1], class(V));
    % coeffIndices = zeros(size(V));
    % 
    % % Calculate the linear indices for the coefficients
    % for idx = 1:numTerms
    %     coeffIndices(idx) = sub2ind([dims, 1], subscripts{idx,:});
    % end
    % 
end
function [T, acc] = calcContractIndex(T)
    % calcContractIndex: calculate index for tensor network to contract
    % acc: remaind index
    
    acc = zeros(1,0);

    % Preallocate result columns
    T.cur = cell(height(T), 1);  % indices in current row that are already in acc
    T.acc = cell(height(T), 1);  % positions in acc matching current
    T.accRank(:)=0;

    for i = 1:height(T)
        cur = T.legList{i};

        % Find contraction positions
        T.cur{i} = find(ismember(cur, acc));
        T.acc{i} =arrayfun(@(x)find(acc==x),cur(T.cur{i}));

        % Remove overlaps
        acc(T.acc{i}) = [];
        cur(T.cur{i}) = [];

        % Add new indices to accumulator
        acc = [acc, cur];
        T.accRank(i) = length(acc);
    end
end
function expr = makeExpressionString(T)
    % Start with the first label as the initial expression
    expr = '1';
    accRank=0;
    % Loop over remaining rows to nest tensorprod calls
    for i = 1:height(T)

        % Convert numeric arrays to string representation
        accStr = mat2str(T.acc{i});
        curStr = mat2str(T.cur{i});
        if isempty(T.acc{i})
            accStr='[]';
            curStr='[]';
        end

        if accRank==0
            expr=[expr '*' T.labels{i}];
            accRank=T.accRank(i);
            continue
        elseif accRank>1
            format='tensorprod(%s,%s,%s,%s)';
        elseif accRank==1
            format='tensorprod(%s,%s,%s,%s,N=1)';
        end
        accRank=T.accRank(i);
        % nest tensorprod
        expr = sprintf(format, expr, T.labels{i}, accStr, curStr);
    end
end
function expr=permDim(expr,acc,ord)
    if ~(length(acc)==length(ord)&&all(ismember(acc,ord)))
        error("Second arguments(order) must be a permutation of remaining tensor indices:%s",mat2str(acc))
    end
    [~, perm] = ismember(ord,acc);

    if isempty(ord)
        perm=[1 2];
    end
    expr=sprintf('permute(%s,%s)',expr,mat2str(perm));
end
function [termTuples, Tjoin] = enumerateTermTuples(T)
    % Enumerates all matching term index tuples [t1 ... tN] from tensors in table T,
    % but only for nonzero elements in T.value. Contracts axes based on matching
    % dummy variable labels in legList.
    %
    % T: table with variables:
    %   - labels (string)         : tensor label (unused here)
    %   - legList (cell array)    : axis dummy variable labels, e.g., [1 2]
    %   - rank (scalar)           : tensor rank (unused)
    %   - value (cell array)      : numeric array for tensor values
    %   - size (cell array)       : numeric vector of axis dimensions
    %
    % Output:
    %   termTuples : (#matches) x N matrix, each row is [t1 ... tN]
    %   Tjoin      : joined table including term_id columns for each tensor

    N = height(T);

    % Collect all unique global dummy variable labels
    allVars = unique(cat(2, T.legList{:}));
    allVars = sort(allVars);

    % Build a table for each tensor containing term_id, contracted variable columns
    X = cell(1,N);
    for i = 1:N
        labels_i = T.legList{i};
        dims_i   = T.size{i};
        assert(numel(labels_i)==numel(dims_i), 'Mismatch between legList and size lengths.');

        % Find nonzero subscripts in value
        nonzeroIdx=find(T.value{i}~=0);
        % 
        [subscripts{1:numel(dims_i)}] = ind2sub(dims_i,nonzeroIdx);
        indexTable=cell2table(subscripts,VariableNames="v"+(1:T.rank(i)));
        X{i} = indexTable;
    end

    % Sequentially inner join tables on shared variable columns
    Tjoin = X{1};
    for i = 2:N
        keys = intersect(Tjoin.Properties.VariableNames, X{i}.Properties.VariableNames);
        keys = keys(startsWith(keys,'v'));
        if isempty(keys)
            error('No common variable labels (keys) between tensor %d and current join.', i);
        end
        Tjoin = innerjoin(Tjoin, X{i}, 'Keys', keys);
    end

    % Extract term index tuples [t1 ... tN]
    termTuples = zeros(height(Tjoin), N);
    for i = 1:N
        col = sprintf('t%d',i);
        termTuples(:,i) = Tjoin.(col);
    end
end
