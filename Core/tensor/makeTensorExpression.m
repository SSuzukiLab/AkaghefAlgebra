function [expression,acc]=makeTensorExpression(str,ord)
    arguments
        str (1,:) char
        ord (1,:) double=[]
    end
    T=parseFormulaArgument(str);
    validateParsedTable(T)
    [T,acc]=calcContractIndex(T);
    expression=makeExpressionString(T);
    expression=permDim(expression,acc,ord);
end

function T = parseFormulaArgument(arg)
    % Match chemical symbols with corresponding indices
    pattern = '(?<label>[A-Za-z\(\)\^\-\d]+)\{(?<indices>[0-9,]+)\}';
    tokens = regexp(arg, pattern, 'names');

    labels = string({tokens.label})';
    indexList={tokens.indices};
    indexList=cellfun(@(str){str2double(strsplit(str, ','))},indexList).';


    T = table(labels, indexList);
    T.rank=cellfun(@length,T.indexList);
end
function validateParsedTable(T)
    % Validate the parsed table where each label is an expression
    % (not just a variable name), which must be evaluable in the caller workspace,
    % and its rank (i.e., ndims) must equal rank(i)
    idx=horzcat(T.indexList{:});
    [vals, ~, grp] = unique(idx);           % Get unique values and group indices
    counts = accumarray(grp, 1)>2;          % Count occurrences
    if any(counts)
        error("contraction indices:%s must appear just 2 times",mat2str(vals(counts)))
    end
    for i = 1:height(T)
        expr = T.labels(i);
        indices = T.indexList{i};

        % Try evaluating the expression in the caller workspace
        try
            value = evalin('base', expr);
        catch
            warning('Expression "%s" cannot be evaluated in the base workspace.', expr);
        end

        % Check that the number of dimensions matches the index list length
        actualRank = ndims(value);
        if isvector(value)
            actualRank=1;
        end
        expectedRank = numel(indices);
        if actualRank ~= expectedRank
            warning('Expression "%s" has rank %d, expected %d.', expr, actualRank, expectedRank);
        end
    end
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
        cur = T.indexList{i};

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