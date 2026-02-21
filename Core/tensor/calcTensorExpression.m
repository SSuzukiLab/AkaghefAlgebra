%{
# Prompt:
Generate `calcTensorExpression` code from a map definition
You are given a multiline **Premise** describing a helper
`calcTensorExpression(exprPattern, ord, Elem2NonzeroIndices)`.
Your job is: from the **Input** (a linear map written in Sweedler form / basis form), output (i) a **Basis form** with explicit indices and (ii) the **Implementation** line
`phi = calcTensorExpression('...', [...]);`.
You don't need to input 'Elem2NonzeroIndices'.
Use the fixed output format. Do NOT generate a handmade script.
we will explain calcTensorExpression below.
%}
function val=calcTensorExpression(exprPattern,ord,Elem2NonzeroIndices)
    % calcTensorExpression Evaluate symbolic tensor expressions specified as indexed patterns.
%{
    Usage:
        phi = calcTensorExpression(exprPattern, indexOrder)
            exprPattern : (char) Tensor contraction pattern using brace-delimited slots.
            indexOrder  : (vector) Permutation describing the external index ordering.
            phi         : Resulting tensor after applying the specified pattern and reindexing.
    All contraction indices appear exactly twice across the tensors in exprPattern and indexOrder.
    This helper interprets expressions such as 'M{3,4,1}C{2,5,6}S{7,6}' and reorders the
    indices according to the provided permutation, enabling compact specification of
    composite tensor contractions.
      No explicit pairing tensor is introduced.
      Contraction is determined solely by repeated index labels.
    %% example
     Φ : U* ⊗ U → U ⊗ U*
      M{i,j,k} represents M_{j,k}^{i}
      C{i,j,k} represents C_{i}^{j,k}
    Sweedler form: Φ(f ⊗ u)=f(e_{i(2)}) e_{i(1)}*u ⊗ e^i
    Basis form:Φ(e^{i1} ⊗ e_{i2})= Σ_i4  C_{i3}^{i4,i1}  M_{i4,i2}^{i5} e_{i5} ⊗ e^{i3}
    Implementation: phi = calcTensorExpression('C{3,4,1}M{4,2,5}', [1,2,5,3]);

## Hard constraints (must follow)

1. **Index labels** are positive integers. Every contraction index must appear **exactly twice** across `exprPattern` and `ord`.
2. **No explicit pairing tensor** is introduced. If a pairing would normally appear, you must encode it by placing the corresponding index into an existing tensor leg (typically into `C{*,*,fIndex}`) consistent with the given conventions.
3. Tensor conventions (fixed):

   * `M{i,j,k}` represents (M_{j,k}^{i}).
   * `C{i,j,k}` represents (C_{i}^{j,k}).
   * (If present) `S{i,j}` represents (S_{j}^{i}).
4. The output **must be compact** and must contain exactly these two blocks:

   * `Basis form: ...`
   * `Implementation: ...`
5. The external index order `ord` must list the indices in the order the resulting tensor is stored, e.g. `[f,u,uOut,fOut]`.

## Procedure you must follow (deterministic)

A. **Choose external indices first**

* Assign indices to each input slot and output slot in the order they are mentioned in the domain/codomain.
* Default assignment:

  * For a map (U^*\otimes U\to U\otimes U^*): input (f\to 1), input (u\to 2), output (u'\to 5), output (f'\to 3).
  * Use new numbers only if the input already specifies them.
* Set `ord` to the external order as they appear in the final tensor (inputs first, then outputs unless user says otherwise).

* Index–Number Consistency Rule.
    The integer label attached to each typed index (e.g. f1, u2, ft9) must coincide exactly with the numeric label used in exprPattern and ord.
    The number uniquely determines the index identity; the prefix (f,u,ft,ut) is only a type annotation.
    No additional numbering systems are allowed.
B. **Translate the basis form into tensors**

* Write (\Phi(e^{i1}\otimes e_{i2})) as a sum over internal indices (e.g. (\Sigma_{i4})).
* Each coefficient must be a product of tensors `M`, `C`, `S`, etc.
* Use your conventions to decide where each index appears in `M{...}` / `C{...}`.

C. **Build ****`exprPattern`**

* Convert each tensor factor into the brace pattern.
* Ensure each internal index appears exactly twice across the full pattern and is not listed in `ord`.
* Ensure each external index appears exactly once in `ord` and at least once in the pattern.
You don't need to input 'Elem2NonzeroIndices'.

D. **Sanity checks**

* Check: number of free indices equals `length(ord)`.
* Check: every contraction index appears exactly twice (once upper/once lower).
* Check: no tensor is needed between a basis and its dual.
* Check: every numeric label appearing in exprPattern appears with the same number in the Basis form (typed version).

E. **tensor of Hopf algebras**

* M{i,j,k} =M_{j,k}^{\,\, i}, C{i,j,k} =C_{i}^{\, j,k},  S{j,i}=S^{j}_{\, i}

## Output format (exact)

Basis form: 
Implementation: phi = calcTensorExpression('', );

## Worked example (must match)

### Input form
input maybe abstract(you need to infer and deduce this input form)
Φ : U* ⊗ U → U ⊗ U*, 
Sweedler form: Φ(f ⊗ u)=f(e_{i(2)}) e_{i(1)}*u ⊗ e^i

### Output form

Basis form: Φ(e^{i1} ⊗ e_{i2}) = Σ_{i4} C_{i3}^{i4,i1} M_{i4,i2}^{i5} e_{i5} ⊗ e^{i3}
Implementation: phi = calcTensorExpression('C{3,4,1}M{4,2,5}', [1,2,5,3]);

## If you need to improve the input (allowed additions)

If the user's input is ambiguous, you may append a short **Assumptions** clause inside the Basis form line (still single-line), but do not add extra paragraphs. Typical allowed assumptions:

* which side of multiplication corresponds to `*u` vs `u*`;
* whether (\Delta) or (\Delta^{cop}) is used;
* whether (S) or (S^{-1}) is used.

If data structure of tensor is ambiguous, ask to user.

    %}
    arguments
        exprPattern (1,:) char
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
    exprPattern=strip(exprPattern);
    T=parseFormulaArgument(exprPattern);
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
    if true
        T=calcContractOrder(T);
        tbl=calcIndexTable(T,ord);
    else
        joinOrd=calcContractOrder2(T);
        tbl=calcIndexTable2(T,joinOrd);
    end
    val=multiplyCoeff(tbl,T,ord);
end
function T = parseFormulaArgument(arg)
    % Match chemical symbols with corresponding indices
    pattern = '(?<label>[A-Za-z\(\)\^\-\d]+)\{(?<indices>[0-9,]+)\}';
    pattern = '(?<label>[A-Za-z\(\)\^\-\d]+(?:\.[A-Za-z\(\)\^\-\d]+)?)\{(?<indices>[0-9,]+)\}';
    pattern = '(?<label>[A-Za-z_\(\)\^\-\d\.]+)\{(?<indices>[0-9,]*)\}';
    tokens = regexp(arg, pattern, 'names');
    
    labels = string({tokens.label})';
    legList={tokens.indices};
    scalar_idx=cellfun(@isempty,legList);
    legList=cellfun(@(str){str2double(strsplit(str, ','))},legList).';
    legList(scalar_idx)={zeros(1,0)};
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
    % tbl=removevars(tbl, "v"+setdiff(uniqueVars,ord));
    numTensors = height(T);
    if numTensors == 0
        val = 1; % If there are no tensors, return 1 as the value
        return
    end
    
    % Initialize dims based on the order of indices
    dims = zeros(1, length(ord));
    for idx = 1:length(ord)
        tensorIdx = cellfun(@(x)ismember(ord(idx),x),T.legList);
        sz = T.value(tensorIdx).size;
        dims(idx) = sz(find(ord(idx) == T.legList{tensorIdx}, 1));
    end
    val=SparseEx;
    if height(tbl)~=0
        % Initialize the value with the first tensor's values
        indices = tbl{:,"i"+T.ID};
        V = T.value(1).val(indices(:,1));
        for idx = 2:numTensors
            % Multiply the values of the tensors together
            V = V .* T.value(idx).val(indices(:,idx));
            % disp("times") % Display a message indicating multiplication
        end
        val.val=V;
        val.key=tbl{:,"v"+ord};
    else
        val.val=[];
    end
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

function ord = calcContractOrder2(T)
    % Greedy pairwise plan over T.indexTable (current-position indexing).
    A = T.indexTable(:); A = A(~cellfun(@isempty,A));
    n = numel(A); ord = zeros(max(n-1,0),2); if n<=1, return; end
    K = cell(n,1); H = zeros(n,1);
    for i=1:n, [K{i},~]=kv(A{i}); H(i)=height(A{i}); end
    
    for s=1:n-1
        m = numel(K); best = inf; bi=1; bj=2;
        for i=1:m-1
            for j=i+1:m
                sh = numel(intersect(K{i},K{j},'stable'));
                % crude but effective: prefer (i) smaller tables, (ii) more shared keys
                sc = log(max(H(i),1)) + log(max(H(j),1)) - sh*log(2) + 0.2*numel(union(K{i},K{j},'stable'));
                if sc < best, best=sc; bi=i; bj=j; end
            end
        end
        ord(s,:) = [bi bj];
        
        % update metadata after "virtual merge"
        newK = union(K{bi},K{bj},'stable');
        other = {};
        for t=1:m, if t~=bi && t~=bj, other = union(other,K{t},'stable'); end, end
        liveK = setdiff(newK, setdiff(newK,other,'stable'), 'stable'); % (= newK \ internal)
        sh = numel(intersect(K{bi},K{bj},'stable'));
        H(bi) = max(1, round(exp(log(max(H(bi),1))+log(max(H(bj),1)) - sh*log(2))));
        K{bi} = liveK; K(bj)=[]; H(bj)=[];
    end
    
    function [keys,val] = kv(tbl)
        vn = tbl.Properties.VariableNames;
        val = '';
        if any(strcmp(vn,'val')) && isnumeric(tbl.val), val='val'; end
        if isempty(val)
            numv = vn(cellfun(@(x)isnumeric(tbl.(x)),vn));
            if isempty(numv), error('No numeric value column.'); end
            val = numv{end};
        end
        keys = setdiff(vn,{val},'stable');
    end
end

function [tbl,T] = calcIndexTable2(T,joinOrd)
    F = T.indexTable(:);
    F = F(~cellfun(@isempty,F));
    
    for s = 1:size(joinOrd,1) % 進度 s/size(joinOrd,1)
        i = joinOrd(s,1); j = joinOrd(s,2);
        if i > j, [i,j] = deal(j,i); end
        F{i} = joinTable(F{i}, F{j});
        F(j) = [];
    end
    if isempty(F)
        tbl=table;
        return
    end
    tbl = F{1};
end

function tbl=joinTable(T1,T2)
    try
        tbl = innerjoin(T1, T2);
    catch % if no contraction btw T1 and T2, calc kronecker product
        C = combinations(1:height(T1), 1:height(T2));
        C = rowfun(@(x,y)[T1(x,:), T2(y,:)], C);
        tbl = C{:,1};
    end
end