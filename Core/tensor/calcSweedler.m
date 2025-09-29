function [ret,ret2]= calcSweedler(expr)
% SweedlerParser Syntax Tree for Sweedler notation
%{
“I want to first decide the plan.
	1.	Analyze the depth of parentheses and parse with * and |
	2.	Parse ev
	3.	Resolve {} indices
	4.	Evaluate values
	5.	Obtain tensor expression

Please suggest method names for each.”

% Example:
%   s = 'fa{1}*ev(gb{1}, ev(S^-1,fa{2}[3])*?*fa{2}[1]) | fa{2}[2]*gb{2}';
output='fa{1,2}Delta{3,1,4}gb{7}Sinv{5,8}Delta{,fa{2}[3])mu{5,4,6}?fa{2}[1]) | fa{2}[2]*gb{2}';
⸻

Here’s a clean naming scheme (MATLAB-style, consistent, verbs for actions):
	1.	parseMulTensor – parse by * and | using parentheses depth.
	2.	parseEv – detect and build Ev nodes from ev(…, …).
	3.	parseIndices – resolve Sweedler indices {…} into integer arrays.
	4.	evalNode – evaluate an AST node to a concrete algebraic object.
	5.	toTensorExpr – convert the evaluated node into a structured tensor expression.
Each function has a clear purpose and consistent naming, making the overall parser easy to understand and maintain.
%}

    arguments
        expr (1,:) char
    end
    expr(isspace(expr))=[];
    st = parseExpr(expr);
    [st,nodes]= parseIndices(st,0);
    atoms_=string({nodes.atom});
    atoms=unique(atoms_);
    value=cell(size(atoms));
    type=strings(size(atoms));
    rank=zeros(size(atoms));
    idxList1=cell(size(atoms));
    idxList2=cell(size(atoms));
    for i=1:length(atoms)
        if ~startsWith(atoms(i),'?')
            value{i}=evalin("caller",atoms(i));
            if isa(value{i},'SparseEx')
                type(i)="sp";
                rank(i)=value{i}.rank;
            else
                rank(i)=0;
            end
        else
            value{i}=atoms(i);
            type(i)="var";
        end
        idxList1{i}=[nodes(ismember(atoms_,atoms(i))).index1];
        idxList2{i}=[nodes(ismember(atoms_,atoms(i))).index2];
        if rank(i)>1
            assert(any(isnan(idxList1{i})) || numel(unique(idxList1{i}))==rank(i), ...
            "invalid index for "+atoms(i))
        end
    end
    % ret=toTensorExpr(st,atoms,value,0);

    ret=st;
    ret2=nodes;
end
function [st, nodes] = parseIndices(st,id0)
    if strcmp(st.type,'leaf')
        % if strcmp(st.value,'?')
        %     nodes=st;
        %     return
        % end
        % Extract indices from a string like "fa{1}{2}" -> fa,1,2
        st.id=id0+1;
        s=st.value;
        % 1) extract leading identifier (function/variable name)
        idx_s=[find(ismember(s,'{([')),length(s)+1]; 
        [st.type1,st.type2]=deal(''); 
        [st.index1,st.index2]=deal(nan); 
        st.atom=s(1:idx_s(1)-1);
        for i=1:length(idx_s)-1
            st.("index"+i)=str2double(s(idx_s(i)+1:idx_s(i+1)-2));
            st.("type"+i)=s(idx_s(i));
        end
        nodes=st;
    else
        [st.left,nodes_l] = parseIndices(st.left,id0);
        [st.right,nodes_r] = parseIndices(st.right,id0+length(nodes_l));
        nodes = [nodes_l,nodes_r];
    end
end
function node = parseExpr(s)
    s =stripOuterParens(s);
    if isempty(s), error('Empty expression at character'); end
    % 1) lowest precedence: split on top-level '|', rightmost for left-assoc
    k = findTopLevelOp(s, '|');
    node=struct(expr=s);
    if k > 0
        node.type='op';
        node.oper='|';
        node.left  = parseExpr(s(1:k-1));
        node.right = parseExpr(s(k+1:end));
        return;
    end
    % 2) next: split on top-level '*', rightmost for left-assoc
    k = findTopLevelOp(s, '*');
    if k > 0
        node.type = 'op';
        node.oper = '*';
        node.left  = parseExpr(s(1:k-1));
        node.right = parseExpr(s(k+1:end));
        return;
    end
    % 3) next: parse ev( ... )
    k=findEval(s);
    if k>0
        node.type='ev';
        node.left=parseExpr(s(2:k-1));
        node.right=parseExpr(s(k+1:end-1));
        return;
    end
    % 4) leaf
    node = struct('type','leaf','value',s);
end

function k = findTopLevelOp(s, op)
% Return the index of the RIGHTMOST occurrence of op at depth==0, else 0.
    idxOp=strfind(s,op);
    depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
    assert(depth(end)==0,'Unbalanced parentheses in "%s"', s);
    idxTopOp=idxOp(depth(idxOp)==0);
    if isempty(idxTopOp), k=0; return; end
    k = idxTopOp(end);
end
function k = findEval(s)
% Return the index of the RIGHTMOST occurrence of op at depth==0, else 0.
    k=0;
    depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
    assert(depth(end)==0,'Unbalanced parentheses in "%s"', s);
    if ~(s(1)=='<'&& s(end)=='>')
        return;
    end
    k=find(s==','& depth==1);
    assert(isscalar(k),'syntax error of ev(...) in "%s"', s);
end

function s = stripOuterParens(s)
% Strip one pair of outer parentheses if they enclose the ENTIRE expression.
    if numel(s) >= 2 && s(1) == '(' && s(end) == ')'
        depth=cumsum(ismember(s,'<{(') - ismember(s,'>})'));
        depth(end)=1;
        if all(depth > 0)
            s = s(2:end-1);
        end
    end
end

function assertNonEmpty(part, s, k, op)
    if isempty(part)
        error('Empty operand around operator %s at position %d in "%s".', op, k, s);
    end
end