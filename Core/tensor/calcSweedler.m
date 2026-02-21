function [ret,ret2,ret3]= calcSweedler(expr)
    % calcSweedler - Sweedler記法の構文解析と評価を行う関数
    %
    % 構文: [ret,ret2] = calcSweedler(expr)
    %
    % 入力:
    %   expr - Sweedler記法で書かれた文字列（テンソル式）
    %
    % 出力:
    %   ret  - 構文解析木（AST: Abstract Syntax Tree）
    %   ret2 - 解析されたノード情報の配列
    %
    % SweedlerParser Syntax Tree for Sweedler notation
    %{
"I want to first decide the plan.tion [ret,ret2]= calcSweedler(expr)
% SweedlerParser Syntax Tree for Sweedler notationZ
“I want to first decide the plan.
	1.	Analyze the depth of parentheses and parse with * and |
	2.	Parse ev
	3.	Resolve {} indices
	4.	Evaluate values
    5.  Calculate zero objects
    6.  resolve coproduct, find ambient space, pad unit
    7.  calculate product and tensor product
	8.	Obtain tensor expression
    9.  evaluate values
    

構文木の葉の式における添字の規則
[d]:duplicate:複製  x_1,x_2,...      rank:N->N
{m}:comultipl:余積  x_1⊗x_2=1⊗x+x⊗1  rank:1->N
[a]:coaction :余作用  x_0⊗g_1=x⊗g    rank:1->N
{i}:indexing :分割アクセス (x⊗y)_1=x,     rank:N->1
{p}:placement:配置 x_2=1⊗x,  R_{2,3} 多重ok  rank:N->N, ambient space内に配置



想定パターン(先頭[d]は除く)
[a], {p}, {i}, {m}, {i}{m}, {i}[a]
Nd(添字の深さ)はそれぞれ，2,2,2,2,3,3
最終パターン
dim:[d]{i}{m}, dia:[d]{i}[a], dp:[d]{p}

{args}の解釈パターン:
1. if start with {            case: add [d]  d_fill=1or2
2. rank>1 & length(args{1})>1 case: {p} multi-index ,d_fill=1
3. Nd=1                       case: {i}{m}で初期化 ,d_fill=1
4. Nd=3                       case: {i} 添字付け {i}{m} or {i}[a],d_fill=2
5. rank=1 & kind(2)='['       case: [a] 余作用と解釈 ,d_fill=0以外は禁則
6. rank=1 & length(args{1})=1 case: {m} 一旦余積と解釈 , d_fill=2
7. rank=1 if error in 2       case: {p} エラーしたら配置と解釈 ,d_fill=1

禁止パターン
[a]:duplicateで始まると混乱 [d][a]にする

d_fill in {0,1,2}:d_fill=0: [d] is specified, 
d_fill=1: fill [d] with 1,2,... 
d_fill=2: fill [d] with 1


Please suggest method names for each.”
sweedler notationには複数の意味が混在していて，文脈により解釈が別れるため，解析が非常に難しい．
特に，idをパディングするのか，余積をとるのかの記法との判別をつける方法を探したい．
例えばR行列など複数同じテンソルを用いる場合は，コピーの何個目かを判別するための添字が必要となる．
そのための手段として，[n]をつけることでコピー要素の判別を行う．
R行列は一括して評価が行われるが，テンソルの式の評価の際に，コピーされて評価されるようにする．
五角関係式ややんバグスター方程式などではコピーする数分だけ変数を用意するのは面倒である．
括弧の種類を区別できたらいいが，通常使うことができる括弧のかずが[].{},(),<>であり，()は数式の優先順位として用いたい，
<>はevaluationで用いたい　という要望から，使える括弧が{}と[]だけになり，括弧の種類が足りていないということになる．
この問題をどうにかできる記法は存在しないだろうか？
文脈で判断させるためには，一回その式の評価を行う必要があり，一貫して同じ評価を行える保証は薄くなる．　これはリスキーである．

[]はコピーの識別
% Example:
%   s = 'fa{1}*ev(gb{1}, ev(S^-1,fa{2}[3])*?*fa{2}[1]) | fa{2}[2]*gb{2}';
output='fa{1,2}Delta{3,1,4}gb{7}Sinv{5,8}Delta{,fa{2}[3])mu{5,4,6}?fa{2}[1]) | fa{2}[2]*gb{2}';
⸻
'fa{1}*<gb{1}, Sn(fa{2}{3})*?*fa{2}{1},-1)> | fa{2}{2}*gb{2}';

pentagon equation
'W[1]{1,2}*W[2]{1,3}*W[3]{2,3}=W[4]{2,3}*W[5]{1,2}'

pentagon equation(coherence)
'<Delta[1]{3},Phi[1]>*<Delta[2]{3},Phi[2]>'
'Phi[1]{2,3,4}*<Delta[1]{2},Phi[2]>*Phi[3]{1,2,3}'

'x{212}'(quasi hopf algebra coproduct)

Here’s a clean naming scheme (MATLAB-style, consistent, verbs for actions):
	1.	parseMulTensor – parse by * and | using parentheses depth.
	2.	parseEv – detect and build Ev nodes from ev(…, …).
	3.	parseIndices – resolve Sweedler indices {…} into integer arrays.
	4.	evalNode – evaluate an AST node to a concrete algebraic object.
	5.	toTensorExpr – convert the evaluated node into a structured tensor expression.
Each function has a clear purpose and consistent naming, making the overall parser easy to understand and maintain.
    %}


    arguments
        expr (1,:) char  % 入力式は1行の文字列でなければならない
    end

    % ステップ1: 空白文字を除去して式を正規化
    expr(isspace(expr))=[]; % remove spaces

    % ステップ2: 式を構文解析して構文木（AST）を構築
    st = parseExpr(expr,0);
    % ステップ3: インデックス情報を解析し、各ノードに添字情報を付加
    [st,nodes]= parseIndices(st,0);

    % ステップ4: ノードから式（atom）を抽出し、ユニークな式のリストを作成
    atoms=unique(string({nodes.atom}));
    values=cell(size(atoms));                   
    for i=1:length(atoms)
        if ~startsWith(atoms(i),'?')
            % 呼び出し元のワークスペースで変数を評価
            values{i}=evalin("caller",atoms(i));
        else
            % '?'で始まる場合は未評価の変数として扱う
            values{i}=atoms(i);
        end
    end
    % ステップ5: インデックス情報を評価（未実装）
    [st,nodesU]=evaluateIndices(st,nodes,atoms,values);
    % ステップ6: 型情報を評価（未実装）
    %構文木の各分岐に対して，実行時に決まる型の情報をzero objとして与える．
    % 
    info=struct(dualzero=nan);
    dict_num2nodeU=getDictOfGroup({nodesU.num});
    st=evaluateTypes(st,nodesU,dict_num2nodeU,info);
    

    % 戻り値を設定
    ret=st;      % 構文木
    ret2=nodesU;  % ノード情報配列
    ret3=nodes;
    % ret2=nodes;  % ノード情報配列
end
function ret=getDictOfGroup(groupCell)
    % getDictOfGroup - ユニークなグループ要素の辞書を作成
    ret=dictionary();
    for i=1:length(groupCell)
        ret(groupCell{i})=i;
    end
end
function st=evaluateTypes(st,nodesU,dict_num2nodeU,info)
    % evaluateTypes - 構文木のノードに型情報を付加する
    % 各分岐に対して，実行時に決まる型の情報をzero objとして与える．
    % ?のノードはpairingのdual objの型情報からもとの情報を復元する操作を経る
    % 入力:
    %   st     - 構文木のノード
    %   nodesU - ノード情報の配列（インデックスが評価されたもの）
    %   info - 追加情報の構造体 ?の型を決めるためにdualzeroを与える
    % 出力:
    %   st     - 型情報が評価された構文木
    %   st.zeroを追加する
    switch st.type
        case 'leaf'
            idx=dict_num2nodeU(st.num); % index of nodesU
            type=nodesU(idx).type;
            switch type
                case  'sp'
                    % not impl
                case 'vec'
                    val=nodesU(idx).value;
                    switch nodesU(idx).itype
                        case 'dim'
                            
                        case 'dp'
                            % not impl
                        case 'dia'
                            % not impl  
                    end

                case 'num'
                    % not impl
            end
            

            
        case 'fun'
            2
        otherwise
            st.left=evaluateTypes(st.left,nodesU,dict_num2nodeU,info);
            st.right=evaluateTypes(st.right,nodesU,dict_num2nodeU,info);
    end

end
function [st,nodesU]=evaluateIndices(st, nodes,atoms, values)
    % evaluateIndices - インデックス情報を評価する関数（未実装）
    %
    % 入力:
    %   st     - 構文木のノード
    %   nodes  - ノード情報の配列
    %   atoms  - ユニークな式のリスト
    %   values  - 各式に対応する値のセル配列
    %
    % 出力:
    %   st     - インデックスが評価された構文木
    %
    % 現在は未実装。将来的にインデックス評価のロジックを追加予定。
    atoms_=string({nodes.atom});  % 全ノードの式を文字列配列に変換
    type=strings(size(atoms));
    itype=cell(size(atoms));
    rank=zeros(size(atoms));
    nodesU=repmat(struct(itype='',index=[],rank=nan,num={}),size(atoms));
    
    for i=1:length(atoms)
        % '?'で始まる式は変数として扱う、それ以外は呼び出し元の環境で評価
        if ~startsWith(atoms(i),'?')
            % 評価結果がSparseExオブジェクトの場合
            if isa(values{i},'SparseEx')
                type(i)="sp";              % 型を"sp"（スパース）に設定
                rank(i)=values{i}.rank;     % テンソルのランクを取得
            elseif isa(values{i},'VectAlg')
                type(i)="vec";
                rank(i)=values{i}.rank;
            else
                type(i)="num";          % それ以外はスカラーとみなす
                rank(i)=0;                 % スカラーの場合はランク0
            end
        else
            % '?'で始まる場合は未評価の変数として扱う
            type(i)="var";
        end
        % この式が使われている全ての箇所のインデックス情報を収集
        idx_atoms = ismember(atoms_,atoms(i));
        Ilist={nodes(idx_atoms).index};
        kindsC={nodes(idx_atoms).kind};
        nums=[nodes(idx_atoms).num];
        Nc=length(Ilist);
        depth=cellfun(@length,kindsC);
        Nd=max(depth);
        Istr=strings(Nc,Nd);
        kinds=char(size(Istr));
        itype{i}='dim';
        for j=1:Nc
            Istr(j,1:depth(j))=string(Ilist{j});
            kinds(j,1:depth(j))=[kindsC{j}{:}];
        end
        if Nd~=0&&any(depth==0)
            error('Inconsistent index depth for atom %s', atoms(i));
        end
        
        if Nd==0||kinds(1,1)~='['
            assert(Nd==0||all(kinds(:,1)~='['), ...
                'Inconsistent copy-index for atom %s', atoms(i))
            kinds=[repmat('[',Nc,1),kinds];
            Istr=[strings(size(Istr,1),1),Istr];
            Nd=Nd+1;
            % Istr(:,1)=string(1:size(Istr,1));
        else
            d_fill=0;
        end
        idxLen=arrayfun(@(x)length(str2num(x)),Istr);
        if any(idxLen(:)>1)
            % {p} placement multi-index
            is_multi_arg=any(idxLen>1   ,1);
            assert(length(is_multi_arg)==2&&isequal(is_multi_arg,[false,true]), ...
                'Currently only support one multi-index at atom %s', atoms(i)); 
            is_multi_valid_rank=rank(i)==idxLen(:,2);
            assert(is_multi_valid_rank, ...
                'multi-index must be equal to rank for atom %s', atoms(i));
            itype{i}='dp';
            d_fill=1;
        elseif Nd==1
            % {i}{m}で初期化
            % assert(rank(i)==0, ...
            %     'rank must be 0 for Nd=1 at atom %s', atoms(i));

            itype{i}='dim';
            d_fill=1;
        elseif Nd==2
            itype{i}='dim';
            d_fill=2;
            Istr(:,3)="1";
        elseif Nd==3
            % {i} 添字付け [d]{i}{m} or [d]{i}[a]
            assert(all(all(kinds==kinds(1,:)|kinds==0)), ...
                'Inconsistent index kinds for atom %s', atoms(i));
            if kinds(1,3)=='['
                itype{i}='dia';
            elseif kinds(1,3)=='{'
                itype{i}='dim';
            else
                error('Invalid index kinds for atom %s', atoms(i));
            end
            d_fill=2;
        elseif rank(i)==1&size(Istr,2)==2 
            % [d]{m}, [d][a] or [d]{p} 
            if kinds(1,2)=='['
                itype{i}='dia';
                Istr(:,3)=Istr(:,2);
                Istr(:,2)="1";
                d_fill=0;
            elseif kinds(1,2)=='{'
                itype{i}='dm';
                try % [d]{m} -->[d]{i}{m}
                    itype{i}='dim';
                    C=getMor(values{i},"coprod");
                    Istr(:,3)=Istr(:,2);
                    Istr(:,2)="1";
                    d_fill=2;
                catch ME % [d]{p} -->[d]{i}{p}
                    itype{i}='dp';
                    d_fill=1;              
                end
            end
        else 
            
        end
        if d_fill==1
            Istr(:,1)=string(1:Nc);
        elseif d_fill==2
            Istr(:,1)="1";
        end
        if itype{i}(2)=='p'
            indexA=zeros(1,1+rank(i));
            for j=1:Nc
                indexA(:,1)=str2num(Istr(j,1));
                indexA(:,2:end)=str2num(Istr(j,2));
            end
        else
            indexA=arrayfun(@str2num,Istr);
        end
        
        nodesU(i).index=indexA;
        nodesU(i).itype=itype{i};
        nodesU(i).rank=rank(i);
        nodesU(i).value=values{i};
        nodesU(i).type=char(type(i));
        nodesU(i).num=nums;
        
    end

    % 将来的な拡張: テンソル式への変換（現在はコメントアウト）
    % ret=toTensorExpr(st,atoms,value,0);

end

function [st, nodes] = parseIndices(st,id0)
    % parseIndices - 構文木のノードからインデックス情報を抽出・解析
    %
    % 入力:
    %   st  - 構文木のノード
    %   id0 - ノードID初期値（再帰呼び出し時に使用）
    %
    % 出力:
    %   st    - インデックス情報が付加された構文木
    %   nodes - 全リーフノードの配列（フラット化された情報）
    %
    % 例:
    %   'fa[3]{2}{3}' -> atom='fa', index={3,2,3}, kind ={'[','{','{'}
    %   'W[3]{1,2}'   -> atom='W', index={3,[1,2]}, kind ={'[','{'}

    if strcmp(st.type,'leaf')  % リーフノード（終端記号）の場合
        dic=dictionary('{' ,'}','[',']');
        s = st.expr;
        % Base atom until first { or [
        pos = find(ismember(s,'{['),1);
        if isempty(pos), pos = length(s)+1; end
        st.atom = s(1:pos-1);
        st.index={};
        st.kind={};
        rest = s(pos:end);
        while ~isempty(rest)
            kind_c = rest(1);
            % Vectorized extraction of all {} and [] groups in one pass
            try
                pos=find(rest==char(dic(kind_c)),1);
            catch
                error('invalid parenthesis in "%s"', s);
            end
            assert(~isempty(pos), 'Unmatched %s in "%s"', kind_c, s);
            s=rest(2:pos-1);
            if contains(s,',')
                s=['[',s,']'];
            end
            st.index{end+1} =s;
            st.kind{end+1} = kind_c;
            rest=rest(pos+1:end);
        end
        nodes=st;  % リーフノードの場合はノード自身を返す
    elseif strcmp(st.type,'fun')  % 関数呼び出しノードの場合
        % 関数呼び出しノードの場合: 各引数を再帰的に処理
        nodes = struct([])
        for i=1:length(st.child)
            [st.child{i}, nodes_i] = parseIndices(st.child{i}, id0 + length(nodes));
            nodes = [nodes, nodes_i];  % 各引数のノードを結合
        end
    else
        % 演算子ノード（非終端記号）の場合: 左右の部分木を再帰的に処理
        [st.left,nodes_l] = parseIndices(st.left,id0);  % 左部分木を処理
        [st.right,nodes_r] = parseIndices(st.right,id0+length(nodes_l));  % 右部分木を処理（IDオフセット付き）
        nodes = [nodes_l,nodes_r];  % 左右のノードを結合してフラットな配列として返す
    end
end

function [node,num] = parseExpr(s,num)
    % parseExpr - 式を構文解析して構文木（AST）を構築
    %
    % 入力:
    %   s - 解析する式（文字列）
    %
    % 出力:
    %   node - 構文木のノード
    %
    % 優先順位（低い順）:
    %   1. '|' （テンソル積）
    %   2. '*' （積）
    %   3. '<...>' （評価演算）
    %   4. 関数呼び出し形式 "fun(arg1, arg2, ...)"
    %   5. リーフ（式）

    s =stripOuterParens(s);  % 外側の不要な括弧を除去
    if isempty(s), error('Empty expression at character'); end

    % 優先順位1: 最も優先度が低い '|' 演算子で分割（右結合）
    k = findTopLevelOp(s, '|');
    num=num+1;
    node=struct(expr=s,num=num);  % 元の式を記録
    if k > 0
        node.type='op';
        node.oper='|';
        [node.left,num]  = parseExpr(s(1:k-1),num);
        [node.right,num] = parseExpr(s(k+1:end),num);
        return;
    end
    % 2) next: split on top-level '*', rightmost for left-assoc
    k = findTopLevelOp(s, '*');
    if k > 0
        node.type = 'op';
        node.oper = '*';
        [node.left,num]  = parseExpr(s(1:k-1),num);
        [node.right,num] = parseExpr(s(k+1:end),num);
        return;
    end
    % 3) next: parse ev( ... )
    k=findEval(s);
    if k>0
        node.type='ev';
        [node.left,num] = parseExpr(s(2:k-1),num);
        [node.right,num] = parseExpr(s(k+1:end-1),num);
        return;
    end
    [tf,fname,arg]=findFunCall(s);
    if tf
        node.fname=fname;
        node.type='fun';
        node.arg=arg;
        child = cell(size(arg));
        for i = 1:numel(arg)
            [child{i}, num] = parseExpr(arg{i}, num);
        end
        node.child=child;
        return;
    end
    % 4) leaf
    node = struct(type='leaf',expr=s,num=num);
end

function k = findTopLevelOp(s, op)
    % findTopLevelOp - トップレベル（括弧の外側）にある中置の二項演算子の位置を検索
    %
    % 入力:
    %   s  - 検索対象の文字列
    %   op - 検索する演算子（'|' または '*'）
    %
    % 出力:
    %   k - 最も右にあるトップレベルの演算子の位置（見つからない場合は0）
    %
    % 右結合のため、最も右にある演算子を返す

    idxOp=strfind(s,op);  % 演算子が出現する全ての位置を取得
    % 括弧の深さを計算（開き括弧で+1, 閉じ括弧で-1）
    depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
    % 括弧のバランスチェック（最終的に深さは0でなければならない）
    assert(depth(end)==0,'Unbalanced parentheses in "%s"', s);
    % 深さ0（トップレベル）にある演算子のみを抽出
    idxTopOp=idxOp(depth(idxOp)==0);
    if isempty(idxTopOp), k=0; return; end  % トップレベルに演算子がない
    k = idxTopOp(end);  % 最も右にある演算子の位置を返す（右結合）
end
function k = findEval(s)
    % findEval - 評価演算 <expr1, expr2> の形式を検出
    %
    % 入力:
    %   s - 検索対象の文字列
    %
    % 出力:
    %   k - カンマの位置（評価演算でない場合は0）
    %
    % '<expr1, expr2>' 形式の場合、カンマの位置を返す

    k=0;
    depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
    % 括弧のバランスチェック
    assert(depth(end)==0,'Unbalanced parentheses in "%s"', s);
    if ~(s(1)=='<'&& s(end)=='>')
        return; % '<' で始まり '>' で終わらない場合は評価演算ではない
    end
    k=find(s==','& depth==1); % 深さ1（最外側の<>内）にあるカンマを検索
    assert(isscalar(k),'syntax error of eval map <arg1,arg2> in "%s"', s);
end
function [tf, fname, arg] = findFunCall(s)
    % findFunCall - 関数呼び出し形式 "fun(arg1, arg2, ...)" を検出
    %
    % 入力:
    %   s - 検索対象の文字列
    %
    % 出力:
    %   tf    - 関数呼び出しの場合true
    %   fname - 関数名
    %   arg   - 引数のセル配列

    arguments
        s (1,:) char
    end

    tf = false; fname = ''; arg = '';

    % '(' の位置を検索
    [~,pos_f]=ismember('(',s);
    if pos_f==0, return; end  % '(' がない場合は関数呼び出しではない
    tf=true;
    fname=s(1:pos_f-1); % 関数名を抽出
    depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
    assert(depth(end)==0||all(depth(pos_f:end-1)>0),'Unbalanced parentheses in "%s"', s);
    pos=[pos_f,find(depth==1&s==','),length(s)]; % 引数の区切り位置を検出
    arg=cell(1,length(pos)-1);
    for i=1:length(arg) % 各引数を抽出
        arg{i}=s(pos(i)+1:pos(i+1)-1);  % カンマ間の文字列を引数として抽出
    end
end

function tf = isValidIdent(name)
    % isValidIdent - MATLAB形式の有効な識別子かどうかをチェック
    %
    % 入力:
    %   name - チェックする文字列
    %
    % 出力:
    %   tf - 有効な識別子の場合true
    %
    % 規則: [A-Za-z_][A-Za-z0-9_]* （最初は文字またはアンダースコア、以降は英数字またはアンダースコア）

    if isempty(name), tf = false; return; end  % 空文字列は無効

    % 最初の文字は英字またはアンダースコアでなければならない
    % 2文字目以降は英数字またはアンダースコア
    tf = (isletter(name(1)) || name(1) == '_') && all(isstrprop(name,'alphanum') | name == '_');
end
function s = stripOuterParens(s)
    % stripOuterParens - 式全体を囲む外側の括弧を除去
    %
    % 入力:
    %   s - 処理対象の文字列
    %
    % 出力:
    %   s - 外側の括弧が除去された文字列
    %
    % 例: '((a+b))' -> 'a+b'
    %     '(a)+(b)' -> '(a)+(b)' （全体を囲んでいないので除去しない）

    % 最も外側の括弧が式全体を囲んでいる場合のみ除去（繰り返し適用）
    while numel(s) >= 2 && s(1) == '(' && s(end) == ')'
        % 括弧の深さを計算
        depth=cumsum(ismember(s,'<{([') - ismember(s,'>})]'));
        depth(end)=1;  % 最後の文字は閉じ括弧なので深さを1に設定

        % すべての位置で深さが正（括弧が式全体を囲んでいる）ならば除去
        if all(depth > 0)
            s = s(2:end-1);  % 最初と最後の文字（括弧）を除去
        end
    end
end

function assertNonEmpty(part, s, k, op)
    % assertNonEmpty - オペランドが空でないことを確認
    %
    % 入力:
    %   part - チェックするオペランド
    %   s    - 元の式
    %   k    - 演算子の位置
    %   op   - 演算子
    %
    % 空のオペランドが検出された場合はエラーを発生

    if isempty(part)
        error('Empty operand around operator %s at position %d in "%s".', op, k, s);
    end
end