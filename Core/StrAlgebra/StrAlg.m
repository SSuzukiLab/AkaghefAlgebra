classdef(InferiorClasses=?sym) StrAlg<IAdditive&matlab.mixin.Heterogeneous
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        cf (:,1)=0
        pw (:,:) cell ={[]}
        % bs (:,:) cell ={Bases.empty}

        ctype NumericType ="D"
        % ZEROにbaseが依存する形式のほうがよくないか？
        base (1,:) Bases
        spec (1,1) SpaceSpec % specify principle of space
        % also, store structure constant 
        ZERO (1,:) StrAlg
        sortedFlag
        inverseSimplifiedFlag
    end
    properties(Dependent)
        term
        deg
        basestr
        dims
        rank

    end
    % properties(Constant,Hidden)
    %     BASE0=Bases(10,"X","default_base")
    % end


    %% binomial operation
    methods
        function [i1,i2]=alignNum(i1,i2)
            % 型をStrAlgにする
            tf=[isnumeric(i1)||isa(i1,"sym") isnumeric(i2)||isa(i2,"sym")];
            if isequal(tf,[1 0])
                i1=casttype(i2,i1);
            elseif isequal(tf,[0 1])
                i2=casttype(i1,i2);
            end
        end
        function ret=casttype(obj,arg)
            if ~isequal(class(obj),class(arg))
                ret=obj.unit;
                ret.cf=arg;
                % assert(isa(arg,ret.ctype))
            end
        end


        function ret=plus(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            assert(i1.rank==i2.rank,'different rank addition error')
            try
                sz=size(zeros(size(i1))+zeros(size(i2)));
            catch
                error("PolAlg:plus","size dimensions must match")
            end
            if any(sz==0)
                ret=repmat(i1,sz);
            else
                % assert(configure(i1,i2))
                ret=arrayfun(@plus_,repmat(i1,sz./size(i1)),repmat(i2,sz./size(i2)));
            end
            % end

            function ret=plus_(i1,i2)
                % PLUS_ supplementary method for scalar
                ret=i1.set_cp([i1.cf;i2.cf],[i1.pw;i2.pw]).calc();
            end
        end
        % function ret=minus(i1,i2)
        %     ret=i1+(-i2);
        % end
        function i1=uminus(i1)
            i1.cf=-i1.cf;
        end

        % function i1=uplus(i1)
        % end
        function ret=eq(i1,i2)
            sub=calc(i1-i2);
            if sub.ctype=="S"
                ret=fold(@and,sub.cf==0.',symtrue);
            else
                ret=all(sub.cf==0);
            end
        end
        %乗算,作用
        function ret=mtimes(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            assert(isequal(i1.base,i2.base),'異なる空間での積エラー')
            R=i1.rank;
            ret=lfun_(i1|i2,@fun);
            ret.base=i1.base;
            ret.ZERO=i1.ZERO;
            ret=ret.calc();
            function [c,p,b]=fun(p,b)
                c=1;
                % p=cellfun(@(P)[p{1} P],i2.pw,UniformOutput=false);
                % b=cellfun(@(B)[b{1} B],i2.bs,UniformOutput=false);
                p=cellfun(@(p1,p2)[p1 p2],p(1:R),p(R+1:end),UniformOutput=false);
                % b=cellfun(@(b1,b2)[b1 b2],b(1:R),b(R+1:end),UniformOutput=false);
            end
        end
        function ret=lb(i1,i2)
            ret=i1*i2-i2*i1;
        end
        %除算
        function [q,r]=mrdivide(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            [q,r]=helpdiv(0,i1,i2);
        end
        %べき
        function ret=mpower(i1,i2)
            %i2がsymの場合の実装をする
            if i1.term>1
                ret=i1.unit;
                for ii=1:i2
                    ret=i1*ret;
                end
            elseif i1.term==1
                % 実装効率に振ってよくない
                if isa(i2,'sym')
                    if ~isempty(symvar(i2))
                        i1.pw=i1.pw*i2;
                        i1.cf=i1.cf^i2;
                        ret=i1;
                        return
                    else
                        i2=double(i2);
                    end
                end
                pw=cellfun(@(pw)repmat(pw,1,i2),i1.pw,UniformOutput=false);
                ret=i1.set_cp(i1.cf^i2,pw);
            end
        end
        function ret=unit(arg,rank)
            if nargin==1, rank=arg.rank; end
            N=rank/arg.rank;
            mustBeInteger(N)
            ret=arg.set_cp(1,cell(1,rank));
            
            % ZEROを設定するときにsetZEROメソッドを呼ぶほうがよい？
            % ret.ZERO=repmat(arg.ZERO,1,N);
            ret.base=repmat(arg.base,1,N);
        end
        %テンソル積
        function o=or(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            o=i1.lfun_(@fun);
            o.base=[i1.base i2.base];
            o.spec=or(i1.spec,i2.spec);
            % o.ZERO=o.ZERO.set_cp(0,cell(1,o.rank));
            % o.ZERO.ZERO=[i1.ZERO i2.ZERO];
            o.ctype=getType([i1.ctype i2.ctype]);
            function [c,p]=fun(p)
                c=i2.cf;
                p=[repmat(p,i2.term,1) i2.pw];
            end
        end
        %副積
        function o=and(i1,i2)
            error('not implemented')
        end
        %副積計算
        function o=not(i1)
            error('not implemented')
        end

        function o=prodeval(arg)
            r=arg.rank/2;
            mustBeInteger(r)
            o=arg.set_cp(arg.cf, ...
                cellfun(@horzcat,arg.pw(:,1:r),arg.pw(:,r+1:end),UniformOutput=false));
            o.base=arg.base(1:r);
        end
        function o=commeval(arg)
            arg=arg.prodeval;
            [pw,idx]=cellfun(@sort,arg.pw,UniformOutput=false);
            o=arg.set_cp(arg.cf,pw);
        end
        function o=commtimes(i1,i2)
            o=commeval(i1|i2);
        end
        function o=split(i1)
            o=repmat(i1.unit,i1.term,i1.rank);
            for i=1:i1.term
                for j=1:i1.rank
                    o(i,j)=i1.set_cp(1,i1.pw(i,j));
                    o(i,j).base=i1.base(j);
                    % o(i,j).ZERO=i1.ZERO(j);
                end
                o(i,1).cf=i1.cf(i);
            end
        end
        %% 作用,表現
        function ret=repMono(obj)
            error('not implemented')
        end
        function ret=act(obj,vec)
            arguments
                obj StrAlg
                vec PolAlg
            end
            dimA=rep(obj).dimV;
            dimV=vec.dim;
            Nfactor=dimV/dimA;
            mustBeInteger(Nfactor)
            if Nfactor==1
                ret=act(repMono(obj),vec);
            else
                ret=act(rep(DeltaN(obj,Nfactor)),vec);
            end
        end
        function result = rep(obj)
            % Split the object into parts
            splitParts = obj.split;

            % Replace each part using the repMono function
            replacedParts = arrayfun(@repMono, splitParts);

            % Combine the results using a logical OR operation
            combinedResult = replacedParts(:, 1);
            for colIndex = 2:size(replacedParts, 2)
                combinedResult = arrayfun(@or, combinedResult, replacedParts(:, colIndex));
            end

            % Sum up the final results
            result = sum(combinedResult);
        end
        function ret=checkRepresentation(obj)
            % 表現(準同型)の検証
            [rel,~,comm,inv]=obj.get2vRelation;
            % rel検証
            reps=arrayfun(@rep,rel);
            ret1=arrayfun(@(x)x==0,reps);
            % comm検証
            ret2=false(1,size(comm,2));
            if isnan(comm(1))
                ret2=true;
            else
                base=obj.base;
                for i=1:size(comm,2)
                    tmp1=obj.set_cp(1,{comm(1,i)});
                    tmp2=obj.set_cp(1,{comm(2,i)});
                    tmp1=rep(tmp1);
                    tmp2=rep(tmp2);
                    tmp=tmp1*tmp2-tmp2*tmp1;
                    ret2(i)=tmp==0;
                end
            end
            % inv検証
            ret3=false(1,size(inv,2));
            if isnan(inv(1))
                ret3=true;
            else
                base=obj.base;
                for i=1:size(inv,2)
                    tmp1=obj.set_cp(1,{inv(1,i)});
                    tmp2=obj.set_cp(1,{inv(2,i)});
                    tmp1=rep(tmp1);
                    tmp2=rep(tmp2);
                    tmp=tmp1*tmp2-1;
                    ret3(i)=tmp==0;
                end
            end
            
            ret=[ret1,ret2,ret3];
        end


        %% 計算基盤

        % 関係式の取得
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel=obj.empty;
                S.comm=[nan;nan];
                S.inv=[nan;nan];
                disp("set relations")
                S=obj.get2vRelation_(S);
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        % 関係式取得の補助メソッド
        function S=get2vRelation_(obj,S)
            S.mlist={};
            for i=1:numel(S.rel)
                S.rel(i)=combineTerm(S.rel(i));
                S.mlist(i)=S.rel(i).pw(end);
            end
            assert(size(S.comm,1)==2)
        end


        % 各項への作用
        function ret=lfun(obj,fun)
            % lfun funを線形作用させる
            ret=lfun_(obj,fun).calc();
        end
        function ret=lfun_(obj,fun)
            % lfun_ 簡約化処理無しの線形作用
            C=cell(obj.term,2);
            if obj.rank==0
                [C{1,:}]=fun({});
            else
                for i=1:obj.term
                    [C{i,:}]=fun(obj.pw(i,:));
                    % disp(obj.set_cp(C{i,:}))
                end
                % Ci=[mat2cell(obj.pw,ones(1,obj.term),obj.rank), ...
                %     mat2cell(obj.bs,ones(1,obj.term),obj.rank)];
                % [C(:,1),C(:,2),C(:,3)]=cellfun(fun,Ci(:,1),Ci(:,2),UniformOutput=false);
                    % disp(obj.set_cp(C{i,:}))
                
            end
            C(:,1)=arrayfun(@(x,y)x{1}*y,C(:,1),obj.cf,UniformOutput=false);
            cf=vertcat(C{:,1});
            pw=vertcat(C{:,2});
            % bs=vertcat(C{:,3});
            ret=obj.set_cp(cf,pw);
        end
        function ret=algID(obj)
            ret=@(p,b)obj.make(1,{p},b);
        end
        function ret=algfun(obj,funs,units)
            % algfun 代数準同型の作用
            % funs,unitsをテンソル階数の分だけ繰り返し入力する
            % funs:具体的にはStrAlg().algIDで返される関数形
            % funs:(power,base)→stralg
            arguments
                obj
            end
            arguments(Repeating)
                funs
                units
            end
            % lfun_ 簡約化処理無しの線形作用
            Rank=obj.rank;
            assert(length(funs)==Rank,'invalid rank of the algebra morphism domain')
            ret=obj.lfun_(@linfun);
            unit=horzcat(units{:});
            ret.base=[unit.base];
            % ret.ZERO=cellfun(@(x)x.set_cp(0,cell(1,x.rank)),units);
            function [c,p]=linfun(p)
                factors=units;
                X=StrAlg().scalar;
                for i=1:Rank
                    for j=1:length(p{i})
                        factors{i}=factors{i}*funs{i}(p{i}(j));
                    end
                    X=X|factors{i};
                end
                c=X.cf;
                p=X.pw;
                % b=X.bs;
            end
        end

        %　簡約化=関係式適用＋同次項括り＋零係数項削除
        function arg=calc(arg)
            arg=replace(arg,30);
            arg=combineTerm(arg);
            arg=removeZero(arg);
        end
        function obj=replace(obj,Ntimes)
            obj.sortedFlag=false(obj.term,obj.rank);
            obj.inverseSimplifiedFlag=false(obj.term,obj.rank);
            for i=1:Ntimes
                for j=1:obj.rank
                    obj=obj.replace2v_(j);
                end
                if all(obj.sortedFlag,"all"), break; end
            end
            % 逆元による相殺
            [~,~,~,inv]=obj.get2vRelation;

            Pinv0=[inv(1,:) inv(2,:)];
            Pinv1=[inv(2,:) inv(1,:)];
            for k=1:numel(obj.pw)
                if ~obj.sortedFlag(k)||isempty(obj.pw{k}), continue; end
                P=obj.pw{k};
                Ccnt0=sum(Pinv0==P',1);
                Ccnt1=sum(Pinv1==P',1);
                Ccnt=min([Ccnt0;Ccnt1]);
                deleteidx=nan(1,sum(Ccnt));
                idx=0;
                for i = 1:numel(Pinv0)
                    if Ccnt(i)==0,continue; end
                    deleteidx(idx+(1:Ccnt(i)))=find(P == Pinv0(i),Ccnt(i));
                    idx=idx+Ccnt(i);
                end
                P(deleteidx)=[];
                obj.pw{k}=P;
                % obj.bs{k}(deleteidx)=[];

            end
        end
        function obj=setSortedFlag(obj,flag)
            obj.sortedFlag=flag;
        end
        function [ret]=replace2v_(obj,Fidx)
            % 2変数の関係式の適用による簡約化をする
            sortedFlag=obj.sortedFlag;
            Z=obj.ZERO(Fidx);
            % 関係式の取得
            [rel,mlist,comm,inv]=Z.get2vRelation;
            R=obj.rank;
            cnt=0;
            ret=obj.lfun_(@fun);
            ret.sortedFlag=sortedFlag;
            if CR.H.replacingDisplay
                disp(obj)
            end
            function [c,p]=fun(p,b)
                % 各項のソートを行う線形関数
                P=p{1,Fidx};
                cnt=cnt+1;

                % if isequal(p{1},[3 1])
                %     disp("stop")
                % end

                % ソート済みをスキップ
                if sortedFlag(cnt,Fidx)||length(P)<=1
                    c=1;
                    sortedFlag(cnt,Fidx)=true;
                    return
                else
                    % いる？
                    sortedFlag(cnt,Fidx)=false;
                end

                % 各因子の反転をおこなう
                for i=1:length(P)-1
                    % flipFlag=false(1,length(mlist));
                    % 交換する組で、かつ順番が入れ替わっているもの
                    if P(i)>P(i+1)&&~isempty(comm)&&any(all(comm==P([i i+1])'))
                        P([i i+1])=P([i+1 i]);
                        c=1;
                        p{Fidx}=P;
                        % b{Fidx}([i i+1])=b{Fidx}([i+1 i]);
                        return
                    end
                    % 関係式の適用を判定
                    len0=length(P)-i+1;
                    lens=cellfun(@length,mlist);
                    for j=1:length(mlist)
                        % flipFlag(j)=isequal(mlist{j},P([i i+1]));
                        if ~(lens(j)<=len0&&isequal(mlist{j},P(i:i+lens(j)-1)))
                            continue
                        end
                        c=-rel(j).cf(1:end-1)/rel(j).cf(end);
                        N=length(c);
                        p=repmat(p,N,1);
                        % b=repmat(b,N,1);
                        p(:,Fidx)=cellfun(@(x)[P(1:i-1) x P(i+lens(j):end)], ...
                            rel(j).pw(1:end-1),UniformOutput=false);
                        % b(:,Fidx)=cellfun(@(x)[b{1}(1:i-1) x b{1}(i+lens(j):end)], ...
                        %     rel(j).bs(1:end-1),UniformOutput=false);
                        sortedFlag=[sortedFlag([1:cnt-1 cnt*ones(1,length(c)) cnt+1:end],:)];
                        sortedFlag(cnt+(1:length(c))-1,Fidx)=false;
                        cnt=cnt+length(c)-1;
                        return
                        % 関係式確認用
                        rel(j)
                    end
                end
                c=1;
                sortedFlag(cnt)=true;
            end
        end

        % 零係数項削除,係数簡約化ステップ
        function obj=removeZero(obj)
            if isequal(obj.ctype,NumericType.S)
                idx=~isAlways(obj.cf==0,Unknown="false");
                % idx=abs(subs(i1.cf,retq(),0.71))>0.000001;
                % idx=1:length(i1.cf);
                % idx=find(i1.cf);
                % s=@(x)x;
                % s=@simplify
                obj=obj.set_cp(simplify(sym(obj.cf(idx))),obj.pw(idx,:));
            else
                idx=abs(obj.cf)>1000*eps(obj.cf);
                obj=obj.set_cp(obj.cf(idx),obj.pw(idx,:));
            end
        end
        function verify(obj)
            if isempty(obj)
                return
            elseif ~isscalar(obj)
                arrayfun(@verify,obj)
                return
            end
            N=obj.term;
            R=size(obj.pw,2);
            L=cellfun(@(x)size(x,2),obj.pw);
            assert(isa(obj.cf,obj.ctype.class), ...
                '係数体エラー\n expected:%s, actual:%s', ...
                class(obj.cf),obj.ctype.class)
            assert(isequal(size(obj.pw),[N R]))
            for i=1:N
                for j=1:R
                    % assert(length(obj.bs{i,j})==L(i,j),'bsサイズエラー')
                    for k=1:L(i,j)
                        obj.verifyBase(obj.pw{i,j}(k))
                    end
                end
            end
        end
        function verifyBase(obj,pw)
            % assert(pw<=bs.dim,'pw範囲外')
        end
    end

    %% objの変更,生成
    methods (Static, Sealed, Access = protected)
        function obj=convertObject(target,arg)
            obj=eval(target);
            obj=obj.set_cp(arg,{[]});
            obj.base=Bases(0,strings(0),"empty");
        end
    end
    methods
        function obj=make(obj,cf,pw,bs)
            arguments
                obj
                cf
                pw (:,1) cell {verifyPW}
                bs (1,1) Bases 
            end
            bs.ZERO=[obj.set_cp(0,{[]})];
            obj.base=bs;
            % degs=cellfun(@length,pw);
            % bsCell=arrayfun(@(n)repmat(bs,1,n),degs,UniformOutput=false);
            obj=obj.set_cp(cf,pw);
        end

        % function prod()
        % コンストラクタ
        function obj=StrAlg(X)
            if nargin==1&&isa(X,"StrAlg")
                obj=obj.set_cp(X.cf,X.pw);
                obj.ZERO=X.ZERO;
                obj.base=X.base;
                return
            end
                % disp(class(obj))
        end
        

        function obj=set_cp(obj,cf,pw)
            try
                if isempty(cf)
                    obj.cf=obj.ctype.zero;
                    obj.pw=cell(1,obj.rank);
                    return
                end
                obj.cf=cf;
                obj.pw=pw;
            catch ME
                mustBeA(pw,"cell")
                rethrow(ME)
            end
        end



        function ret=testFunc(i1,i2,i3)
            ret=i1.lfun_(@fun);
            function [c,p]=fun(p)
                c=i2.cf;
                p=cellfun(@(x)[p{1} x],i2.pw,UniformOutput=false);
            end
        end

        %% 表示
        function ret=pol(arg)
            if arg.dim~=0
                vars=sym(arg.base.string);
                G=prod(vars.^arg.pw,2);
                ret=arg.cf.*G;
            else
                ret=sym(arg.cf);
            end
        end

        function disp(i1)
            if isempty(i1)
                disp("Empty StrAlg: ("+join(string(size(i1)),",")+")")
            else
                for ii=1:numel(i1)
                    disp_(i1(ii))
                end
            end
        end
        function disp_(i1)
            try
            feval("disp"+AlgebraConfig.H.disp_StrAlg,i1);
            catch ME
                warning(ME.identifier,'cannot display properly:\n %s',ME.message)
                i1.disp0
            end
        end
        function disp0(arg)
            builtin("disp",arg)
        end

        function ret=convertBaseString(obj,pw,arg)
            bsstrArr=string(arg);
            ret=bsstrArr(pw);
            if isempty(ret)
                ret="1";
            end
        end
        % テーブル形式表示
        function disp1(arg)
            % Example string array with duplicates
            bsstr=strings(arg.term,arg.rank);
            for i=1:arg.rank
                bc=arg.base(i);
                for j=1:arg.term
                    bsstr(j,i)=join(arg.convertBaseString(arg.pw{j,i},bc)," ");
                end
            end
            % bsstr=cellfun(@arg.convertBaseString,arg.pw,arg.bs,UniformOutput=false);
            if isempty(bsstr)
                disp(table(arg.cf,VariableNames="coeff"))
                return
            end
            base=categorical(join(bsstr," ⊗ ",2));
            % base(ismissing(base))=categorical("1");
            coeff=arg.cf;


            try
                disp(table(coeff,base));
            catch
                arg.verify
                arg.disp0
            end
        end
        % 数式の形式の表示
        function disp2(arg)
            disp(join(arg.convertTermToSym," + "))
        end
        function ret=string(arg)
            ret=join(arg.convertTermToSym," + ");
            % c_s=string(obj.c);
            % vn=obj.s.vname;
            % baseStr=cellfun(@(b,p){"*"+join(vn(b)+"^"+p,'*')},obj.b,obj.p);
            % emptyBase=cellfun(@isempty,obj.b);
            % baseStr(emptyBase)={""};
            % ret=join("("+c_s(:)+")"+[baseStr{:}]',"+");
        end
        function ret=latex(arg)
            str=join(arg.convertTermToSym," + ");
            str=strrep(str,"|","\otimes ");
            str=strrep(str,"*","");
            ret=str;
        end
        % 多数項の表示
        function disp3(i1)
            assert(numel(i1)==1)
            disp()
        end
        function ret=convertTermToSym(arg)
            % 数式表示の
            bs=num2cell(repmat(arg.base,size(arg.pw,1),1));
            bsstr=cellfun(@arg.convertBaseString,arg.pw,bs,UniformOutput=false);
            if isempty(bsstr)
                disp(table(arg.cf,VariableNames="coeff"))
                ret=string(arg.cf);
                return
            end
            base=join(cellfun(@(bsstr)join(bsstr,"*"),bsstr),"|",2);
            % Term="zzzz"+(1:arg.term)';
            % format=sum(sym(arg.cf).*(1+sym(Term)));
            % formatstr=string(format);
            % for i=1:arg.term
            %     formatstr=strrep(formatstr,"(zzzz"+i+" + 1)",base(i));
            % end
            % ret=formatstr;
            % return
            % formatの変数zzzzを置き換える処理がうまくいかないので下記を実行する
            % if isnumeric(arg.cf)
            %     coeff=arrayfun(@(x)sprintf("%+g",x),arg.cf);
            % else
            %     coeff=string(arg.cf);
            % end
            coeff="("+string(arg.cf)+")";
            if arg.rank>1
                base="("+base+")";
            end
            ret=coeff+"*"+base;

        end
        % convert to sym
        function ret=sym(obj)
            if true
                F=@simplify;
            else
                F=@(x)x;
            end
            ret=arrayfun(@(x)F(sum(x.pol)),obj);
        end


        %% additional function
        function validate(obj)
            assert(all(size(obj.pw)==[length(obj.cf),obj.rank]))
        end
        function ret=ones(obj)

        end
        function z = zeros(obj,varargin)
            if nargin==1
                z=obj.set_cp([]);
            elseif any([varargin{:}] <= 0)
                z = PolAlg.empty(varargin{:});
            else
                z = repmat(PolAlg,varargin{:});
            end
        end


        %複製
        function ret=matrix(obj,i1)
            ret=repmat(obj,i1);
        end
        % function obj=setBase(obj,algB,Zero)
        %     obj.base=algB;
        %     obj.ZERO=Zero;
        % end
        % function ret=get.base(obj)
        %     if isempty(obj.ZERO)
        %         ret=obj.base;
        %     else
        %         ret=[obj.ZERO.base];
        %     end
        % end
        function ret=get.ZERO(obj)
            ret=[obj.base.ZERO];
        end
        function ret=subs(obj,varargin)
            ret=obj;
            ret.cf=subs(obj.cf,varargin{:});
            ret.pw=subs(obj.pw,varargin{:});
        end
        function ret=get.term(obj)
            ret=length(obj.cf);
        end
        function ret=get.deg(obj)
            ret=max(cellfun(@length,obj.pw));
        end
        function ret=get.dims(obj)
            ret=[obj.base.dim];
        end
        function ret=get.rank(obj)
            ret=length(obj.base);
        end
        function ret=scalar(obj)
            ret=obj.set_cp(1,{});
            % ret.ZERO=StrAlg.empty;
            ret.base=Bases.empty;
        end
    end
    methods(Hidden)
        function ret=dimP(obj)
            ret=size(obj.pw,2);
        end
    end
    %zeros
    methods (Static)

        function F=tensorMor(func,dim)
            arguments(Repeating)
                func
                dim
            end
            N=numel(func);
            F=@fun;
            dims=cell2mat(dim);
            idx=arrayfun(@(x,y){x-(0:y-1)},cumsum(dims),dims);
            function [c,p]=fun(p)
                C=cell(1,N);
                P=cell(1,N);
                for i=1:N
                    [C{i},P{i}]=func{i}(p(idx{i}));
                end
                terms=cellfun(@numel,C);
                idx0=terms==0;
                C{idx0}=0;
                P{idx0}=zeros(1,0);
                terms(idx0)=1;
                combi=arrayfun(@(x){1:x},terms);
                T=combinations(combi{:});
                h=height(T);
                w=width(T);
                % cの初期化？
                TC=cell(h,w);
                TP=cell(h,w);
                for i=1:w
                    TC{:,i}=C{i}(T{:,i});
                    TP{:,i}=P{i}(T{:,i},:);
                end
                c=zeros(h,1,CR.H.cft);
                p=zeros(h,1,CR.H.pft);
                for i=1:h
                    for j=1:w
                        c(i)=c(i)+TC{i,j};
                    end
                    p(i,:)=hotzcat(TP{i,:});
                end
            end
        end
    end

end


%% Local method
% 同類項括り
function ret=combineTerm(arg)
    if isa(arg,'strEndV')&&false
        % strEndVが同類項を間違えてくくってしまうため，追加
        ret=arg;
        return
    end

    cf=arg.cf;
    pw=arg.pw;
    % bs=arg.bs;
    t=size(arg.pw,1);
    Rank=size(arg.pw,2);
    if isempty(pw)
        ret=arg.set_cp(sum(cf),pw);
        return
    end
    maxDegs = nan(1,Rank); % Maximum array length
    % Prepare padded power arrays for sorting
    paddedArrays=cell(size(arg.pw));
    for i=1:Rank
        maxDegs(i)=max([0;cellfun(@length,pw(:,i))]);
        paddedArrays(:,i) = cellfun(@(x) [-Inf(1, maxDegs(i) - length(x)),x], pw(:,i), 'UniformOutput', false);
    end
    sortingMatrix = cell2mat(paddedArrays);

    % Sort terms by their padded power arrays (length first, then dictionary order)
    [sortedPw,sortIdx]=sortrows(sortingMatrix);
    groupedTerms=cell(t,2);
    cf=cf(sortIdx);
    % bs=bs(sortIdx,:);
    % Group terms with identical power arrays
    groupIdx=1;
    currentPw=sortedPw(1,:);
    for ii=1:numel(sortIdx)
        if ~isequal(currentPw,sortedPw(ii,:))
            groupIdx=groupIdx+1;
            currentPw=sortedPw(ii,:);
        end
        groupedTerms{groupIdx,1}(end+1)=cf(ii);
        groupedTerms{groupIdx,2}=currentPw;
        % groupedTerms{groupIdx,3}(end+1,:)=bs(ii,:);
    end
    Nterm=groupIdx;
    groupedTerms=groupedTerms(1:Nterm,:);
    % Reduce grouped terms (combine coefficients and pick one base)
    reducedTerms=cell(Nterm,1);
    for i=1:Nterm
        reducedTerms{i}=reduceGroup(groupedTerms(i,:));
    end
    groupedTerms=vertcat(reducedTerms{:});
    % Extract reduced powers, remove padding, and combine all terms
    pw=cellfun(@(factors)mat2cell(factors,1,maxDegs),groupedTerms(:,2),UniformOutput=false);
    pw=vertcat(pw{:});
    % bs=cell(size(pw));
    for i=1:Nterm
        for j=1:Rank
            pw{i,j}(pw{i,j}==-Inf)=[];
        end
    end
    % bs=vertcat(groupedTerms{:,3});
    cf=vertcat(groupedTerms{:,1});
    ret=arg.set_cp(cf,pw);

    function termGroup=reduceGroup(termGroup)
        %  Helper function to combine terms in a group
        termGroup{1}=sum(termGroup{1});
        % termGroup{3}=termGroup{3}(1,:);
    end


end
function [q,r,flag]=helpdiv(q,r,d)
    arguments
        q PolAlg
        r PolAlg
        d PolAlg
    end
    flag=false;
    if length(r.cf)<length(d.cf)
        return
    elseif length(r.cf)==length(d.cf)
        flag=true;
    end
    Q=d;
    Q=d.set_cp(r.cf(1)/d.cf(1),r.pw(1,:)-d.pw(1,:));
    q=q+Q;
    r=r-Q*d;
    if flag&&(length(r.cf)==length(d.cf))
        disp("infinite loop")
        return
    end
    [q,r]=helpdiv(q,r,d);
end


function cp=set_cptype(cp,typ)
    % SET_CP c,p型設定用ヘルパー
    if isa(cp,typ), return; end
    try
        cp=feval(typ,cp);
    catch ME
        warning cannot_convert
        rethrow(ME)
    end
end

function verifyPW(arg)
    tf=cellfun(@(x)isnumeric(x)&&size(x,1)<=1,arg);
    assert(all(tf),"invalid power: %s",join(string(find(~tf))))
end
function verifyBS(arg)
    tf=cellfun(@(x)isa(x,'Bases')&&size(x,1)==1,arg);
    assert(all(tf),"invalid bases: %s",join(string(find(~tf))))
end
