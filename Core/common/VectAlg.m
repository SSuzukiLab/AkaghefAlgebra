classdef(InferiorClasses=?sym) VectAlg<IAdditive&matlab.mixin.indexing.RedefinesBrace
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        cf  %coefficient
        % 1,2,3の階数方向に係数を立体的に並べる．rank=1なら縦ベクトルなので注意

        bs (1,:) Bases %basis
        ZERO (1,:) cell %zero element of each vector space
        SC TypeParam=TypeParam([])%structure constants
    end
    properties(Dependent)
        dim %dimension of the vector space
        dims %dimension of each vector space
        rank % Tensor rank
    end
    %% binomial operation
    methods
        function obj=setBase(obj,base)
            % setBase 基底の設定
            obj.bs=base;
            obj.cf=Czeros(obj,base.dim,1);
            obj.ZERO={obj};
        end
        function ret=getSC(obj,arg)
            ret=obj.SC.get([obj.bs.name arg]);
        end

        function [i1,i2]=alignNum(i1,i2)
            % 型をstrAlgにする
            if ~isequal(class(i1),class(i2))||~isequal(i1.bs,i2.bs)
                tf=[isa(i1,"VectAlg") isa(i2,"VectAlg")];
                if isequal(tf,[1 0])
                    i2=castCtype(i1,i2);
                elseif isequal(tf,[0 1])
                    i1=castCtype(i2,i1);
                elseif isequal(tf,[1 1])
                    try
                        i2=i1.casttype(i2);
                    catch
                        try
                            i1=i2.casttype(i1);
                        catch
                            error("symp:alignNum","invalid input")
                        end
                    end

                end
            end
        end
        function ret=casttype(obj,arg)
            % casttype 型を合わせる
            if arg.rank==0
                ret=obj.unit;
                ret.cf=arg*ret.cf;
            else
                error not_impl
            end
        end
        function ret=castCtype(obj,arg)
            % casttype 型を合わせる
            ret=obj.unit;
            ret.cf=arg*ret.cf;
        end
        function ret=identifier(obj)
            ret=obj.bs.name;
        end


        function ret=plus(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            assert(i1.dim==i2.dim,'different dim addition error')
            try
                sz=size(zeros(size(i1))+zeros(size(i2)));
            catch
                error("symp:plus","size dimensions must match")
            end
            if any(sz==0)
                ret=repmat(i1,sz);
            else
                % assert(configure(i1,i2))
                ret=arrayfun(@plus_,repmat(i1,sz./size(i1)),repmat(i2,sz./size(i2)));
            end
            % end

            function i1=plus_(i1,i2)
                % PLUS_ supplementary method for scalar
                i1.cf=i1.cf+i2.cf;
            end
        end
        % function ret=minus(i1,i2)
        %     ret=i1+(-i2);
        % end
        function i1=uminus(i1)
            i1.cf=-i1.cf;
        end

        function ret=eq(i1,i2)
            sub=calc(i1-i2);
            if sub.bs.getCtype=="S"
                ret=fold(@and,all(sub.cf==0),symtrue);
            else
                ret=all(sub.cf==0);
            end
        end
        %乗算,作用
        function ret=mtimes(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            assert(isequal(i1.bs,i2.bs),'異なる空間での積エラー')
            ret=i1;
            R=i1.rank;
            cf=tensorprod(i1.cf,i2.cf,Num=R);
            for k=1:R
                M=i1.ZERO{k}.getSC('_μ');
                cf=tensorprod(cf,M,[1 R+2-k],[1 2]);
            end
            ret.cf=cf;
        end
        function ret=times(i1,i2)

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
            ret=i1.unit;
            for ii=1:i2
                ret=i1*ret;
            end
        end
        function ret=unit(arg)
            eta=1;
            for k=1:arg.rank
                eta_=arg.ZERO{k}.getSC('_η');
                eta=tensorprod(eta,eta_,Num=k);
            end
            ret=arg.set_c(shiftdim(eta,1));
        end
        %テンソル積
        function o=or(i1,i2)
            try
                [i1,i2]=alignNum(i1,i2);
            catch ME
                if ~(isa(i1,"VectAlg")&&isa(i2,"VectAlg"))
                    % (g|g)|gのような元に対処する
                    rethrow(ME)
                end
            end
            o=i1;
            o.bs=[i1.bs i2.bs];
            o.ZERO=[i1.ZERO i2.ZERO];
            o.bs=[i1.bs i2.bs];
            o.cf=tensorprod(i1.cf,i2.cf,Num=i1.rank);
        end
        %副積
        function o=and(i1,i2)
            error('not implemented')
        end
        %副積計算
        function o=not(i1)
            error('not implemented')
        end

        function ret=Delta(obj)
            ret=obj|obj;
            ret.cf(:)=0;
            C=obj.SC.get([obj.identifier '_Δ']);
            ret.cf=tensorprod(C,obj.cf,1,1);
            % for k1=1:obj.dim
            %     for k2=1:obj.dim
            %         for k3=1:obj.dim
            %             ret.cf(k2,k3)=ret.cf(k2,k3)+C(k1,k2,k3)*obj.cf(k1);
            %         end
            %     end
            % end
        end
        function ret=counit(obj)
            ep=obj.SC.get([obj.identifier '_ε']);
            ret=obj.cf.'*ep;
            ep=1;
            for k=1:arg.rank
                eta_=arg.ZERO{k}.getSC('_ε');
                eta=tensorprod(eta,eta_,Num=k);
            end
            ret=arg.set_c(shiftdim(eta,1));

        end
        function ret=S(obj)
            % S: Hopf algebra antipode
            S=obj.SC.get([obj.identifier '_S']);
            ret=obj;
            ret.cf(:)=0;
            for k1=1:obj.dim
                for k2=1:obj.dim
                    ret.cf(k2)=ret.cf(k2)+S(k1,k2)*obj.cf(k1);
                end
            end
        end
        %% 計算基盤

        % 各項への作用
        function ret=lfun(obj,fun)
            % lfun funを線形作用させる
            ret=lfun_(obj,fun).calc();
        end
        function ret=lfun_(obj,fun)
            % lfun_ 簡約化処理無しの線形作用
            C=cell(obj.term,3);
            if obj.rank==0
                [C{1,:}]=fun({},{});
            else
                for i=1:obj.term
                    [C{i,:}]=fun(obj.pw(i,:),obj.bs(i,:));
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
            bs=vertcat(C{:,3});
            ret=obj.set_cp(cf,pw,bs);
        end
        function ret=algID(obj)
        end
        function ret=algfun(obj,funs,units)
            % algfun 代数準同型の作用
            % funs,unitsをテンソル階数の分だけ繰り返し入力する
            % funs:具体的にはstrAlg().algIDで返される関数形
            % funs:(power,base)→stralg
        end

        %　簡約化=関係式適用＋同次項括り＋零係数項削除
        function arg=calc(arg)
            % arg=replace(arg,30);
            % arg=combineTerm(arg);
            % arg=removeZero(arg);
        end
    end
    %% objの変更,生成
    methods
        function obj=make(obj,cf,idx)
            arguments
                obj
                cf
            end
            arguments(Repeating)
                idx
            end
            obj.cf(:)=0;
            idx=sub2ind([obj.dims 1],idx{:});
            obj.cf(idx)=cf;
        end

        % function prod()
        % コンストラクタ
        function obj=VectAlg(X)
            if nargin==1
                error("vect:cast","implicit casting not allowed:%s",string(X))
            end
        end
        function obj=setSC(obj,identifier,mu,eta,Delta,ep,S)
            % setSC ホップ代数の構造定数を設定する
            obj.SC.insert([identifier '_μ'],mu);
            obj.SC.insert([identifier '_η'],eta);
            obj.SC.insert([identifier '_Δ'],Delta);
            obj.SC.insert([identifier '_ε'],ep);
            obj.SC.insert([identifier '_S'],S);
            Si=S^-1;
            obj.SC.insert([identifier '_Si'],Si);
            P=tensorprod(tensorprod(tensorprod(mu,S,3,1),Delta,3,1),S,[3 1],[1 2]);
            [U,K,V] = svd(P);
            Ir=U(:,1);
            Cr=K(1,1)*V(:,1);
            Cr=Cr/(Ir.'*Cr);
            Cl=Si*Cr;
            Il=S.'*Ir;
            obj.setIntegrals(Ir,Cr,Il,Cl);
        end
        function setIntegrals(obj,Ir,Cr,Il,Cl)
            % Ir,Cr,Il,Cl (D,1) size
            obj.SC.insert([obj.identifier '_intr'],Ir);
            obj.SC.insert([obj.identifier '_cointr'],Cr);
            obj.SC.insert([obj.identifier '_intl'],Il);
            obj.SC.insert([obj.identifier '_cointl'],Cl);
            obj.SC.insert([obj.identifier '_CrIl'],dot(Cr,Il));         
        end
        function ret=verify(obj)
            try
                assert(isequal(size(obj.cf),obj.dims)||(obj.rank==1&&numel(obj.cf)==obj.dim))
            catch ME
                obj.disp0
                rethrow(ME)
            end
        end
        function ret=verifyHopf(obj,isequal_)
            % verifyHopf ホップ代数の構造定数を確認する
            D=obj.dim;
            I=eye(D);
            mu=obj.getSC(['_μ']);
            eta=obj.getSC(['_η']);
            Delta=obj.getSC(['_Δ']);
            ep=obj.getSC(['_ε']);
            S=obj.getSC(['_S']);
            Si=obj.getSC(['_Si']);

            % size
            assert(isequal(size(mu),[D D D]),"invalid size:μ")
            assert(isequal(size(eta),[D 1]),"invalid size:η")
            assert(isequal(size(Delta),[D D D]),"invalid size:Δ")
            assert(isequal(size(ep),[D 1]),"invalid size:ε")
            assert(isequal(size(S),[D D]),"invalid size:S")
            assert(isequal(size(Si),[D D]),"invalid size:Si")

            % unitality
            tmp=tensorprod(eta,mu,1,1,Num=1);
            assert(isequal_(tmp,I),"unitality error")
            tmp=tensorprod(eta,mu,1,2,Num=1);
            assert(isequal_(tmp,I),"unitality error")
            % counitality
            tmp=tensorprod(ep,Delta,1,2,Num=1);
            assert(isequal_(tmp,I),"counitality error")
            tmp=tensorprod(ep,Delta,1,3,Num=1);
            assert(isequal_(tmp,I),"counitality error")

            % associativity
            tmp=tensorprod(mu,mu,3,1);
            tmp2=permute(tensorprod(mu,mu,2,3),[1 3 4 2]);
            assert(isequal_(tmp,tmp2),"associativity error")
            % coassociativity
            tmp=permute(tensorprod(Delta,Delta,2,1),[1 3 4 2]);
            tmp2=tensorprod(Delta,Delta,3,1);
            assert(isequal_(tmp,tmp2),"coassociativity error")

            % bialgebra property
            %  Δ∘μ = (μ⊗μ)∘τ23∘(Δ⊗Δ)
            tmp = tensorprod(Delta,Delta,Num=3);
            tmp2 = tensorprod(mu,mu,Num=3);
            tmp = tensorprod(tmp,tmp2,[2 5 3 6],[1 2 4 5]);
            tmp2= tensorprod(mu,Delta,3,1);
            assert(isequal_(tmp,tmp2),"bialgebra property error")


            % antipode property
            % μ ∘ (S ⊗ id) ∘ Δ = η ∘ ε  and  μ ∘ (id ⊗ S) ∘ Δ = η ∘ ε
            tmp2=tensorprod(ep,eta,2,2,Num=2);

            tmp=tensorprod(S,mu,2,1);
            tmp=tensorprod(Delta,tmp,[2 3],[1 2]);
            assert(isequal_(tmp,tmp2),"antipode left inverse error")

            tmp=tensorprod(Delta,S,3,1);
            tmp=tensorprod(tmp,mu,[2 3],[1 2]);
            assert(isequal_(tmp,tmp2),"antipode right inverse error")

            % integral property
            %  δ_r*u=ε(u)δ_r, u*δ_l=ε(u)δ_l,
            % ∫_r(u_1)u_2=∫_r(u), u_1∫_l(u_2)=∫_l(u)

            intr=obj.getSC(['_intr']);
            cointr=obj.getSC(['_cointr']);
            intl=obj.getSC(['_intl']);
            cointl=obj.getSC(['_cointl']);
            tmp=tensorprod(cointr,mu,1,1,Num=1);
            tmp2=tensorprod(ep,cointr,Num=1);
            assert(isequal_(tmp,tmp2),"right cointegral property error")
            tmp=tensorprod(intr,Delta,1,2,Num=1);
            tmp2=tensorprod(intr,eta,Num=1);
            assert(isequal_(tmp,tmp2),"right integral property error")
            % tmp=tensorprod(intl,cointr,1,2,Num=1);
            % tmp2=tensorprod(intl,cointr,Num=1);
            % leftはまだ未実装
            tmp=[intr.'*cointr, intl.'*cointr; ...
                 intr.'*cointl, intl.'*cointl];
            assert(isequal_(tmp([1 2 4]),[1 1 1]),"integral normalization error")

            disp("Confirmed to be a Hopf algebra")
        end
        function verifyInt(obj)
            isequal_=@(x,y)all(eqD(x,y,1e-5),"all");
            D=obj.dim;
            I=eye(D);
            mu=obj.getSC(['_μ']);
            eta=obj.getSC(['_η']);
            Delta=obj.getSC(['_Δ']);
            ep=obj.getSC(['_ε']);
            S=obj.getSC(['_S']);
            Si=obj.getSC(['_Si']);

            intr=obj.getSC(['_intr']);
            cointr=obj.getSC(['_cointr']);
            intl=obj.getSC(['_intl']);
            cointl=obj.getSC(['_cointl']);
            tmp=tensorprod(cointr,mu,1,1,Num=1);
            tmp2=tensorprod(ep,cointr,Num=1);
            assert(isequal_(tmp,tmp2),"right cointegral property error")
            tmp=tensorprod(intr,Delta,1,2,Num=1);
            tmp2=tensorprod(intr,eta,Num=1);
            assert(isequal_(tmp,tmp2),"right integral property error")
            % tmp=tensorprod(intl,cointr,1,2,Num=1);
            % tmp2=tensorprod(intl,cointr,Num=1);
            % leftはまだ未実装
            tmp=[intr.'*cointr, intl.'*cointr; ...
                 intr.'*cointl, intl.'*cointl];
            assert(isequal_(tmp([1 2 4]),[1 1 1]),"integral normalization error")
            
        end
        function dispInt(obj)
            intr=obj.getSC(['_intr']);
            cointr=obj.getSC(['_cointr']);
            intl=obj.getSC(['_intl']);
            cointl=obj.getSC(['_cointl']);
            T=table([intr cointr intl cointl].',VariableNames="coordinates",RowNames=["intr" "cointr" "intl" "cointl"]);
            disp(T)
        end

        function obj=set_c(obj,cf)
            if obj.rank==1
                cf=reshape(cf,[],1);
            end
            if isequal(size(obj.cf),size(cf))
                obj.cf=cf;
            else
                error("VectAlg:set_c","size mismatch")
            end
        end



        function ret=testFunc(i1,i2,i3)
            ret=i1.lfun_(@fun);
            function [c,p,b]=fun(p,b)
                c=i2.cf;
                p=cellfun(@(x)[p{1} x],i2.pw,UniformOutput=false);
                b=cellfun(@(x)[b{1} x],i2.bs,UniformOutput= false);
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
                disp("Empty strAlg: ("+join(string(size(i1)),",")+")")
            else
                for ii=1:numel(i1)
                    disp_(i1(ii))
                end
            end
        end
        function disp_(i1)
            feval("disp"+CR.H.displayRule,i1);
        end
        function disp0(arg)
            builtin("disp",arg)
        end

        % テーブル形式表示
        function disp1(arg)

            bsnames=fliplr(arrayfun(@string,arg.bs,UniformOutput=false));
            if isempty(bsnames)
                disp(table(arg.cf,VariableNames="coeff"))
                return
            end
            T=combinations(bsnames{:});
            base=categorical(join(fliplr(T{:,:})," ⊗ ",2));
            coeff=arg.cf(:);
            if CR.H.vectD1removeZero
                idx=arg.removeZero;
                coeff=coeff(idx);
                base=base(idx);
                if ~any(idx)
                    coeff=0;
                    base=categorical("-");
                end
            end
            disp(table(coeff,base));
        end
        function ret=removeZero(arg)
            ret=abs(arg.cf)>1e-10;
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
            bsstr=string(arg.bs);
            if isempty(bsstr)
                disp(table(arg.cf,VariableNames="coeff"))
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
        %% bracing
    end

    methods (Access=protected)
        function varargout = braceReference(obj,indexOp)
            idx=[0,0];
            [idx(1),idx(2)]=ind2sub(obj.dims,indexOp.Indices{1});

            Z=obj.ZERO{indexOp.Indices{2}};
            Z.ZERO={Z};
            Z.cf(idx(indexOp.Indices{2}))=obj.cf(indexOp.Indices{1});
            varargout={Z};
            % [varargout{1:nargout}] = obj.cf.(indexOp);
        end

        function obj = braceAssign(obj,indexOp,varargin)
            error non_impl
        end

        function n = braceListLength(obj,indexOp,indexContext)
            n = 1;
        end
    end
    methods

        %% additional function
        function ret=ones(obj)

        end
        function z = zeros(obj,varargin)
            if nargin==1
                z=obj.set_cp([]);
            elseif any([varargin{:}] <= 0)
                z = symp.empty(varargin{:});
            else
                z = repmat(symp,varargin{:});
            end
        end
        function ret=Czeros(obj,varargin)
            ret=zeros(varargin{:},obj.bs.getCtype.type);
        end


        %複製
        % function ret=matrix(obj,i1)
        %     ret=repmat(obj,i1);
        % end
        % function obj=setBase(obj,algB,Zero)
        %     obj.algbase=algB;
        %     obj.ZERO=Zero;
        % end
        % function ret=subs(obj,varargin)
        %     ret=obj;
        %     ret.cf=subs(obj.cf,varargin{:});
        %     ret.pw=subs(obj.pw,varargin{:});
        % end
        function ret=get.dims(obj)
            ret=obj.bs.dims;
        end
        function ret=get.dim(obj)
            ret=obj.bs.dim;
        end
        function ret=get.rank(obj)
            ret=length(obj.bs);
        end
        % function ret=scalar(obj)
        %     ret=obj.set_cp(1,{},{});
        %     ret.ZERO={};
        %     ret.algbase=Bases.empty;
        % end
    end
    methods(Hidden)

    end
    %zeros
    methods (Static)
        function F=tensorMor(func,dim)
        end
    end

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
