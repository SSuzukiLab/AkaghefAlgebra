classdef(InferiorClasses=?sym) VectAlg<IAdditive&matlab.mixin.indexing.RedefinesBrace
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        cf  %coefficient
        % 1,2,3の階数方向に係数を立体的に並べる．rank=1なら縦ベクトルなので注意

        bs (1,:) Bases %basis
        % ZERO (1,:) cell %zero element of each vector space
        ZERO 
        spec (1,1) SpaceSpec % specify principle of space
        % also, store structure constant 
        sparse SparseEx=SparseEx
    end
    properties(Dependent)
        dim %dimension of the vector space
        dims %dimension of each vector space
        rank % Tensor rank
    end
    %% binomial operation
    methods
        % attempt: switch storage to SparseEx
        function obj=set.cf(obj,value)
            obj.sparse=SparseEx(value);
        end
        function value=get.cf(obj)
            value=obj.sparse.toMatrix();
        end
        function obj=set.ZERO(obj,value)
            obj.bs.ZERO=value;
        end
        function value=get.ZERO(obj)
            value=[obj.bs.ZERO];
        end
        
        function obj=setBase(obj,base)
            % setBase 基底の設定
            obj.bs=base;
            base.ZERO=obj.zeros;
            obj.spec.base=base;
            % obj.spec=SpaceSpec(base); %issue: space specをここで決めるか？
            obj.cf=Czeros(obj,base.dim,1);
            
            % obj.ZERO={obj};
        end
        function varargout=getSC(obj,arg)
            varargout=cell(1,length(string(arg)));
            try
            [varargout{:}]=obj.spec.SC{arg};
            catch ME
                
                if strcmp(ME.identifier, 'MATLAB:dictionary:ScalarKeyNotFound')
                    error('The specified key "%s" was not found in the dictionary.',join(arg));
                else
                    rethrow(ME);
                end
            end
        end

        function [i1,i2]=alignNum(i1,i2)
            % 型をStrAlgにする
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
                            error("PolAlg:alignNum","invalid input")
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
                error("PolAlg:plus","size dimensions must match")
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
            Mu=struct;
            sp1=i1.sparse;
            sp2=i2.sparse;
            expr=sprintf("sp1{%s}sp2{%s}", ...
                join(string(1:R),","),join(string(R+1:2*R),","));
            for k=1:R
                M=i1.bs(k).ZERO.getSC('prod');
                Mu.("t"+k)=M;
                expr=expr+sprintf("Mu.t%d{%d,%d,%d}",k,k,k+R,k+2*R);
            end
            sp=calcTensorExpression(expr,(1:R)+2*R);
            ret.sparse=sp;
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
                eta_=arg.bs(k).ZERO.getSC('unit');
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
            % o.ZERO=[i1.ZERO i2.ZERO];
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
            C=obj.getSC(['coprod']);
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
            ep=obj.SC.get([obj.identifier 'counit']);
            ret=obj.cf.'*ep;
            ep=1;
            for k=1:arg.rank
                eta_=arg.bs(k).ZERO.getSC('counit');
                eta=tensorprod(eta,eta_,Num=k);
            end
            ret=arg.set_c(shiftdim(eta,1));

        end
        function ret=S(obj)
            % S: Hopf algebra antipode
            S=obj.spec.SC{'antipode'};
            ret=obj;
            ret.cf(:)=0;
            % for k1=1:obj.dim
            %     for k2=1:obj.dim
            %         ret.cf(k2)=ret.cf(k2)+S(k1,k2)*obj.cf(k1);
            %     end
            % end
            ret.cf=S*obj.cf; %250816 Sの定義のテンソルを転置した
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
            % funs:具体的にはStrAlg().algIDで返される関数形
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
            obj.spec=SpaceSpec;
            if nargin==1
                error("vect:cast","implicit casting not allowed:%s",string(X))
            end
        end
        function obj=setSC(obj,identifier,mu,eta,Delta,ep,S)
            % setSC ホップ代数の構造定数を設定する
            mu=SparseEx(mu);
            eta=SparseEx(eta);
            Delta=SparseEx(Delta);
            ep=SparseEx(ep);
            S=SparseEx(S);
            obj.spec.SC{'prod'}=mu;
            obj.spec.SC{'unit'}=eta;
            obj.spec.SC{'coprod'}=Delta;
            obj.spec.SC{'counit'}=ep;
            obj.spec.SC{'antipode'}=S;
            S_=S.toMatrix();
            Si=SparseEx(S_^-1);
            obj.spec.SC{'antipode_inv'}=Si;
            % Refer "Tensor_on_Vertex.jpg"
            % issue: U^*=AとA^*=Uどっちの演算？
            Tp=calcTensorExpression('mu{3,5,1}Delta{4,5,2}',[1,2,3,4]);
            Tm=calcTensorExpression('mu{3,6,1}S{6,5}Delta{4,5,2}',[1,2,3,4]);
            % Tp=TP(De,Mu,2,2); %250406
            % % Tm=PE(TP(TP(De,An,2,2),De,3,2),[2 1 4 3]);
            % Tm=eval(makeTensorExpression('De{1,4,0}An{4,5}Mu{3,5,2}',[0 1 2 3]));
            obj.spec.SC{'Tp'}=Tp;
            obj.spec.SC{'Tm'}=Tm;
            if ~isa(mu,'double'), return; end
            P=calcTensorExpression('mu{1,2,3}S{3,4}Delta{4,5,6}S{5,1}',[2,6]);
            % P=tensorprod(tensorprod(tensorprod(mu,S,3,1),Delta,3,1),S,[3 1],[1 2]);
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
            Ir=SparseEx(Ir);
            Cr=SparseEx(Cr);
            Il=SparseEx(Il);
            Cl=SparseEx(Cl);
            obj.spec.SC{'intr'}=Ir;
            obj.spec.SC{'cointr'}=Cr;
            obj.spec.SC{'intl'}=Il;
            obj.spec.SC{'cointl'}=Cl;
            obj.spec.SC{'CrIl'}=Cr.toMatrix().'*Il.toMatrix();
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
            I=SparseEx(eye(D));
            mu=obj.getSC('prod');
            eta=obj.getSC('unit');
            Delta=obj.getSC('coprod');
            ep=obj.getSC('counit');
            S=obj.getSC('antipode');
            Si=obj.getSC('antipode_inv');
            % size
            assert(isequal(mu.size,[D D D]),"invalid size:μ")
            assert(isequal(eta.size,[D]),"invalid size:η")
            assert(isequal(Delta.size,[D D D]),"invalid size:Δ")
            assert(isequal(ep.size,[D]),"invalid size:ε")
            assert(isequal(S.size,[D D]),"invalid size:S")
            assert(isequal(Si.size,[D D]),"invalid size:Si")

            % unitality
            % tmp=tensorprod(eta,mu,1,1,Num=1);
            tmp=calcTensorExpression('eta{1}mu{1,2,3}',[2,3]);
            assert(isequal_(tmp,I),"left unitality error")
            % tmp=tensorprod(eta,mu,1,2,Num=1);
            tmp=calcTensorExpression('eta{2}mu{1,2,3}',[1,3]);
            assert(isequal_(tmp,I),"right unitality error")
            % counitality
            % tmp=tensorprod(ep,Delta,1,2,Num=1);
            tmp=calcTensorExpression('ep{2}Delta{1,2,3}',[1,3]);
            assert(isequal_(tmp,I),"left counitality error")
            % tmp=tensorprod(ep,Delta,1,3,Num=1);
            tmp=calcTensorExpression('ep{3}Delta{1,2,3}',[1,2]);
            assert(isequal_(tmp,I),"right counitality error")

            % associativity
            % tmp=tensorprod(mu,mu,3,1);
            tmp=calcTensorExpression('mu{1,2,3}mu{3,4,5}',[1,2,4,5]);
            % tmp2=permute(tensorprod(mu,mu,2,3),[1 3 4 2]);
            tmp2=calcTensorExpression('mu{1,2,3}mu{4,5,2}',[1,4,5,3]);
            assert(isequal_(tmp,tmp2),"associativity error")
            % coassociativity
            % tmp=permute(tensorprod(Delta,Delta,2,1),[1 3 4 2]);
            tmp=calcTensorExpression('Delta{1,2,3}Delta{2,4,5}',[1,4,5,3]);
            % tmp2=tensorprod(Delta,Delta,3,1);
            tmp2=calcTensorExpression('Delta{1,2,3}Delta{3,4,5}',[1,2,4,5]);
            assert(isequal_(tmp,tmp2),"coassociativity error")

            % bialgebra property
            %  Δ∘μ = (μ⊗μ)∘τ23∘(Δ⊗Δ)
            % tmp = tensorprod(Delta,Delta,Num=3);
            % tmp2 = tensorprod(mu,mu,Num=3);
            % tmp = tensorprod(tmp,tmp2,[2 5 3 6],[1 2 4 5]);
            % tmp2= tensorprod(mu,Delta,3,1);
            % assert(isequal_(tmp,tmp2),"bialgebra property error")
            tmp=calcTensorExpression( ...
                'Delta{1,3,4}Delta{2,5,6}mu{3,5,7}mu{4,6,8}',[1,2,7,8]);
            tmp2=calcTensorExpression( ...
                'mu{1,2,3}Delta{3,7,8}',[1,2,7,8]);
            dif=tmp-tmp2;
            assert(isequal_(dif.val,0),"bialgebra property error")


            % antipode property
            % μ ∘ (S ⊗ id) ∘ Δ = η ∘ ε  and  μ ∘ (id ⊗ S) ∘ Δ = η ∘ ε
            % tmp2=tensorprod(ep,eta,2,2,Num=2);
            tmp2=calcTensorExpression('ep{1}eta{2}',[1,2]);

            % tmp=tensorprod(S,Delta,2,2);
            % tmp=tensorprod(tmp,mu,[1,3],[1,2]);
            tmp=calcTensorExpression('S{2,3}Delta{1,3,4}mu{2,4,5}',[1,5]);
            assert(isequal_(tmp,tmp2),"antipode left inverse error")

            % tmp=tensorprod(Delta,S,3,2);
            % tmp=tensorprod(tmp,mu,[2 3],[1 2]);
            tmp=calcTensorExpression('Delta{1,2,3}S{4,3}mu{2,4,5}',[1,5]);
            assert(isequal_(tmp,tmp2),"antipode right inverse error")

            % integral property
            %  δ_r*u=ε(u)δ_r, u*δ_l=ε(u)δ_l,
            % ∫_r(u_1)u_2=∫_r(u), u_1∫_l(u_2)=∫_l(u)

            intr=obj.getSC('intr').toMatrix();
            cointr=obj.getSC('cointr').toMatrix();
            intl=obj.getSC('intl').toMatrix();
            cointl=obj.getSC('cointl').toMatrix();
            % tmp=tensorprod(cointr,mu,1,1,Num=1);
            tmp=calcTensorExpression('cointr{1}mu{1,2,3}',[2,3]);
            % tmp2=tensorprod(ep,cointr,Num=1);
            tmp2=calcTensorExpression('ep{1}cointr{2}',[1,2]);
            assert(isequal_(tmp,tmp2),"right cointegral property error")
            % tmp=tensorprod(intr,Delta,1,2,Num=1);
            tmp=calcTensorExpression('intr{2}Delta{1,2,3}',[1,3]);
            % tmp2=tensorprod(intr,eta,Num=1);
            tmp2=calcTensorExpression('intr{1}eta{2}',[1,2]);
            assert(isequal_(tmp,tmp2),"right integral property error")
            % tmp=tensorprod(intl,cointr,1,2,Num=1);
            % tmp2=tensorprod(intl,cointr,Num=1);
            % issue:leftはまだ未実装
            tmp=[intr.'*cointr, intl.'*cointr; ...
                 intr.'*cointl, intl.'*cointl];
            assert(isequal_(tmp([1 2 4]),[1 1 1]),"integral normalization error")
            % intl*cointrのみCrIlになる．他はnormalizeされる
            disp("Confirmed to be a Hopf algebra")
            function assertT(expr1,ord1,expr2,ord2,msg)
                % assertT: assert with tensor order
                val1=calcTensorExpression(expr1,ord1);
                val2=calcTensorExpression(expr2,ord2);
                diff=val1-val2;
                
                if ~isequal(size(expr1),ord1) || ~isequal(size(expr2),ord2)
                    error("VectAlg:verifyHopf",msg)
                end
            end
        end
        function verifyInt(obj)
            isequal_=@(x,y)all(eqD(x,y,1e-5),"all");
            D=obj.dim;
            I=eye(D);
            mu=obj.getSC('prod');
            eta=obj.getSC('unit');
            Delta=obj.getSC('coprod');
            ep=obj.getSC('counit');
            S=obj.getSC('antipode');
            Si=obj.getSC('antipode_inv');

            intr=obj.getSC('intr');
            cointr=obj.getSC('cointr');
            intl=obj.getSC('intl');
            cointl=obj.getSC('cointl');
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
            intr=obj.getSC('intr');
            cointr=obj.getSC('cointr');
            intl=obj.getSC('intl');
            cointl=obj.getSC('cointl');
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
                disp("Empty StrAlg: ("+join(string(size(i1)),",")+")")
            else
                for ii=1:numel(i1)
                    disp_(i1(ii))
                end
            end
        end
        function disp_(i1)
            try
            feval("disp"+AlgebraConfig.H.disp_VectAlg,i1);
            catch ME
                warning(ME.identifier,'cannot display properly:\n %s',ME.message)
                i1.disp0
            end
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

            Z=obj.bs.ZERO{indexOp.Indices{2}};
            % Z.ZERO={Z};
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
            z=obj.set_c(obj.cf*0);
            if nargin==1
                varargin={1,1};
            end
            z = repmat(z,varargin{:});
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
