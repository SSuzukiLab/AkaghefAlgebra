classdef(InferiorClasses=?sym) StrEndV<StrAlg
    %UNTITLED このクラスの概要をここに記述
    %
    %{
      issue:bsがpropertyにあるとmethodの拡張が面倒であるため，できるかぎり同一のメソッドを用いるようにしたい
      その関係で，pwのformatを書き換える．cell {type,param} をNterm×Nrankだけ用意する
      今bsというパラメータが残っているが，それをpwに組み込むようにする
　　　　overrideする必要があるmethodはget2vrel,verify,verifyBase,pol,subs,combineTermである
　　　　combine termはHomの場合は無くても対応できそう．
　　　　StrEndV represents a linear operator on a tensor-product vector space.

DATA STRUCTURE
------------------------------------------------------------

1. cf : coefficient vector
   - Size: Nterm
   - cf[i] is the scalar coefficient of the i-th operator term.
   - Field type K (double, symbolic, etc.)

2. Operator descriptor (elementary operator)
   Each elementary operator is encoded as:

       { opKind , opArgs }

   where:
     - opKind : integer or enum identifying the operator kind
                (e.g., 1 = d/dx_i, 2 = multiply-by-x_i, ...)
     - opArgs : 1×M array storing operator-specific arguments
                (variable index, degree, flags, etc.)
       (exact format TBD)

   This pair is called an "OpNode".

3. A : operator layout matrix
   - A is a 2D array of size (Nterm × Nrank)
   - A[i][j] is an OpNode applied on the j-th tensor factor
     in the i-th term.
   - Thus, each row of A encodes a tensor-product operator term.

4. Overall representation:
       (cf , A)
   The program interprets this as:
       sum_i  cf[i] *  ⊗_{j=1…Nrank}  OpNode_ij

   Linearity is applied by expanding the operator on basis elements.

------------------------------------------------------------
TBD ITEMS
------------------------------------------------------------
- List and enumeration of opKind values (final operator catalog)
- Standardized structure and length of opArgs
- Memory layout of A (row-major vs. compressed form)
- Rules for simplifying/merging operator terms
- Definition of coefficient field K

    %}
    % .

    properties(Hidden,Constant)
        types=["id","d","x","xd","dd","th","qth"]
        NargRequired=[0 2 1 4 1 1 2]
    end
    properties
        bs
        varname=@(x)x
        dimV (1,:) =0
    end
    properties(Access=private,Transient)
        dispContents
    end

    methods

        function ret = StrEndV()
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
            ret.ctype="S";
        end
        function arg=castV(obj,arg)
            if isa(arg,'sym')||isnumeric(arg)
                arg=symp(arg(:),zeros(numel(arg),obj.dimV));
            end
        end
        function ret=or(i1,i2)
            dimV=[i1.dimV i2.dimV];
            ret=or@StrAlg(i1,i2);
            ret.dimV=dimV;
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
                ret=i1.set_cpb([i1.cf;i2.cf],[i1.pw;i2.pw],[i1.bs;i2.bs]).calc();
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
                obj=obj.set_cpb(obj.cf(idx),obj.pw(idx,:),obj.bs(idx,:));
            end
        end
        function obj=replace(obj,~)
        end
        function ret=getActMono(obj,idxTens,idxFactor)
            % ["id","d","x","xd","dd","th","qth"]
            opKind=obj.pw{idxTens}(idxFactor);
            opArgs=obj.bs{idxTens}{idxFactor};
            switch opKind
                case 1, ret=@id;
                case 2, ret=@d;
                case 3, ret=@m;
                case 4, ret=@xd;
                case 5, ret=@dd;
                case 6, error("unimplemented");
                case 7, ret=@qth;
                case 8
            end
            % 恒等演算子
            function [c,p]=id(p)
                c=1;
            end
            % q微分演算子
            function [c,p]=d(p)
                c=qNumS(opArgs{2}).n(p(opArgs{1}));
                p(opArgs{1})=p(opArgs{1})-1;
            end
            % 微分演算子
            function [c, p] = dd(p)
                c = p(opArgs{1});
                p(opArgs{1}) = p(opArgs{1}) - 1;
            end


            % M_ij の演算子
            function [c, p] = xd(p)
                c = qNumS(opArgs{4}).n(p(opArgs{2})) * opArgs{4}^(opArgs{3} * p.');
                p(opArgs{1}) = p(opArgs{1}) + 1;
                p(opArgs{2}) = p(opArgs{2}) - 1;
            end

            % 乗算演算子
            function [c, p] = m(p)
                c = 1;
                if length(opArgs)==2
                    power=opArgs{2};
                else
                    power=1;
                end
                p(opArgs{1}) = p(opArgs{1}) + power;
            end
            function [c, p] = qth(p)
                c = opArgs{3}^(opArgs{2}*p(opArgs{1}));
            end

            % 固有演算子の作用
            function [c, p] = th(p, f, eig)
                c = f(eig * p.'); % p, eig は横ベクトル
            end


        end
        function ret=act(obj,arg)
            arg=obj.castV(arg);
            if obj.term==0
                ret=0;
            elseif obj.term==1
                ret=arg;
                ret.cf=obj.cf*arg.cf;
                if obj.rank==1
                    Nfactor=length(obj.pw{1});
                    for i=Nfactor:-1:1
                        fun=obj.getActMono(1,i);
                        ret=ret.lfun(fun);
                        % disp(ret)
                    end
                else
                    Nfactor=cellfun(@length,obj.pw);
                    Nmax=max(Nfactor);
                    rank=obj.rank;
                    funCell=cell(Nmax,rank);
                    for i=1:Nmax
                        for j=1:rank
                            if Nfactor(j)>=i
                                funCell{i,j}=obj.getActMono(j,i);
                            else
                                funCell{i,j}=@(p)deal(1,p);
                            end
                        end
                    end
                    for i=Nmax:-1:1
                        C=[funCell(i,:);num2cell(obj.dimV)];
                        % fun=ret.tensorMor(sym(0),0,C{:});
                        % 指数をsymで扱う場合
                        fun=ret.tensorMor(sym(0),sym(0),C{:});
                        ret=ret.lfun(fun);
                    end
                end
            else
                ret=0;
                % tmp=obj;
                for i=1:obj.term
                    tmp=obj.set_cpb(obj.cf(i),obj.pw(i,:),obj.bs(i,:));
                    acted=tmp.act(arg);
                    ret=acted+ret;
                end
            end
        end
        function ret=unit(obj)
            ret=obj.make("id");
            ret=ret.set_cpb(1,repmat(ret.pw,1,obj.rank),repmat(ret.bs,1,obj.rank));
        end
        function ret=mtimes(i1,i2)
            ret=mtimes@StrAlg(i1,i2);
            ret.dimV=ret.dimV(1:end/2);
        end

        function ret=casttype(obj,arg)
            ret=obj.unit;
            ret.cf=arg;
        end
        function obj=updateDispContents(obj)
            contents=strings(obj.term,obj.rank);
            for j=1:numel(contents)
                opArgs=obj.bs{j};
                pw=obj.pw{j};
                V=obj.varname;
                content="";
                if isempty(pw)
                    content="id";
                end
                for i=1:length(pw)
                    opArgs_=opArgs{i};
                    switch pw(i)
                        case 1, tmp="id";
                        case 2
                            q=getQstr(opArgs_{2});
                            tmp="D"+q+"_"+V(opArgs_{1});
                        case 3, tmp="X_"+V(opArgs_{1});
                        case 4
                            q=getQstr(opArgs_{4});
                            eigv=opArgs_{3};
                            idx=find(eigv);
                            eigs=string(eigv);
                            if numel(idx)==0
                                qth="";
                            elseif isscalar(idx)
                                qth=q+"^{"+eigs(idx)+"θ_"+V(idx)+"}";
                            else
                                if q=="q", q2="q^"; else, q2=q; end
                                qth=q2+"{"+join(eigs(idx)+"θ_"+V(idx),"+")+"}";
                            end
                            tmp="X_"+V(opArgs_{1})+"D"+q+"_"+V(opArgs_{2})+qth;
                        case 5, tmp="∂_"+V(opArgs_{1});
                        case 6, tmp="θ_"+V(opArgs_{1});
                        case 7
                            q=getQstr(opArgs_{3});
                            tmp=sprintf(q+"^{%sθ_%s}",string(opArgs_{2}),string(V(opArgs_{1})));
                        case 8
                        otherwise
                            tmp="?";
                    end
                    content=content+" "+tmp;
                end
                contents(j)=content;
            end
            obj.dispContents=contents;

            function q=getQstr(arg)
                if isnumeric(arg)
                    q="q";
                else
                    q=string(arg);
                end
            end
        end
        function disp1(obj)
            % 今のところdisp1のみ対応させている．
            obj=updateDispContents(obj);
            obj.pw=num2cell(reshape(1:numel(obj.pw),size(obj.pw)));
            obj.disp1@StrAlg()
        end
        function ret=convertTermToSym(obj)
            % disp2はまだちゃんとできていない
            obj=updateDispContents(obj);
            ret=obj.dispContents;

        end

        function disp_(i1)
            try
                feval("disp"+AlgebraConfig.H.disp_StrEndV,i1);
            catch ME
                warning(ME.identifier,'cannot display properly:\n %s',ME.message)
                i1.disp0
            end
        end
        function content=convertBaseString(obj,pw,bs)
            % ["id","d","x","xd","dd","th","qth"]
            content=obj.dispContents(pw);
        end
        function obj=set_cpb(obj,varargin)
            obj=obj.set_cp(varargin{1:end-1});
            obj.bs=varargin{end};
        end
        function ret=make(obj,type,arg)
            arguments
                obj
                type (1,1) {mustBeMember(type,["id","d","x","xd","dd","th","qth"])}
            end
            arguments(Repeating)
                arg
            end
            pw=getType(type);
            assert(obj.NargRequired(pw)<=length(arg),'Too few arguments')
            ret=obj.set_cp(1,{[{pw},arg]});
        end
    end
    methods(Static)
        function ret=makeV(dimV)
            ret=StrEndV;
            ret.dimV=dimV;
            ret.ZERO={ret};
            ret.cf=1;
            ret.bs={{}};

        end
    end
end


function ret=getType(typ)
    ret=find(["id","d","x","xd","dd","th","qth"]==string(typ));
end