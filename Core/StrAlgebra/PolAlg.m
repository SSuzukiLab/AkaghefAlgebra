classdef (InferiorClasses={?sym})PolAlg
    % PolAlg 多項式クラス
    properties(SetAccess=protected,Dependent)
        term
        dim
    end

    properties
        cf (:,1)
        pw
        base (1,:) Bases
        spec (1,1) SpaceSpec % specify principle of space
        % also, store structure constant 
        ctype NumericType ="D"
        ptype NumericType ="D"
    end
    methods

        % コンストラクタ
        function obj=PolAlg(i1,i2)
            if nargin==0
                return
            elseif isa(i1,"PolAlg")
                % error("convert!")
                tmp=i1;
                % obj.ctype=tmp.ctype;
                % obj.ptype=tmp.ptype;
                i1=tmp.cf;
                i2=tmp.pw;
            elseif nargin==1
                i1=sum(i1);
                i2=zeros(1,0);

            end
            obj.base=Bases(size(i2,2),"p");
            % obj.dim=size(i2,2);
            if isa(i1,'sym')
                obj.ctype="S";
            end
            if isa(i2,'sym')
                obj.ptype="S";
            end
            obj.cf=obj.ctype.zero+i1;
            obj.pw=obj.ptype.zero+i2;
            obj=obj.C();
        end
        %% cast type

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


        %% binomial operation
        function ret=plus(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            try
                sz=size(zeros(size(i1))+zeros(size(i2)));
            catch
                error("PolAlg:plus","size dimensions must match")
            end
            if any(sz==0)
                ret=repmat(i1,sz);
            else
                % assert(configure(i1,i2))
                assert(isequal(i1.base,i2.base))
                ret=arrayfun(@(x,y)plus_(x,y),repmat(i1,sz./size(i1)),repmat(i2,sz./size(i2)));
            end
            % end
        end
        function obj=set.base(obj,b)
            if obj.dim<b.dim
                obj.pw(:,obj.dim+1:b.dim)=obj.ptype.zero;
            end
            obj.base=b;
            obj.ctype=getType([obj.base.ctype]);
            obj.ptype=getType([obj.base.ptype]);

            if b.dim==0
                obj.base=eval(class(b)+".empty(1,0)");
            end
            % obj.dim=b.dim;
        end
        function ret=plus_(i1,i2)
            % PLUS_ supplementary method for scalar
            [i1,i2]=alignNum(i1,i2);
            ret=i1.set_cp([i1.cf;i2.cf],[i1.pw;i2.pw]).C();
        end
        function ret=minus(i1,i2)
            ret=i1+(-i2);
        end
        function i1=uminus(i1)
            i1.cf=-i1.cf;
        end

        function i1=uplus(i1)
        end
        function ret=eq(i1,i2)
            sub=C(i1-i2);
            ret=all(sub.cf==0);
        end
        %乗算,作用
        function ret=mtimes(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            N=i1.term;
            M=i2.term;
            D=i1.dim;
            C=repmat(i1.ctype.zero,N*M,1);
            P=repmat(i1.ptype.zero,N*M,D);
            for ii=1:N
                for jj=1:M
                    C((ii-1)*M+jj)=i1.cf(ii)*i2.cf(jj);
                    P((ii-1)*M+jj,:)=i1.pw(ii,:)+i2.pw(jj,:);
                end
            end
            ret=i1.set_cp(C,P).C;
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
                ret=i1.set_cp(i1.cf^i2,i1.pw*i2);
            end
        end
        function ret=unit(i1)
            ret=i1.set_cp(1,zeros(1,i1.base.dim));
        end
        %テンソル積
        function o=or(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            o=i1.lfun(@fun);
            % o.dim=i1.dim+i2.dim;
            o.base=[i1.base i2.base];
            function [c,p]=fun(p)
                c=i2.cf;
                p_=p;
                p=[repmat(p,i2.term,1) i2.pw];
            end
            o.ctype=getType([i1.ctype i2.ctype]);
        end
        %副積
        function o=and(i1,i2)
            o=~(i1|i2);
        end
        %副積計算
        %▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
        %▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
        %▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
        %▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
        function o=not(i1)
            o=i1.prodeval;
        end

        function o=prodeval(i1)
            b=i1.base(1:end/2);
            d=b.dim;
            o=i1;
            o.pw=o.pw(:,1:d)+o.pw(:,d+1:2*d);
            o.base=b;
        end

        % 各項への作用
        function ret=lfun(obj,fun)
            CC=obj.ctype.zeros(0,1);
            PP=obj.ptype.zeros(0,0);
            for ii=1:obj.term
                [cc,pp]=fun(obj.pw(ii,:));
                CC=[CC;cc*obj.cf(ii)];
                PP=[PP;pp];
            end
            ret=obj.set_cp(CC,PP).C();
        end
        %　簡約化=同次項括り＋零係数項削除
        function obj=C(obj)
            obj=combineTerm(obj);
            obj=obj.removeZero();
            % obj.term=numel(obj.cf);
        end
        % 零係数項削除,係数簡約化ステップ
        function obj=removeZero(obj)
            try
                idx=~isAlways(obj.cf==0,Unknown="false");
                % idx=abs(subs(i1.cf,retq(),0.71))>0.000001;
                % idx=1:length(i1.cf);
                % idx=find(i1.cf);
                % s=@(x)x;
                % s=@simplify
                obj=obj.set_cp(simplify(obj.cf(idx)),obj.pw(idx,:));
            catch
                idx=abs(obj.cf)>1000*eps(obj.cf);
                obj=obj.set_cp(obj.cf(idx),obj.pw(idx,:));
            end
        end
        function ret=split(arg,dims)
            assert(length(dims)==length(arg.base))
            ret=cellfun(@(p)arg.set_cp(1,p),mat2cell(arg.pw,ones(1,arg.term),dims));
            for i=1:size(ret,1)
                for j=1:size(ret,2)
                    ret(i,j).base=arg.base(j);
                end
                ret(i,1).cf=arg.cf(i);
            end
        end

        %% display function
        % objの変更
        function i1=set_cp(i1,i2,i3)
            % assert(isa(i2,i1.ctype))
            % assert(isa(i3,i1.ptype))
            i1.cf=i2;
            i1.pw=i3;
            if isempty(i2)
                i1.cf=i1.ctype.zero;
                i1.pw=i1.ptype.zeros(1,i1.dim);
            end
        end

        function ret=testFunc(i1,i2,i3)
            [ret,ret1]=helpdiv(i1,i2,i3)
        end
        function ret=pol(arg)
            if arg.dim~=0
                vars=sym(arg.base.string);
                G=prod(vars.^arg.pw,2);
                ret=arg.cf.*G;
            else
                ret=sym(arg.cf);
            end
        end
        %% 表示
        function disp(i1)
            if isempty(i1)
                disp("Empty PolAlg: ("+join(string(size(i1)),",")+")")
            else
                for ii=1:numel(i1)
                    disp_(i1(ii))
                end
            end
        end
        function disp_(i1)
            try
            feval("disp"+AlgebraConfig.H.disp_PolAlg,i1);
            catch ME
                warning(ME.identifier,'cannot display properly:\n %s',ME.message)
                i1.disp0
            end
        end
        function disp0(i1)
            builtin("disp",i1)
        end


        % テーブル形式表示
        function disp1(i1)
            N=numel(i1);
            for ii=1:N
                if isempty(i1(ii).cf)
                    disp(table(0,VariableNames="coeff"))
                    return
                end
                % Example string array with duplicates
                strArray = i1(ii).base.string;

                % Find unique elements and their counts
                [uniqueElements, ~, idx] = unique(strArray, 'stable');
                counts = histcounts(idx, 1:numel(uniqueElements) + 1);

                % Initialize a cell array to store modified elements
                modifiedStrArray = strings(size(strArray));

                % Loop through each unique element and add index to duplicates
                for i = 1:numel(uniqueElements)
                    element = uniqueElements(i);
                    occurrences = find(strcmp(strArray, element));  % Indices of occurrences

                    for j = 1:numel(occurrences)
                        if counts(i) > 1  % If the element is duplicated
                            modifiedStrArray(occurrences(j)) = sprintf('%s_%d', element, j);
                        else
                            modifiedStrArray(occurrences(j)) = element;  % Keep original if no duplicates
                        end
                    end
                end
                pp=i1(ii).pw;
                try
                    disp(cell2table([num2cell(i1(ii).cf),num2cell(i1(ii).pw)], ...
                        VariableNames=["coeff" modifiedStrArray]));
                catch ME
                    warning(ME.identifier,'%s',ME.message)
                    i1(ii).disp0
                end
            end
        end
        % 数式の形式の表示
        function disp2(i1)
            Name=inputname(1);
            if isempty(Name)
                Name='ans';
            end
            if numel(i1)==1
                disp(Name==sum(i1.pol))
            else
                disp(sum(i1.pol))
            end
        end
        % 多数項の表示
        function disp3(i1)
            assert(numel(i1)==1)
            disp()
        end
        function ret=convertTermToSym(obj)
            assert(numel(obj)==1)
            try
                ret=cellfun(@(x,y)x*prod(sym(string(obj.base)).^y),num2cell(sym(obj.cf)),mat2cell(obj.pw,ones(1,obj.term),obj.dim));
            catch ME
                error("PolAlg:disp","cannot convert coefficient to sym")
            end
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
        %複製
        function ret=matrix(obj,i1)
            ret=repmat(obj,i1);
        end
        function ret=subs(obj,varargin)
            ret=obj;
            ret.cf=subs(obj.cf,varargin{:});
            ret.pw=subs(obj.pw,varargin{:});
        end


        function ret=get.dim(obj)
            ret=size(obj.pw,2);
        end
        function ret=get.term(obj)
            ret=length(obj.cf);
        end
        function ret=getB(obj)
            ret=obj.base;
        end
    end
    methods(Hidden)
        function ret=dimP(obj)
            ret=size(obj.pw,2);
        end
    end
    %zeros
    methods (Static)
        function z = zeros(varargin)
            if (nargin == 0)
                z = PolAlg;
            elseif any([varargin{:}] <= 0)
                z = PolAlg.empty(varargin{:});
            else
                z = repmat(PolAlg,varargin{:});
            end
        end
        function F=tensorMor(c0,p0,func,dim)
            arguments
                c0
                p0
            end
            arguments(Repeating)
                func
                dim
            end
            N=numel(func);
            F=@fun;
            dims=cell2mat(dim);
            idx=mat2cell(1:sum(dims),1,dims);
            function [c,p]=fun(p)
                C=cell(1,N);
                P=cell(1,N);
                for i=1:N
                    [C{i},P{i}]=func{i}(p(idx{i}));
                end
                terms=cellfun(@numel,C);
                idx0=terms==0;
                [C{idx0}]=deal(0);
                [P{idx0}]=deal(zeros(1,0));
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
                c=repmat(c0,h,1);
                p=repmat(p0,h,sum(dims));
                for i=1:h
                    c(i)=prod([TC{i,:}]);
                    p(i,:)=horzcat(TP{i,:});
                end



            end
        end
    end

end


%% Local method
% 同類項括り
function ret=combineTerm(i1)
    c=i1.cf;
    p=i1.pw;
    if isempty(p)
        ret=i1.set_cp(sum(c),p);
        return
    end
    % [u,uidx]=unique(p,"rows");
    % [v,vidx]=sortrows(p);
    % jj=0;
    % N=length(uidx);
    % S=zeros(N,1,i1.ctype);
    % for ii=1:length(vidx)
    %     if jj<N && uidx(jj+1)==vidx(ii)
    %         jj=jj+1;
    %     end
    %     S(jj)=S(jj)+c(vidx(ii));
    % end
    % ret=i1.set_cp(S,u);
    [v,vidx]=sortrows(p);
    S=repmat(i1.ctype.zero,numel(vidx),1);
    c=c(vidx);
    uidx=zeros(1,numel(vidx));
    N=1;
    v_=v(1,:);
    for ii=1:numel(vidx)
        if ~isequal(v_,v(ii,:))
            N=N+1;
            v_=v(ii,:);
        end
        uidx(N)=ii;
        S(N)=S(N)+c(ii);
    end
    uidx(uidx==0)=[];
    ret=i1.set_cp(S(1:N),v(uidx,:));

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

function [i1,i2]=alignNum(i1,i2)
    arguments
        i1 PolAlg
        i2 PolAlg
    end
    tf=logical([numel(i1.base),numel(i2.base)]);
    if all(~tf)

    elseif all(tf)
        assert(isequal(i1.base,i2.base))
    elseif tf(2)
        i1.base=i2.base;
        i1.pw=repmat(i1.ptype.zero,i1.term,i1.base.dim);
    elseif tf(1)
        i2.base=i1.base;
        i2.pw=repmat(i2.ptype.zero,i2.term,i2.base.dim);

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
