classdef(InferiorClasses=?sym) StrCnmodalg2<StrAlg
    properties(Constant,Hidden)
        B=TypeParam(@(l)Bases(2*l,["x_m"+(l:-1:1) "x_"+(1:l)],"StrC"+l+"alg"))
    end
    properties
        l
    end
    properties(Hidden,Dependent)
        dict
    end
    methods
        function o=times(i1,i2)
            o=commeval(i1|i2);
        end
        function obj=StrCnmodalg2(l)
            % コンストラクタ　C_l型の代数を生成する
            obj.priority=1:2*l;
            obj.l=l;
            obj.ctype="S";
        end
        function ret=get.dict(obj)
            l=obj.l;
            ret=dictionary([-l:-1 1:l],1:2*l);
        end
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            % get2vRelation 二次の交換関係を実装
            persistent S
            try
                rel=S{obj.l}.rel;
                mlist=S{obj.l}.mlist;
                comm=S{obj.l}.comm;
                inv=S{obj.l}.inv;
            catch
                if isempty(S)
                    S={};
                end
                ll=obj.l;
                D=obj.dict;
                q=sym('q');
                S(ll)={Struct};
                S{ll}.rel=obj.empty;
                O=StrCnmodalg2.getGenerator(ll);
                % x_i,x_j q-commutative
                for i=-ll:ll
                    for j=-ll:ll
                        if i==-j||i>=j||i==0||j==0
                            continue
                        end
                        S{ll}.rel(end+1)=O.make([-1 q],{D([j i]) D([i j])});
                    end
                end
                % x_i x_{-i}
                for i=1:ll
                    [cf,pw]=Omega(i+1);
                    S{ll}.rel(end+1)=O.make([-1 q^2 q^2*(q-q^-1)*cf],[{D([i -i]) D([-i i])} pw]);
                end
                S{ll}.comm=[nan;nan];
                S{ll}.inv=[nan;nan];
                S{ll}=obj.get2vRelation_(S{ll});

                [rel,mlist,comm,inv]=get2vRelation(obj);
            end
            function [cf,pw]=Omega(i)
                cf=sym.empty;
                pw=cell(0);
                for k=i:ll
                    cf(end+1)=q^(k-i);
                    pw(end+1)={D([-k k])};
                end
            end
        end

        function arg=casttype(obj,arg)
            arg=obj.make(arg,{[]});
        end
        
        function obj=make(obj,cf,pw)
            obj=obj.make@StrAlg(cf,pw,StrCnmodalg2.B.get(obj.l));
        end
    end
    methods(Static)
        function [O,varargout]=getGenerator(l,opt)
            arguments
                l {mustBeInteger}
                opt String {mustBeMember(opt,["all","dict"])}="dict"
            end
            O=StrCnmodalg2(l).make(0,{[]});
            X=cell(0);
            for i=1:2*l
                X{end+1}=O.make(1,{i});
            end
            if opt=="all"
                varargout=X;
            else
                X=[num2cell([-l:-1 1:l]);X];
                varargout={dictionary(X{:})};
            end
        end
    end
end
