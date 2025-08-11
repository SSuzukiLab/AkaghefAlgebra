classdef(InferiorClasses=?sym) symWeylXQ<symp
    properties(Constant,Hidden)
        B=TypeParam(@(N)Bases(3*N,reshape(["X","qt" "th"]+(1:N)',1,3*N),"weyl_"+N))
    end
    properties
        Nvar
        q=sym('q')
    end

    %% generation
    methods(Static)
        function v=getGenerator(l)
            v=rowfun(@(x)symWeylXQ(1,x),table(eye(3*l)));
            v=dictionary((1:3*l)',v.Var1);
        end
    end
    methods
        function obj=symWeylXQ(varargin)
            obj@symp(varargin{:})
            obj.Nvar=obj.dim/3;
            mustBeInteger(obj.Nvar)
            B=obj.B.get(obj.Nvar);
            B.ctype="S";
            % B.ptype="S";
            obj.base=B;
        end
        function ret=mtimes(i1,i2)
            expr=i1|i2;
            Nvar=expr.Nvar;
            P=expr.pw(:,Nvar+(1:Nvar)).*expr.pw(:,3*Nvar+(1:Nvar));
            ret=expr.set_cp(expr.cf.*expr.q.^(sum(P,2)), ...
                expr.pw(:,1:end/2)+expr.pw(end/2+1:end));
            ret.base=ret.base(1:end/2);
        end
        function obj=unit(obj)
            obj=1;
        end
        function ret=mrdivide(i1,i2)
            assert(isa(i2,"symWeylXQ"))
            assert(i2.term==1,'除算が定義されません')
            Nvar=i2.Nvar;
            i2=i2.set_cp(i2.cf,-i2.pw);
            ret=i1*i2;
        end
        function ret=or(i1,i2)
            ret=or@symp(i1,i2);
            try
                ret.Nvar=[i1.Nvar i2.Nvar];
            catch ME
                if ~isa(i1,"symWeylXQ")
                    ret.Nvar=[i2.Nvar i2.Nvar];
                elseif ~isa(i2,"strWeylXQ")
                    ret.Nvar=[i1.Nvar i1.Nvar];
                else
                    rethrow(ME)
                end
            end
        end
        function ret=act(i1,i2)
            if ~isa(i2,"symp")
                i2=symp(i2,zeros(1,i1.Nvar));
            else
                assert(i2.dim==i1.Nvar);
            end
            assert(i2.dim==sum(i1.Nvar))
            splitted1=split(i1,3*i1.Nvar);
            splitted2=split(i2,i1.Nvar);
            

        end
    end
end