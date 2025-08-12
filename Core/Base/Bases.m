classdef Bases<matlab.mixin.Heterogeneous&handle
    properties(Access=protected)
        dim_ (1,1) {mustBeInteger} 
        var_ (1,:) string %names of each basis vector
    end
    properties
        name (1,:) char='basis' %name of basis
        ZERO strAlg
    end
    properties
        ctype NumericType =NumericType.D 
        ptype NumericType =NumericType.D
    end
    methods
        function obj=Bases(dim,var,name)
            % GVBC (dim_,varname)
            if nargin==0, return; end
            obj=obj.setVar(dim,var);
            if nargin==3
                obj.name=name;
            end
        end
    end
    methods(Sealed)
        function ret=dim(obj)
            ret=sum([obj.dim_]);
        end
        function ret=dims(obj)
            ret=[obj.dim_];
        end
        function obj=setVar(obj,dim,var)
            obj.dim_=dim;
            if numel(var)==1&&dim>1
                var=var+"_"+(1:dim);                                                                                                                                                                                                                                                        
            end
            obj.var_=var;
        end
        function ret=getCtype(obj)
            ret=getType([obj.ctype]);
        end
        function ret=validate(obj,arg)
            ret=isa(arg,getCtype(obj).type)||isequal(class(arg),'double');
        end
        function ret=string(obj)
            ret=horzcat(obj.var_);
            if isempty(ret)
                ret="";
            end
        end
        function disp(obj)
            if numel(obj)==0
                disp("Empty bases")
                return
            end
            disp("Basis:["+join({obj.name},", ")+"]")
            disp("["+join(string(obj),", ")+"]")
        end
    end
end