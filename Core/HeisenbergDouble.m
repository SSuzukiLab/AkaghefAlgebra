classdef HeisenbergDouble<VectAlg
    % HEISENBERGDOUBLE heisenberg double of Hopf algebra
    properties
        H1
        H2
    end
    methods(Static)
        function [g,x]=getGenerator(H1,H2)
            Z=HeisenbergDouble();
            Z.ZERO={Z};
            Z.H1=H1;
            Z.H2=H2;
            Z.setConst;
            g=Z.make(1,Z.bs,2);
            x=Z.make(1,Z.bs,3);
        end
    end
    methods
        function obj=HeisenbergDouble()
           
        end
        funnction ret=casttype(obj,arg)
            if ~isa(arg,'VectAlg')
                
            else
                error not_impl
            end
        end
        function obj=casttype(obj,arg)
            if isa(arg,'double')&&isscalar(arg)
                obj=obj.set_c([arg;0;0;0]);
            else
                error not_impl
            end
        end
    end
end