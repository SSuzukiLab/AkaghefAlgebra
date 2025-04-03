classdef(InferiorClasses=?sym) SweedlerAlg<VectAlg
    % SWEEDLERALG 4-dimensional Sweedler Hopf algebra
    properties(Constant)
        bs0=Bases(4,["1" "g" "x" "x*g"])
    end
    methods(Static)
        function [g,x]=getGenerator()
            Z=SweedlerAlg();
            Z.ZERO={Z};
            Z.setConst;
            g=Z.make(1,2);
            x=Z.make(1,3);
        end
    end
    methods
        function obj=SweedlerAlg()
            obj.cf=zeros(4,1);
            obj.bs=obj.bs0;
        end
        function obj=casttype(obj,arg)
            if isa(arg,'double')&&isscalar(arg)
                obj=obj.set_c([arg;0;0;0]);
            else
                error not_impl
            end
        end
        function setConst(obj)
            % Sweedler Hopf algebra H4 structure constants

            % Basis: e1 = 1, e2 = g, e3 = x, e4 = gx

            % Multiplication tensor M(i,j,k): e_i * e_j = sum_k M(i,j,k) * e_k
            M = zeros(4,4,4);
            M(1,1,1) = 1;
            M(1,2,2) = 1; M(2,1,2) = 1;
            M(1,3,3) = 1; M(3,1,3) = 1;
            M(1,4,4) = 1; M(4,1,4) = 1;
            M(2,2,1) = 1;
            M(2,3,4) = -1; M(3,2,4) = 1;
            M(2,4,3) = -1; M(4,2,3) = 1;

            % Comultiplication tensor C(i,j,k): Δ(e_i) = sum_{j,k} C(i,j,k) * (e_j ⊗ e_k)
            C = zeros(4,4,4);
            C(1,1,1) = 1;
            C(2,2,2) = 1;
            C(3,3,1) = 1; C(3,2,3) = 1;
            C(4,4,2) = 1; C(4,1,4) = 1;

            % Counit vector: epsilon(i) = ε(e_i)
            epsilon = [1; 1; 0; 0];

            % Counit vector: η(1) = η^i(e_i)
            eta = [1; 0; 0; 0];

            % Antipode matrix: S(i,j) = coefficient of e_j in S(e_i)
            S = zeros(4,4);
            S(1,1) = 1;
            S(2,2) = 1;
            S(3,4) = -1;
            S(4,3) = -1;
            
         
            
            obj.setSC(class(obj),M,eta,C,epsilon,S);
        end
    end
end
