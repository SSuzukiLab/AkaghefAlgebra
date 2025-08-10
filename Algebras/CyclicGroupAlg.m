classdef(InferiorClasses=?sym) CyclicGroupAlg<VectAlg
    % SWEEDLERALG 4-dimensional Sweedler Hopf algebra
    properties(Constant)
        bs0= TypeParam(@makeBase)
    end
    properties
        N
    end
    methods(Static)
        function [Z,One]=getGenerator(N)
            Z=CyclicGroupAlg();
            Z=Z.setBase(Z.bs0.get(N));
            Z.N=N;
            Z.setConst;
            One=Z.make(1,2);
        end
    end
    methods
        % function obj=SweedlerAlg()
        % end
        function obj=casttype(obj,arg)
            if isa(arg,'double')&&isscalar(arg)
                obj=obj.set_c([arg;0;0;0]);
            else
                error not_impl
            end
        end
        function setConst(obj)
            % Sweedler Hopf algebra H4 structure constants

            % Basis: K^i E^j (i,j) i,j=0,1,...,N-1
            % E^N=0, K^N=1

            % Multiplication tensor M(i,j,k): e_i * e_j = sum_k M(i,j,k) * e_k

            N=obj.N;
            M = zeros([N N N]);
            for i=0:N-1
                for j=0:N-1
                    M(i+1,j+1,mod(i+j,N)+1)=1;
                end
            end

            % Comultiplication tensor C(i,j,k): Δ(e_i) = sum_{j,k} C(i,j,k) * (e_j ⊗ e_k)
            % Δ(g)=g⊗g
            C = zeros([N N N]);
            for i=0:N-1
                C(i+1,i+1,i+1) = 1;
            end

            % Counit vector: epsilon_i = ε(e_i)
            % ε(e_n)=1
            epsilon = ones(N,1);

            % Unit vector: η(1) = η^i(e_i)
            eta = zeros(N,1);
            eta(1) = 1;
            % Antipode matrix: S(i,j) = coefficient of e_j in S(e_i)
            % S(e_n)=e_{-n}
            S = zeros([N N]);
            for i=0:N-1
                S(i+1,mod(-i,N)+1)=1;
            end
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
            integral=[1;zeros(N-1,1)];
            cointegral=deal(ones(N,1));
            obj.setIntegrals(integral,cointegral,integral,cointegral);
        end
    end
end
function ret=makeBase(N)
    name=['Z/' num2str(N) 'Z'];
    var="e"+(0:N-1);
    ret=Bases(N,var,name);
end