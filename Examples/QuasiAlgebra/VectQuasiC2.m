classdef VectQuasiC2<VectQuasiHopfAlg

    properties
        N=2
        bs2=Bases(2,["p0","p1"],"Z/2Z")
    end
    methods(Static)
        function [Z,g]=getGenerator()
            Z=VectQuasiC2();
            Z=Z.setBase(makeBase(2));
            Z.setConst;
            g=Z.make(1,2);
        end
        
    end
    methods

        function setConst(obj)
            % Sweedler Hopf algebra H4 structure constants
            % g:generator g^N=1
            % 
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
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
            % associator
            Phi=zeros(N,N,N);
            Phi_inv=zeros(N,N,N);
            for i=0:N-1
                for j=0:N-1
                    for k=0:N-1
                        Phi(i+1,j+1,k+1)=exp(2*pi*1i*i*j*k/N);
                        Phi_inv(i+1,j+1,k+1)=exp(-2*pi*1i*i*j*k/N);
                    end
                end
            end
            alpha=zeros(N,1);
            beta=zeros(N,1);
            alpha(2)=1;
            beta(1)=1;
            obj.setSCQuasi(associator=Phi,associator_inv=Phi_inv,alpha_S=alpha,beta_S=beta);

        end
    end
end
function ret=makeBase(N)
    name=['Z/' num2str(N) 'Z'];
    var="e"+(0:N-1);
    ret=Bases(N,var,name);
end