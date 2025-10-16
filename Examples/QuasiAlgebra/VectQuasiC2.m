classdef VectQuasiC2<VectQuasiHopfAlg

    properties(Constant)
        N=2
        bs2=Bases(2,["p+","p-"],"Z/2Z")
    end
    methods(Static)
        function [Z,g]=getGenerator(type)
            arguments
                type {mustBeMember(type,1:2)}
            end
            Z=VectQuasiC2();
            if type==1
                Z=Z.setBase(makeBase(2));
                Z.setConst(1);
                g=Z.make(1,2);
            else
                Z=Z.setBase(makeBase3());
                Z.setConst(2);
                g=[Z.make(1,1),Z.make(1,2)];
            end
        end
        
    end
    methods

        function setConst(obj,type)
            % Sweedler Hopf algebra H4 structure constants
            % g:generator g^N=1
            % 
            N=obj.N;
            if type==1
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
                % associator
                Psi=reshape([3, 1, 1, -1, 1, -1, -1, 1]/4,[N,N,N]);
                Psi_inv=Psi;
                % for i=0:N-1
                %     for j=0:N-1
                %         for k=0:N-1
                %             Psi(i+1,j+1,k+1)=(-1)^(i*j*k)/4;
                %             Psi_inv(i+1,j+1,k+1)=(-1)^(-i*j*k)/4; % まちがい
                %         end
                %     end
                % end
                alpha=zeros(N,1);
                beta=zeros(N,1);
                alpha(2)=1;
                beta(1)=1;
            else
                M = zeros([N N N]);
                M(1,1,1)=1;
                M(2,2,2)=1;
    
                % Comultiplication tensor C(i,j,k): Δ(e_i) = sum_{j,k} C(i,j,k) * (e_j ⊗ e_k)
                % Δ(g)=g⊗g
                C = zeros([N N N]);
                C(1,1,1)=1;
                C(1,2,2)=1;
                C(2,2,1)=1;
                C(2,1,2)=1;
    
                % Counit vector: epsilon_i = ε(e_i)
                % ε(e_n)=1
                epsilon = [1;0];
    
                % Unit vector: η(1) = η^i(e_i)
                eta = [1;1];
                % Antipode matrix: S(i,j) = coefficient of e_j in S(e_i)
                % S(e_n)=e_{-n}
                S=eye(N);

                obj.setSC(obj.identifier,M,eta,C,epsilon,S);
                % associator
                Psi=ones(N^3,1);
                Psi(end)=-1;
                Psi=reshape(Psi,[N,N,N]);
                Psi_inv=Psi;

                alpha=[1,-1];
                beta=[1,1];
            end
            obj.setSCQuasi(associator=Psi,associator_inv=Psi_inv,alpha_S=alpha,beta_S=beta);

        end
    end
end
function ret=makeBase(N)
    name=['Z/' num2str(N) 'Z'];
    var="e"+(0:N-1);
    ret=Bases(N,var,name);
end
function ret=makeBase3()
    name=['h2'];
    var=["p+","p-"];
    ret=Bases(2,var,name);
end