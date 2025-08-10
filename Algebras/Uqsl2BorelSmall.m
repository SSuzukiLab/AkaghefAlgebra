classdef(InferiorClasses=?sym) Uqsl2BorelSmall<VectAlg
    % SWEEDLERALG 4-dimensional Sweedler Hopf algebra
    properties(Constant)
        bs0= TypeParam(@makeBase)
    end
    properties
        N
        zN
    end
    methods(Static)
        function [Z,K,E]=getGenerator(N)
            Z=Uqsl2BorelSmall();
            Z=Z.setBase(Z.bs0.get(N));
            Z.N=N;
            Z.zN=exp(2*pi*1i/N);
            Z.setConst;
            K=Z.make(1,2);
            E=Z.make(1,N+1);
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
            z=obj.zN;
            M = zeros(N*[1 1 1 1 1 1]);
            for i1=0:N-1
                for i2=0:N-1
                    for j1=0:N-1
                        for j2=0:N-1
                            if j1+j2>=N
                                continue
                            end
                            M(i1+1,j1+1,i2+1,j2+1,mod(i1+i2,N)+1,j1+j2+1) = z^(-j1*i2);
                        end
                    end
                end
            end
            M=reshape(M,N^2*[1 1 1]);

            qN=qNumA(z);
            % Comultiplication tensor C(i,j,k): Δ(e_i) = sum_{j,k} C(i,j,k) * (e_j ⊗ e_k)
            % Δ(E)=1⊗E+E⊗K, Δ(K)=K⊗K, 
            C = zeros(N*[1 1 1 1 1 1]);
            for i=0:N-1
                for j=0:N-1
                    for k=0:j
                        C(i+1,j+1,i+1,k+1,mod(i+k,N)+1,j-k+1) = qN.nchoosek(j,k)*z^(-k*(j-k));
                    end
                end
            end
            C=reshape(C,N^2*[1 1 1]);

            % Counit vector: epsilon_i = ε(e_i)
            % ε(E)=0, ε(K)=1
            epsilon = [ones(N,1); zeros(N^2-N,1)];

            % Counit vector: η(1) = η^i(e_i)
            eta = zeros(N^2,1);
            eta(1) = 1;
            % Antipode matrix: S(i,j) = coefficient of e_j in S(e_i)
            % S(E)=-EK^-1, S(K)=K^-1
            S = zeros(N*[1 1 1 1]);
            for i=0:N-1
                for j=0:N-1
                    S(i+1,j+1,mod(-i-j,N)+1,j+1)=(-1)^j*z^(i*j+j*(j+1)/2);
                end
            end
            S=reshape(S,N^2*[1 1]);
            
            
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
            [Ir,Cr,Il,Cl]=deal(zeros(N,N)*z);
            Ir(2,N)=z^-1;
            Cr(:,N)=z.^(0:N-1);
            Il(1,N)=z^1;
            Cl(:,N)=z.^(2*(0:N-1)-1);
            Ir=reshape(Ir,[N^2 1]);
            Cr=reshape(Cr,[N^2 1]);
            Il=reshape(Il,[N^2 1]);
            Cl=reshape(Cl,[N^2 1]);
            % Cl=(S^-1)*Cr;
            % Il=S.'*Ir;
            obj.setIntegrals(Ir,Cr,Il,Cl);
            % obj.setIntegrals()
        end
    end
end
function ret=makeBase(N)
    name=['Uqsl2BorelSmall_ζ' num2str(N)];
    var=reshape(["" "K^"+(1:N-1)].'+["" "E^"+(1:N-1)],[N^2 1]);
    var(1)="1";
    ret=Bases(N^2,var,name);
end