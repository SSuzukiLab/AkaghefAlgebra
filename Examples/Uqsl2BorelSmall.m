classdef(InferiorClasses=?sym) Uqsl2BorelSmall<VectAlg
    % uqsl2borelsmall borel subalgebra of Uq(sl2)
    % style=R: Δ(E)=1⊗E+E⊗K, style=L: Δ(E)=E⊗1+K⊗E
    properties(Constant)
        bs0= TypeParam(@makeBase)
    end
    properties
        N
        M
        q
        style {mustBeMember(style,["L","R"])}="R"
    end
    methods(Static)
        function [Z,K,E]=getGenerator(ratio,style,arg)
            arguments
                ratio % =M/N 
                style {mustBeMember(style,["L","R"])}
                arg.qtype {mustBeMember(arg.qtype,["complex","sym"])} = "complex"
                arg.q
            end
            [M,N]=rat(ratio);
            Z=Uqsl2BorelSmall();
            Z=Z.setBase(Z.bs0.get(N));
            Z.N=N;
            Z.M=M;
            if isfield(arg,'q')
                Z.q = arg.q; % Set q if provided
            elseif strcmp(arg.qtype,"complex")
                Z.q=exp(2*pi*1i*M/N);
            else
                Z.q=sym("z"+N)^M;
                assume(sym("z"+N)^N==1);
            end
            Z.setConst(style);
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
        function setConst(obj,style)
            % Set structure constants based on the chosen style
            % First, define with style=R: Δ(E)=1⊗E+E⊗K
            
            % Basis: K^i E^j (i,j) i,j=0,1,...,N-1
            % E^N=0, K^N=1, KE=qEK
            % Multiplication tensor M(i,j,k): e_i * e_j = sum_k M(i,j,k) * e_k
            N=obj.N;
            q=obj.q;
            M = zeros(N*[1 1 1 1 1 1],'like',q);
            for i1=0:N-1
                for i2=0:N-1
                    for j1=0:N-1
                        for j2=0:N-1
                            if j1+j2>=N
                                continue
                            end
                            M(i1+1,j1+1,i2+1,j2+1,mod(i1+i2,N)+1,j1+j2+1) = q^(-j1*i2);
                        end
                    end
                end
            end
            M=reshape(M,N^2*[1 1 1]);

            qN=qNumA(q);
            % Comultiplication tensor C(i,j,k): Δ(e_i) = sum_{j,k} C(i,j,k) * (e_j ⊗ e_k)
            % Δ(E)=1⊗E+E⊗K, Δ(K)=K⊗K,
            C = zeros(N*[1 1 1 1 1 1],'like',q);
            for i=0:N-1
                for j=0:N-1
                    for k=0:j
                        C(i+1,j+1,i+1,k+1,mod(i+k,N)+1,j-k+1) = ...
                            qN.nchoosek(j,k)*q^(-k*(j-k));
                    end
                end
            end
            C=reshape(C,N^2*[1 1 1]);

            % Counit vector: epsilon_i = ε(e_i)
            % ε(E)=0, ε(K)=1
            epsilon = [ones(N,1); zeros(N^2-N,1,'like',q)];

            % Counit vector: η(1) = η^i(e_i)
            eta = zeros(N^2,1,'like',q);
            eta(1) = 1;
            % Antipode matrix: S(j,i) = coefficient of e_j in S(e_i)
            % S(E)=-EK^-1, S(K)=K^-1
            S = zeros(N*[1 1 1 1],'like',q);
            for i=0:N-1
                for j=0:N-1
                    S(mod(-i-j,N)+1,j+1,i+1,j+1)=(-1)^j*q^(i*j+j*(j+1)/2);
                end
            end
            S=reshape(S,N^2*[1 1]);
            if style=="L"
                C=permute(C,[1,3,2]);
                S=S^-1;
            end
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
            
            [Ir,Cr,Il,Cl]=deal(zeros(N,N)*q);
            % redefine structure constants
            if style=="R"
                Ir(2,N)=1;
                Cr(:,N)=q.^((0:N-1)-1);
                Il(1,N)=1;
                Cl(:,N)=1;
            elseif style=="L"
                Ir(1,N)=q;
                Cr(:,N)=q.^((0:N-1)-1);
                Il(2,N)=q^1;
                Cl(:,N)=q^-1;
            end
            Ir=reshape(Ir,[N^2 1]);
            Cr=reshape(Cr,[N^2 1]);
            % Il=reshape(Il,[N^2 1]);
            % Cl=reshape(Cl,[N^2 1]);
            Cl=(S^-1)*Cr;
            Il=S.'*Ir;
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