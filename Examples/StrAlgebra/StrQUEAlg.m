classdef(InferiorClasses=?sym) StrQUEAlg<StrAlg&HopfAlg
    properties(Constant,Hidden,Abstract)
        % % B
        % B            %=TypeParam(@(N)Bases(5*N,reshape(["E","F","K","Ki" "H"]+(1:N)',1,5*N),"Uqsl_"+(N+1)))
        % CM           %=TypeParam(@createDelta)
        % RVL          %=TypeParam(@createRVL)
        % DeltaStorage %=TypeParam(@createDelta)
        % RepStorage   %=TypeParam(@createRep)
        % RelStorage   %=TypeParam(@createRel)
        % Vtype        %='symAnmodlg'
    end
    properties
        CD CartanData
        
    end
    properties(Abstract,Dependent)
        dimV
    end
    properties
        rank
        q=sym('q')
    end
    
    %% generation
    methods
        % function obj=StrAnalg(rank)
        %     obj.rank=rank;
        %     obj.algbase=obj.B.get(rank);
        % end
        % function obj=unit(obj)
        %     obj=unit@StrAlg(obj);
        % end
        function obj=make(obj,cf,pw)
            obj=obj.make@StrAlg(sym(cf),pw,obj.algbase);
        end
    end
    methods
        % function [obj,E,F,K,Ki,H]=getGenerator(obj,rank)
        %     % Template for getGenerator at subclass
        %     obj.rank=rank;
        %     obj.CD=CartanData("X"+rank);
        %     % obj.q=sym('q');
        %     [obj,E,F,K,Ki,H]=getGenerator@StrQUEAlg(obj,rank,base);
        %     obj.setRelation();
        %     obj.setDelta("R");
        % end
        function [obj,E,F,K,Ki,H]=getGenerator(obj,rank,base)
            arguments
                obj
                rank {mustBeInteger}
                base (1,1) Bases
            end
            obj=obj.set_cp(sym(0),{[]});
            obj.ctype="S";
            obj.rank=rank;
            obj.base=base;
            obj.spec.base=base;
            obj=obj.make(0,{[]});
            C=num2cell(1:rank);
            C(2:6,:)=reshape(arrayfun(@(i){obj.make(1,{i})},1:5*rank),rank,5).';
            E=dictionary(C{[1 2],:});
            F=dictionary(C{[1 3],:});
            K=dictionary(C{[1 4],:});
            Ki=dictionary(C{[1 5],:});
            H=dictionary(C{[1 6],:});
        end
    end
    methods(Access=protected)
        function setDelta(obj,kind)
            arguments
                obj
                % kind of coproduct: K is put at right, middle, or left of returned tensor of generator "E"
                kind {mustBeMember(kind,["R","M","L"])}
            end
            rank=obj.rank;
            [O,E,F,K,Ki,H]=obj.getGenerator(rank,obj.base);
            I_=O.unit;
            switch kind
                case "R"
                    % K (and Ki) on the right: Δ(E)=E⊗K+1⊗E, Δ(F)=F⊗1+Ki⊗F
                    val=[arrayfun(@(n)(E(n)|K(n))+(I_|E(n)),1:rank),...
                        arrayfun(@(n)(F(n)|I_)+(Ki(n)|F(n)),1:rank)];
                case "L"
                    % K (and Ki) on the left: Δ(E)=K⊗E+E⊗1, Δ(F)=1⊗F+F⊗Ki
                    val=[arrayfun(@(n)(K(n)|E(n))+(E(n)|I_),1:rank),...
                        arrayfun(@(n)(I_|F(n))+(F(n)|Ki(n)),1:rank)];
                otherwise
                    error("Unsupported kind of coproduct");
            end
            val=[val,arrayfun(@(n)K(n)|K(n),1:rank),...
                arrayfun(@(n)Ki(n)|Ki(n),1:rank),...
                arrayfun(@(n)(H(n)|I_)+(I_|H(n)),1:rank),...
                I_];
            obj.spec.SC{"Delta"}=val;
        end
    end
    methods
        function ret=convertBaseString(obj,pw,bs)
            ret=convertBaseString@StrAlg(obj,pw,bs);
        end
        %% relation
        function setRelation(obj)
            % create Relation of Quantum Groups from Cartan data
            % should have set q, CartanData, rank
            % =========================================================================
            % QUANTUM GROUP U_q(g) DEFINING RELATIONS (Drinfeld-Jimbo Algebra)
            % These relations define the generators K_i, E_i, F_i for the quantum group
            % associated with the Cartan matrix A = (a_ij) of a Lie algebra g.
            %
            % Notations:
            %   q_i = q^(d_i), where d_i = (\alpha_i, \alpha_i) / 2 is the symmetrizer.
            %   [n]_q! is the q-factorial, and [n choose k]_q is the quantum binomial coefficient.
            %   a_ij is the entry of the Cartan matrix A.
            % -------------------------------------------------------------------------
            % 1. Commutation of Cartan Generators (K):
            %   K_i * K_j = K_j * K_i
            % -------------------------------------------------------------------------
            % 2. Action of Cartan on Step Operators (E, F):
            %   K_i * E_j = q^((α_i, α_j)) * E_j * K_i
            %   K_i * F_j = q^(-(α_i, α_j)) * F_j * K_i
            %
            %   NOTE: The exponent is the inner product (α_i, α_j).
            % -------------------------------------------------------------------------
            % 3. Commutation of Opposite Step Operators:
            %   [E_i, F_j] = delta_ij * (K_i - K_i^(-1)) / (q_i - q_i^(-1))
            % -------------------------------------------------------------------------
            % 4. Quantum Serre Relations (i != j):
            %
            %   (E-operators):
            %   sum_{k=0}^{1-a_ij} (-1)^k * [1-a_ij choose k]_{q_i} * E_i^k * E_j * E_i^(1-a_ij-k) = 0
            %
            %   (F-operators):
            %   sum_{k=0}^{1-a_ij} (-1)^k * [1-a_ij choose k]_{q_i} * F_i^k * F_j * F_i^(1-a_ij-k) = 0
            %
            % =========================================================================
            r=obj.rank;
            R=dictionary(["E","F","K","Ki","H"],r*(0:4));
            q=obj.q;
            S=struct;
            CD=obj.CD;
            a=CD.SimpleRoots;
            CM=CD.CartanMatrix;
            d=CD.Symmetrizer;
            
            % [K_i,K_j]=0 Commutation of Cartan Generators (K):
            T=combinations(1:r,1:r);
            comm_K=[T.Var1+R("K"), T.Var2+R("K");
                T.Var1+R("Ki"),T.Var2+R("Ki");
                T.Var1+R("K"), T.Var2+R("Ki");
                T.Var1+R("Ki"),T.Var2+R("K");].';
            
            %  Action of Cartan on Step Operators (E, F):
            comm_KE=[];
            comm_KF=[];
            rel_KE=obj.empty;
            rel_KF=obj.empty;
            for i=1:r
                for j=1:r
                    if CM(i,j)==0
                        comm_KE([1,2],end+[1,2])=[i;j]+[R("K"),R("E");R("Ki"),R("E")].';
                        comm_KF([1,2],end+[1,2])=[i;j]+[R("K"),R("F");R("Ki"),R("F")].';
                    else
                        rel_KE(end+1)=obj.make([1 -q^(+dot(a(i),a(j)))],{[i+R("K"), j+R("E")] ...
                            ,[j+R("E"), i+R("K")]});
                        rel_KE(end+1)=obj.make([1 -q^(-dot(a(i),a(j)))],{[i+R("Ki"),j+R("E")] ...
                            ,[j+R("E"), i+R("Ki")]});
                        rel_KF(end+1)=obj.make([1 -q^(-dot(a(i),a(j)))],{[i+R("K"), j+R("F")] ...
                            ,[j+R("F"), i+R("K")]});
                        rel_KF(end+1)=obj.make([1 -q^(+dot(a(i),a(j)))],{[i+R("Ki"),j+R("F")] ...
                            ,[j+R("F"), i+R("Ki")]});
                    end
                end
            end
            
            % [E_i, F_j] = delta_ij * (K_i - K_i^(-1)) / (q_i - q_i^(-1))
            rel_EF=arrayfun(@(i)obj.make([1 -1 -q^(d(i))*[1, -1]/(q^(2*d(i))-1)], ...
                {[i+R("E") i+R("F")] [i+R("F") i+R("E")] i+R("K") i+R("Ki")}),1:r);
            T2=T;
            T2(T2.Var1==T2.Var2,:)=[];
            comm_EF=[T2.Var1+R("E"), T2.Var2+R("F")];
            
            % Quantum Serre Relations (i != j):
            rel_SerreE=obj.empty;
            rel_SerreF=obj.empty;
            for i=1:r
                for j=1:r
                    if i==j, continue; end
                    kMax=1-CM(i,j);
                    coeffs=arrayfun(@(k)(-1)^k*qNumS.binom(q^d(i),kMax,i),0:kMax);
                    
                    pw_=repmat({(i+R("E"))*ones(1,kMax+1)},1,kMax+1);
                    for k=0:kMax
                        pw_{k+1}(k+1)=j+R("E");
                    end
                    rel_SerreE(end+1)=obj.make(coeffs,pw_);
                    
                    pw_=repmat({(i+R("F"))*ones(1,kMax+1)},1,kMax+1);
                    for k=0:kMax
                        pw_{k+1}(k+1)=j+R("F");
                    end
                    rel_SerreF(end+1)=obj.make(coeffs,pw_);
                end
            end
            
            S.rel =[rel_KE,rel_KF,rel_EF,rel_SerreE,rel_SerreF];
            S.comm=[comm_K,comm_KE,comm_KF,comm_EF];
            S.inv =[R("K");R("Ki")]+(1:r);
            S=obj.get2vRelation_(S);
            obj.spec.SC{"Relation"}=S;
        end
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            S=obj.spec.SC{"Relation"};
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        %% rep
        function ret = repMono(obj)
            arr=obj.RepStorage.get(obj.Ntype);
            I=arr(end);
            converted=obj.algfun(@fun,I);
            ret=I.set_cp(converted.cf,converted.pw,converted.bs);
            ret.dimV=obj.rank+1;
            function ret=fun(p,b)
                % assert(b==StrStrUqsl2.B)
                ret=arr(p);
            end
        end
        function v=getModGenerator(obj)
            N=obj(1).rank;
            v=eval([obj.Vtype '.getGenerator(N)']);
            % v0=repmat(v0,[1 N]);
            % v=cellfun(@(p)v0.set_cp(1,p),mat2cell(eye(N),ones(1,N),N));
        end
        %% hopf
        function ret = Delta(obj)
            arr=obj.spec.SC{"Delta"};
            I=arr(end);
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p,b)
                % assert(b==StrStrUqsl2.B)
                ret=arr(p);
            end
        end
        
        % Counit
        function ret = counit(obj)
            ret=obj.lfun(@fun);
            ret=sum(ret.cf);
            function [c,p,b]=fun(p,b)
                c=all(p{1}==0);
                p={[]};
                b={[]};
            end
        end
        
        % Antipode
        function ret = antipode(obj)
            % NG
            ret=obj.lfun(@fun);
            function [c,p,b]=fun(p,b)
                c=(-1)^length(p{1});
            end
        end
    end
    methods(Static)
        function ret=getName()
            ret=eval(mfilename);
        end
    end
end
function ret = createRVL(Ntype)
    % createRVL: Returns root vector length of representation
    ret = ones(1,Ntype);
end

function ret = createRep(Ntype)
    % createRep: Returns a weyl algebra representation based on rank
    q=obj.q;
    [O,x,qth]=StrWeylXQ.getGenerator(Ntype+1);
    I_=O.unit;
    NanVal=I_;
    NanVal.cf=nan;
    ret=[arrayfun(@(k)(x(k)/x(k+1))*(qth(k+1)-1/qth(k+1))*(q-q^-1)^-1,1:Ntype), ...
        arrayfun(@(k)(x(k+1)/x(k))*(qth(k)-1/qth(k))*(q-q^-1)^-1,1:Ntype), ...
        arrayfun(@(k)qth(k)/qth(k+1),1:Ntype), ...
        arrayfun(@(k)qth(k+1)/qth(k),1:Ntype), ...
        repmat(NanVal,1,Ntype), ...
        I_];
end

function ret = createRel(rank)
    % createRel: Returns relationships based on rank
    error("not implemented")
end
