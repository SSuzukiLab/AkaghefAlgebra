classdef HeisenbergDouble<VectAlg
    % HEISENBERGDOUBLE heisenberg double of Hopf algebra H1#H2
    properties
        rdim %  dimension of H1
        H1 % 1st Hopf algebra
        H2 % 2nd Hopf algebra
    end
    methods(Static)
        function Z=getGenerator(H1,H2,name)
            % getGenerator H(H2)=H1#H2 smashproduct
            arguments
                H1 (1,1) VectAlg
                H2 (1,1) VectAlg
                name (1,:) char=''
            end
            if isempty(name)
                name=sprintf("H(%s)",H1.bs.name);
            end
            D=H1.dim;
            Z=HeisenbergDouble();
            Z.cf=zeros([D D],class(H1.cf));
            Z.bs=TensorBases([H1.bs H2.bs],name);
            
            Z.rdim=D;
            Z.H1=H1;
            Z.H2=H2;
            Z.ZERO={Z};
            
            % Z.setConst;
        end

        function Z=getGenerator1(H1,name)
            % getGenerator1 H(H1)=H1#(H1^*)
            arguments
                H1 (1,1) VectAlg
                name (1,1) string
            end
            dualobj=DualAlg.getGenerator(H1);
            Z=HeisenbergDouble.getGenerator(H1,dualobj,name);
        end
        function Z=getGenerator2(H2,name)
            % getGenerator2 H(H2)=(H2^*)#H2
            arguments
                H2 (1,1) VectAlg
                name (1,1) string
            end
            dualobj=DualAlg.getGenerator(H2);
            Z=HeisenbergDouble.getGenerator(dualobj,H2,name);
        end
    end
    methods
        function setConst(obj)
            H2=obj.H2;
            D=H2.dim;
            M=H2.SC.get([H2.identifier '_μ']);
            C=H2.SC.get([H2.identifier '_Δ']);
            MH=zeros(D*ones(1,6),class(H1.cf));
            T=repmat({1:D},1,9);
            T=combinations(T{:});
            for k=1:height(T)
                i=T{k,:};
                MH(i(2),i(1),i(4),i(3),i(5),i(6))=MH(i(2),i(1),i(4),i(3),i(5),i(6))+...
                    C(i(5),i(1),i(7))*M(i(7),i(9),i(3))*C(i(2),i(9),i(8))*M(i(8),i(4),i(6));
            end
            obj.SC.insert([Z.identifier '_μ'],MH);
        end
        function setConst2(obj,MH)
            D=obj.rdim;
            MH2 = reshape(permute(MH, [2 1 4 3 5 6]), D^2, D^2, D^2);
            obj.SC.insert([obj.identifier '_μ'],MH2);
        end
        % function ret=mtimes(i1,i2)
        % 
        % end
        function ret=times(i1,i2)
            % product of Heisenberg double
            MH=i1.SC.get([i1.identifier '_μ']);
            D=i1.rdim;
            cf1=reshape(i1.cf,[D D]);
            cf2=reshape(i2.cf,[D D]);
            cf3=zeros([D D],class(cf1));
            T=repmat({1:D},1,6);
            T=combinations(T{:});
            for k=1:height(T)
                i=T{k,:};
                cf3(i(5),i(6))=cf3(i(5),i(6))+...
                    MH(i(2),i(1),i(4),i(3),i(5),i(6))*cf1(i(1),i(2))*cf2(i(3),i(4));
            end
            ret=i1;
            ret.cf=reshape(cf3,[],1);
            
        end

    end
    methods
        % function obj=HeisenbergDouble()

        % end

        function ret=casttype(obj,arg)
            if ~isa(arg,'VectAlg')

            else
                error not_impl
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
end