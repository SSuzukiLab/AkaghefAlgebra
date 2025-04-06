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
            D=H1.dim;
            Z=HeisenbergDouble();
            Z.rdim=D;
            Z.cf=H1.Czeros([D D]);
            Z.H1=H1;
            Z.H2=H2;
            Z=Z.setBase(TensorBases([H1.bs H2.bs],name));
            Z.bs.helperHD;
            Z.setConst;
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
            % set constants of Heisenberg double
            % MH: multiplication tensor MH\indices{^1_2^3_{45}^6}
            % Δ: comultiplication tensor \Delta\indices{_1^{23}}
            % μ: multiplication tensor \mu\indices{_{12}^3}
            % c: coefficient
            % c\indices{_1^2}c\indices{_3^4}(e^1\#e_2)*(e^3\#e_4)
            % =MH\indices{^1_2^3_{45}^6}c\indices{_1^2}c\indices{_3^4}e^5\#e_6
            % =MH\indices{^1_2^3_{45}^6}c\indices{_1^2}c\indices{_3^4}*
            % Δ\indices{_5^{17}}Δ\indices{_2^{98}}μ\indices{_{79}^3}μ\indices{_{84}^6}e^5e_6
            % 
            H2=obj.H2;
            D=H2.dim;
            M=H2.getSC('_μ');
            C=H2.getSC('_Δ');
            eta=H2.getSC('_η');
            ep=H2.getSC('_ε');
            MH2=permute(tensorprod(tensorprod(tensorprod(C,C) ...
                ,M,[3 5],[1 2]),M,4,1),[2 3 4 5 1 6]);
            MH = reshape(MH2, D^2, D^2, D^2);
            etaH=reshape(ep*eta.',D^2,1);
            obj.SC.insert([obj.identifier '_μ'],MH);
            obj.SC.insert([obj.identifier '_η'],etaH);
        end

        function [G,W,Wi]=getGW(obj)
            D=obj.rdim;
            Z=obj.H2;
            S=Z.getSC('_S');
            mu=Z.getSC('_μ');
            Delta=Z.getSC('_Δ');
            eta=Z.getSC('_η');
            ep=Z.getSC('_ε');
            G=obj;
            G.cf=reshape(tensorprod(tensorprod(tensorprod(Delta,S^-1,3,1) ...
                ,S^2,3,1),mu,[2 3],[1 2]),[D^2 1]);
            W=obj|obj;
            Wi=obj|obj;
            W.cf=reshape(tensorprod(tensorprod(ep,eye(D),Num=1),eta),[D^2 D^2]);
            Wi.cf=reshape(tensorprod(tensorprod(ep,S,Num=1),eta),[D^2 D^2]);
        end
        function ret=castype(obj,arg)
            if isequal(obj.H1.bs,arg.bs)
                eta=obj.H2.getSC('_η');
                ret=obj.set_c(arg.cf*eta.');
            elseif isequal(obj.H2.bs,arg.bs)
                ep=obj.H2.getSC('_ε');
                ret=obj.set_c(ep*arg.cf.');
            else
                error("fail to cast")
            end
        end
        function ret=rep(obj)
            % representation of Heisenberg double in matrix form

            A=reshape(obj.cf,obj.rdim*[1 1]);
            mu=obj.H2.getSC('_μ');
            Delta=obj.H2.getSC('_Δ');
            ret=tensorprod(tensorprod(A,mu,2,2),Delta,[1 3],[2 1]);
        end
        function ret=act(obj,arg)
            %  act on H^*=H1 
            A=rep(obj);
            ret=A*arg;
        end
    end
    methods
        % function obj=HeisenbergDouble()

        % end

        % function ret=casttype(obj,arg)
        %     if ~isa(arg,'VectAlg')
        % 
        %     else
        %         error not_impl
        %     end
        % 
        %     function obj=casttype(obj,arg)
        %         if isa(arg,'double')&&isscalar(arg)
        %             obj=obj.set_c([arg;0;0;0]);
        %         else
        %             error not_impl
        %         end
        %     end
        % end
    end
end