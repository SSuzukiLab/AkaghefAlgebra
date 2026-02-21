classdef VectDrinfeldDouble<VectAlg
    % VECTDRINFELDDOUBLE Drinfeld double D(H)=(H^*)^cop⋈H of finite dimensional Hopf algebra H
    % 
    properties % D(H)=H1⋈H2
        rdim %  dimension of H1
        H0 % H1⊗H2
        H1 % 1st Hopf algebra
        H2 % 2nd Hopf algebra
    end
    methods(Static)
        function Z=getGenerator(H1,H2,name)
            % getGenerator H(H2)=H1⋈H2 bicrossed product of Hopf algebras H1 and H2
            arguments
                H1 (1,1) VectAlg
                H2 (1,1) VectAlg
                name (1,:) char=''
            end
            D=H1.dim;
            Z=VectDrinfeldDouble();
            Z.rdim=D;
            Z.cf=H1.Czeros([D^2,1]);
            Z.H0=H1|H2;
            Z.H1=H1;
            Z.H2=H2;
            Z=Z.setBase(TensorBases([H1.bs H2.bs],name));
            Z.bs.helperHD;
            Z.setConst;
        end

        function Z=getGenerator1(H1,name)
            % getGenerator1 H(H1)=H1⋈(H1^*)
            arguments
                H1 (1,1) VectAlg
                name (1,1) string
            end
            dualobj=DualAlg.getGenerator(H1);
            Z=VectDrinfeldDouble.getGenerator(H1,dualobj,name);
        end
        function Z=getGenerator2(H2,name)
            % getGenerator2 H(H2)=(H2^*)^cop#H2
            arguments
                H2 (1,1) VectAlg
                name (1,1) string
            end
            dualobj=DualAlg.getGenerator(H2);
            Z=VectDrinfeldDouble.getGenerator(dualobj,H2,name);
        end
    end
    methods
        function setConst(obj)
            % set constants of Drinfeld double D(H)=(H*^cop)⋈H
            % MD: multiplication tensor MD\indices{^1_2^3_{45}^6}
            % Δ: comultiplication tensor of H \Delta\indices{_1^{23}}
            % μ: multiplication tensor of H \mu\indices{_{12}^3}
            % c: coefficient
            % (f⋈a)*(g⋈b) = <g_3,a_1> <g_1,S^-1(a_3)fg_2⋈a_2b
            % c1\indices{_1^2}c2\indices{_3^4}(e^1⋈e_2)*(e^3⋈e_4)
            % =MD\indices{^1_2^3_{45}^6}c1\indices{_1^2}c2\indices{_3^4}e^5⋈e_6
            % =c1\indices{_1^2}c2\indices{_3^4}e^5⋈e_6*
            %   μ\indices{_{78}^3}μ\indices{_{9,10}^8}Δ\indices{_2^{10,11}}Δ\indices{_11^{12,13}}*
            %   (S^{-1})\indices{^7_{13}}μ\indices{_{12,4}^6}Δ\indices{_5^{1,9}}
            %
            isQuasi=isa(obj.H1,"VectQuasiHopfAlg")||isa(obj.H2,"VectQuasiHopfAlg");
            
            H2=obj.H2;
            D=H2.dim;
            M=H2.getSC('prod');
            C=H2.getSC('coprod');
            eta=H2.getSC('unit');
            ep=H2.getSC('counit');
            S_inv=H2.getSC('antipode_inv');
            if ~isQuasi
                MD2=calcTensorExpression('M{7,8,3}M{9,10,8}C{2,10,11}C{11,12,13}S_inv{7,13}M{12,4,6}C{5,1,9}',1:6);
                MD = reshape(MD2, [D^2, D^2, D^2]);
                etaD=reshape(ep*eta.',D^2);
            else
                error("not impl")
                % Psi=H2.getSC('associator');
                % MH2=calcTensorExpression( ...
                %     'M{10,7,1}C{5,10,11}M{11,8,12}M{12,13,3}M{9,14,15}M{15,4,6}C{2,13,14}Psi{7,8,9}',1:6);
                % MH = reshape(MH2, [D^2, D^2, D^2]);
                % etaH=reshape(ep*eta.',[D^2,1]);
            end
            SC=obj.spec.SC;
            SC{'prod'}=MD;
            SC{'unit'}=etaD;
            obj.spec.SC=SC;
        end

        function [G,W,Wi]=getGW(obj)
            D=obj.rdim;
            Z=obj.H2;
            S=Z.getSC('antipode');
            mu=Z.getSC('prod');
            Delta=Z.getSC('coprod');
            eta=Z.getSC('unit');
            ep=Z.getSC('counit');
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
                eta=obj.H2.getSC('unit');
                ret=obj.set_c(arg.cf*eta.');
            elseif isequal(obj.H2.bs,arg.bs)
                ep=obj.H2.getSC('counit');
                ret=obj.set_c(ep*arg.cf.');
            else
                error("fail to cast")
            end
        end
        function ret=rep(obj)
            % representation of Heisenberg double in matrix form

            A=reshape(obj.cf,obj.rdim*[1 1]);
            mu=obj.H2.getSC('prod');
            Delta=obj.H2.getSC('coprod');
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