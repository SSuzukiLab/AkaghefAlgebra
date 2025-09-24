classdef DualAlg<VectAlg
    %DUALALG Class representing a dual Hopf algebra of VectAlg
    properties
        % Store coefficients as a map from basis elements to values
        dualobj
    end
    methods(Static)
        function obj=getGenerator(dualobj)
            % getGenerator Create a dual algebra generator
            obj=DualAlg();
            obj=setDualObj(obj,dualobj);
            obj.setConst()
        end
    end
    methods
        function ret=get.dualobj(obj)
            ret=obj.spec.data.dual.base.ZERO;
        end
        % function obj = DualAlg(basis, coeffs)
        %     obj.basis = basis;
        %     obj.coeffs = coeffs;
        % end
        function obj = setDualObj(obj,dualobj)
            % setDualObj Set the dual object
            dualobj=dualobj.zeros();
            bs=dualobj.bs;
            obj.cf=dualobj.cf;
            obj=obj.setBase(Bases(bs.dim,"φ_{"+string(bs)+"}","dual_"+bs.name));
            obj.spec.data=struct(dual=dualobj.spec);
            
            % obj.bs.ZERO=obj;
        end
        function setConst(obj)
            SCdual=dictionary;
            SC=obj.spec.data.dual.SC;
            SCdual{'prod'}=permute(SC{'coprod'},[2,3,1]);
            SCdual{'coprod'}=permute(SC{'prod'},[3,1,2]);
            SCdual{'unit'}=SC{'counit'};
            SCdual{'counit'}=SC{'unit'};
            SCdual{'antipode'}=SC{'antipode'}.';
            obj.spec.SC=SCdual;
            obj.setHP();
        end


        % % Convolution product as Hopf algebra multiplication
        % function ret = mtimes(i1, i2)
        %     % Convolution product as Hopf algebra multiplication
        %     [i1,i2]=alignNum(i1,i2);
        %     assert(isequal(i1.bs,i2.bs),'異なる空間での積エラー')
        %     z=i1.ZERO{1};
        %     C=z.SC.get([z.identifier '_Δ']);
        %     ret=i1;
        %     ret.cf(:)=0;
        %     for k1=1:i1.dim
        %         for k2=1:i2.dim
        %             for k3=1:i1.dim
        %                 ret.cf(k3)=ret.cf(k3)+C(k1,k2,k3)*i1.cf(k1)*i2.cf(k2);
        %             end
        %         end
        %     end
        % end

        % % Comultiplication
        % function mapOut = Delta(obj)
        %     mapOut = containers.Map('KeyType','char','ValueType','double');
        %     for i = 1:length(obj.basis)
        %         bi = obj.basis{i};
        %         val = obj.coeffs(bi);
        %         % diagonal-like comultiplication
        %         key = [bi '_' bi];
        %         mapOut(key) = val;
        %     end
        % end

        % % Counit
        % function val = epsilon(obj)
        %     val = 0;
        %     for i = 1:length(obj.basis)
        %         val = val + obj.coeffs(obj.basis{i});
        %     end
        % end

        % % Antipode (simple involution-like example)
        % function sObj = S(obj)
        %     sCoeffs = containers.Map('KeyType','char','ValueType','double');
        %     for i = 1:length(obj.basis)
        %         bi = obj.basis{i};
        %         sCoeffs(bi) = -obj.coeffs(bi);
        %     end
        %     sObj = DualAlg(obj.basis, sCoeffs);
        % end

        function setHP(obj,hp)
            % Set the structure constants for the Hopf algebra pairing
            % This is a placeholder; actual implementation will depend on the specific algebra
            if nargin==1
                hp=SparseEx(eye(obj.dim));
            end
            obj.spec.SC{'HP'}=hp;
        end
        function ret=HP(i1,i2)
            % Hopf algebra pairing
            % i1: StrAlg
            % i2: StrAlg
            % M_ij=<dualobj_i,obj_j>
            assert(isa(i1,"VectAlg")&&isa(i2,"VectAlg"))
            dualpair=[(isa(i1,"DualAlg")&&i1.dualobj.spec==i2.spec), ...
                (isa(i2,"DualAlg")&&i2.dualobj.spec==i1.spec)];
            if dualpair(1)
                M=i1.getSC('HP');
            elseif dualpair(2)
                M = i2.getSC('HP').';
            else
                error("must be dual pair")
            end
            s1=i1.sparse;
            s2=i2.sparse;
            ret=toMatrix(calcTensorExpression('s1{1}M{1,2}s2{2}',[]));
        end
    end
end