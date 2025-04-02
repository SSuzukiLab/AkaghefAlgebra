classdef DualAlg<VectAlg
    %DUALALG Class representing a dual Hopf algebra of VectAlg
    properties
        % Store coefficients as a map from basis elements to values
        dualobj
    end
    
    methods
        % function obj = DualAlg(basis, coeffs)
        %     obj.basis = basis;
        %     obj.coeffs = coeffs;
        % end
        function obj = setDualObj(obj,dualobj)
            % setDualObj Set the dual object    
            obj.dualobj = dualobj;
            obj=DualAlg;
            bs=dualobj.bs;
            obj.bs=Bases(bs.dim,"φ_{"+string(bs)+"}","dual_"+bs.name);
            obj.cf=dualobj.cf;
            Z=obj;
            Z.cf(:)=0;
            obj.ZERO={Z};        
        end


        % Convolution product as Hopf algebra multiplication
        function ret = mtimes(i1, i2)
            % Convolution product as Hopf algebra multiplication
            [i1,i2]=alignNum(i1,i2);
            assert(isequal(i1.bs,i2.bs),'異なる空間での積エラー')
            z=i1.ZERO{1};
            C=z.SC.get([z.identifier '_Δ']);
            ret=i1;
            ret.cf(:)=0;
            for k1=1:i1.dim
                for k2=1:i2.dim
                    for k3=1:i1.dim
                        ret.cf(k3)=ret.cf(k3)+C(k1,k2,k3)*i1.cf(k1)*i2.cf(k2);
                    end
                end
            end
        end
        
        % Comultiplication
        function mapOut = Delta(obj)
            mapOut = containers.Map('KeyType','char','ValueType','double');
            for i = 1:length(obj.basis)
                bi = obj.basis{i};
                val = obj.coeffs(bi);
                % diagonal-like comultiplication
                key = [bi '_' bi];
                mapOut(key) = val;
            end
        end
        
        % Counit
        function val = epsilon(obj)
            val = 0;
            for i = 1:length(obj.basis)
                val = val + obj.coeffs(obj.basis{i});
            end
        end
        
        % Antipode (simple involution-like example)
        function sObj = S(obj)
            sCoeffs = containers.Map('KeyType','char','ValueType','double');
            for i = 1:length(obj.basis)
                bi = obj.basis{i};
                sCoeffs(bi) = -obj.coeffs(bi);
            end
            sObj = DualAlg(obj.basis, sCoeffs);
        end

        function setHP(obj)
            % Set the structure constants for the Hopf algebra pairing
            % This is a placeholder; actual implementation will depend on the specific algebra
            hp=eye(obj.dim);
            obj.SC.insert([obj.identifier '_HP'],hp);
            end
        end


        function ret=HP(i1,i2)
            % Hopf algebra pairing
            % i1: strAlg
            % i2: strAlg
            % [i1,i2]=alignNum(i1,i2);
            assert(isequal(i1.bs,i2.bs),'異なる空間での積エラー')
            if isa(i1,'DualAlg')
                z=i1.ZERO{1};
            else
                z=i2.ZERO{1};
            end            
            M=z.SC.get([z.identifier '_HP']);
            ret=i1;
            ret.cf(:)=0;
            for k1=1:i1.dim
                for k2=1:i2.dim
                    for k3=1:i1.dim
                        ret.cf(k3)=ret.cf(k3)+M(k1,k2,k3)*i1.cf(k1)*i2.cf(k2);
                    end
                end
            end
        end
    end
end