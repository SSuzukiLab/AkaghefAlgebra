classdef VectHomAlg <VectAlg
    %VECTHOMALG : Homomorphism between vector spaces
    %   sparse: tensor representation of linear map
    %   type_in/type_out: input/output type (Bases)
    %  **sparse must be arranged in the order of type_in, type_out**

    properties 
        type_in (1,:) Bases
        type_out (1,:) Bases
        spec_in (1,1) SpaceSpec
        spec_out (1,1) SpaceSpec
    end
    methods(Static)
        function obj = get(tensor,type_in,type_out,spec_in,spec_out)
            %get(tensor,type_in,type_out) get instance of hom between vector spaces
            %   tensor: tensor representation of linear map
            %   type_in/type_out: input/output type (Bases)
            arguments
                tensor (1,1) SparseEx
                type_in (1,:) Bases
                type_out (1,:)       Bases
                spec_in (1,1) SpaceSpec
                spec_out (1,1) SpaceSpec
            end
            assert(tensor.rank==length(type_in)+length(type_out),"invalid tensor rank")
            assert(isequal(tensor.size,[type_in.dims,type_out.dims]),"invalid tensor dimension")
            obj=VectHomAlg();
            obj.sparse=tensor;
            obj.type_in=type_in;
            obj.type_out=type_out;
            obj.spec_in=spec_in;
            obj.spec_out=spec_out;
            obj.bs=[type_in,type_out];
        end
        function idObj=id(arg)
            if isa(arg,'VectAlg')
                type=arg.bs;
                spec=arg.spec;
            elseif isa(arg,'Bases')
                type=arg;
                spec=makeSpec(type);    
            end
            idObj=VectHomAlg.get(SparseEx(eye(type.dim)),type,type,spec,spec);
        end
    end
    methods
        function obj = VectHomAlg()
            % Constructor for VectHomAlg
        end
        function ret=or(i1,i2)
            %OR(i1,i2) direct sum of homomorphisms
            %   i1,i2: VectHomAlg instances
            arguments
                i1 (1,1) VectHomAlg
                i2 (1,1) VectHomAlg
            end
            ret=or@VectAlg(i1,i2);
            ret.type_in=[i1.type_in,i2.type_in];
            ret.type_out=[i1.type_out,i2.type_out];
            Nio=[length(i1.type_in),length(i2.type_in),length(i1.type_out),length(i2.type_out)];
            ret.sparse=premute(ret.sparse, ...
            [1:Nio(1),Nio(1)+Nio(3)+(1:Nio(2)),Nio(1)+Nio(2)+(1:Nio(3)),Nio(1)+Nio(2)+Nio(3)+(1:Nio(4))]);
            ret.spec_in=or(i1.spec_in,i2.spec_in);
            ret.spec_out=or(i1.spec_out,i2.spec_out);
        end

        function ret=eval(hom,vec)
            %EVAL(hom,vec) Evaluate the homomorphism on a vector
            %   hom: VectHomAlg instance
            %   vec: input tensor
            arguments
                hom (1,1) VectHomAlg
                vec (1,1) VectAlg
            end
            assert(isequal(hom.type_in,vec.bs),"invalid input type")
            InputStr=join(string(1:length(hom.type_in)),',');
            OutputIdx=length(hom.type_in)+(1:length(hom.type_out));
            OutputStr=join(string(OutputIdx),',');
            expr=sprintf('hom.sparse{%s}vec.sparse{%s}',InputStr+","+OutputStr,InputStr);
            s=calcTensorExpression(expr,OutputIdx);
            ret=vec;
            ret.sparse=s;
            ret.bs=hom.type_out;
            ret.spec=hom.spec_out;
        end
    end
end
