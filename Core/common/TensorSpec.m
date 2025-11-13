classdef TensorSpec < ISpaceSpec
    % TensorSpec: value classとして実装
    methods
        function obj = TensorSpec(i1, i2)
            if nargin == 2
                obj.base = [i1.base i2.base];
                obj.kind = "tensor";
                obj.data = [i1 i2];
            end
        end
        function disp(obj)
            name = join(["" obj.base.name]);
            id = mod(keyHash(obj), 1000);
            fprintf("tensor spec:%s, id:%d\n", name, id)
        end
    end
end
