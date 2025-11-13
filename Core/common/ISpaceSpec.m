classdef (Abstract,HandleCompatible) ISpaceSpec 
    properties
        base
        kind {mustBeMember(kind,["tensor" "scalar" "vector" "matrix"])} = "scalar"
        data
        SC %structure constant 
    end
    methods (Abstract)
        disp(obj)
    end
    methods
        function ret = or(i1, i2)
            % ISpaceSpec同士のor演算。派生クラスでオーバーライド推奨
            if isequal(i1, i2)
                ret = i1;
            else
                % デフォルトはSpaceSpecを返すが、派生クラスで適切にオーバーライドすること
                ret = SpaceSpec([i1.base i2.base]);
            end
        end
    end
end
