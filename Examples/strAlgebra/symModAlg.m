classdef symModAlg<symp
    properties
        rule= TypeParam(@(x)nan)
        type char
        idxV
    end
    methods
        function obj=symModAlg(varargin)
            obj@symp(varargin{:})
            obj.ctype="S";
            obj.ptype="S";
            obj.base.ctype="S";
            obj.base.ptype="S";
        end
        function ret=mtimes(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            type=i1.type;
            if isempty(type)
                type=i2.type;
            end
            assert(~isempty(type))
        end
    end
end