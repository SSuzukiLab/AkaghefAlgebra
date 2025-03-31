classdef NumericType
    %NUMERIC このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        zero
        priority
        class
    end
    
    enumeration
        D (0,0)
        S (sym(0),1)
        P (sym(0),1)
    end
    methods
        function obj=NumericType(zero,priority)
            obj.zero=zero;
            obj.priority=priority;
            obj.class=class(zero);
        end
        function ret=getType(obj)
            m=max([obj.priority]);
            obj1=obj([obj.priority]==m);
            ret=obj1(1);
            assert(all(obj1==ret),"優先順位が混乱しています:%s",join(unique(string(obj1)),","))
        end
        function ret=zeros(obj,varargin)
            ret=repmat(obj.zero,varargin{:});
        end
    end
    
    
end

