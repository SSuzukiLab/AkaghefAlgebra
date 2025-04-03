classdef CR<handle& matlab.mixin.CustomCompactDisplayProvider ...
            &dynamicprops& matlab.mixin.SetGet
    %CR このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant)
        H CR=CR
    end
    properties
        name char='empty'
        displayRule=0
        cft='double'
        pft='double'
        iszero 
        vectD1removeZero=true
    end
    properties(Hidden)
        replacingDisplay=false
    end
    
    methods
        function obj=CR(name)
            if nargin==1
                obj.name=name;
            end
        end

        function  ret=string(obj)
            ret=reshape(string({obj.name}),size(obj));
        end

        function rep = compactRepresentationForColumn(obj,displayConfiguration,~)
            rep = fullDataRepresentation(obj,displayConfiguration);
        end
    end
end

