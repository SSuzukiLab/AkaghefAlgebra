classdef CR<handle&dynamicprops& matlab.mixin.SetGet
    %CR このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant)
        H CR=CR
        ProjectPath = fileparts(mfilename('fullpath'));
    end
    properties
        name char='default'
        displayRule=0
        cft='double'
        pft='double'
        iszero 
        vectD1removeZero=true
        Fp_default=5
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
    end
end

