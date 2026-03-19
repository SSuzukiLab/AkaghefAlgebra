classdef loopProp
    %UNTITLED2 このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        prop
    end

    methods
        function obj = loopProp()
            %UNTITLED2 このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj.prop=loopProp2();
            obj.prop.prop=obj;
        end
    end
end