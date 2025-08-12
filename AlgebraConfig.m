classdef AlgebraConfig < dynamicprops & matlab.mixin.SetGet
    %TOPOLOGYCONFIG Configuration manager class

    properties(Constant)
        % Singleton instance
        H AlgebraConfig = AlgebraConfig();
        % Path of this class file
        ProjectPath = fileparts(mfilename('fullpath'));
    end

    properties


    end

    methods
        function obj = AlgebraConfig()
            persistent flag
            if ~isempty(flag)
                error(['TopologyConfig is a singleton. ' ...
                    'Use TopologyConfig.H to access the instance.']);
            end
            flag=1;
        end
        function add(obj,prop,value)
            arguments
                obj
            end
            arguments(Repeating)
                prop
                value
            end
            % Add a new property to the object
            for i=1:length(prop)
                if ~isprop(obj, prop{i})
                    addprop(obj, prop{i});
                end
                obj.(prop{i}) = value{i}; % Set the value
            end
        end
    end
end
