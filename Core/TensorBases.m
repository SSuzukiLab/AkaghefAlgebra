classdef TensorBases<Bases
    % TensorBases - A class to handle tensor of bases 
    
    properties
        bases (1,:) Bases % Array of bases
    end
    
    methods
        % Constructor
        function obj = TensorBases(bases,name)
            % Constructor for TensorBases
            obj.bases = bases;
            obj.var_ = string(bases);
            if nargin > 1
                obj.name = name;
            else
                obj.name = join({bases.name}, "⊗");
            end
        end
        
        % Define other methods here
    end
end