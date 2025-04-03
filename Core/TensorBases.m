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
            vars=fliplr({bases.var_});
            Tvar=combinations(vars{:});
            obj.dim_=prod([bases.dim_]);
            obj.var_ =join(fliplr(Tvar{:,:})," ⊗ ",2);
            if nargin > 1
                obj.name = name;
            else
                obj.name = join({bases.name}, "⊗");
            end
        end
        
        % Define other methods here
    end
end