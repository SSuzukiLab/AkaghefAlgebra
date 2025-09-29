
classdef(Abstract) ICompare
% ICOMPARE Interface for comparison operations
%
% OBJECTIVE:
%   Define a minimal interface for objects that support relational comparisons.
%   Concrete subclasses must implement ge (>=) and eq (==). Default implementations
%   are provided for some other relational operators (e.g., gt using ge and eq).
%
% USAGE:
%   - Subclass ICompare and implement the following abstract methods:
%       ret = ge(objA, objB)   % Return true where objA >= objB
%       ret = eq(objA, objB)   % Return true where objA == objB
%
% NOTES:
%   - Default implementations in the interface may implement operators using
%     ge and eq (for example: gt = ge & ~eq). Subclasses should ensure their
%     ge and eq implementations are efficient and vectorized where appropriate.
%   - Handle type mismatch and invalid comparisons by returning logical false
%     or throwing informative errors according to class design.
    methods(Abstract)
        ret=ge(arg1,arg2)
        ret=eq(arg1,arg2)
    end
    methods
        function ret = gt(arg1, arg2)
            ret = ge(arg1, arg2) & ~eq(arg1, arg2);
        end
        
        function ret = le(arg1, arg2)
            ret = ~gt(arg1, arg2);
        end
        
        function ret = lt(arg1, arg2)
            ret = ~ge(arg1, arg2);
        end
        
        function ret = ne(arg1, arg2)
            ret = ~eq(arg1, arg2);
        end
    end
end