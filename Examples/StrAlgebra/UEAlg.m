classdef UEAlg<HopfAlg
    %% Hopf
    methods
        % comult
        function ret = Delta(obj)
            I=obj.unit;
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p)
                % primitive-like
                X=I.set_cp(1,{p});
                ret=(X|I)+(I|X);
            end
        end

        % Counit
        function ret = counit(obj)
            len=cellfun(@length,obj.pw);
            ret=sum(obj.cf(len==0));
        end

        % Antipode
        function obj = antipode(obj)
            len=cellfun(@length,obj.pw);
            obj.cf=(-1).^len.*obj.cf;
        end
    end
end