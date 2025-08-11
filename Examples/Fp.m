classdef Fp<double
    %FP このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        p
    end

    
    methods
        function obj = Fp(z,p)
            %FP このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj@double(z);
            if nargin==1
                p=CR.H.Fp_default;
            end
            obj.p=p;
        end
        function [p,arg1,arg2]=align(arg1,arg2)
            if isa(arg1,"Fp")
                p=arg1.p;
                if isa(arg2,"Fp")
                    assert(p==arg2.p,"inconsistent mod p")
                end
            else
                p=arg2.p;
                if isa(arg1,"Fp")
                    assert(p==arg1.p,"inconsistent mod p")
                end
            end
        end
        function ret=plus(arg1,arg2)
            p=align(arg1,arg2);
            val=plus@double(arg1,arg2);
            ret=Fp(mod(val,p),p);
        end
        function ret=s(arg)
            ret=Fp(mod(double(arg),arg.p),arg.p);
        end
        % function ret=minus(arg)
        function outputArg = add(obj,inputArg)
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
            outputArg = obj.Property1 + inputArg;
        end
    end
end

