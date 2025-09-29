classdef(HandleCompatible) TensorSpaceSpec
    %SPACESPEC このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant)
        INIT =TensorSpaceSpec() %初期化されたオブジェクトを
        % space specをdelegateするクラスに渡そうと考えたか，初期化プロセスは起動時に
        % 実行されるため，ほかのクラスのプロパティの初期値として利用することは不可能？.
    end
    properties
        base (1,:) Bases
        kind {mustBeMember(kind,["tensor" "scalar" "vector" "matrix"])} = "scalar"
        data
        SC dictionary %structure constant 
    end
    methods(Static)
        function obj=tensorSpec(i1,i2)
            obj=TensorSpaceSpec([i1.base i2.base]);
            obj.kind="tensor";
            obj.data=[i1 i2];
        end
    end
    methods
        function obj = TensorSpaceSpec(varargin)
            %SPACESPEC このクラスのインスタンスを作成
            %   詳細説明をここに記述
            cellfun(@(c)assert(isa(c,'SpaceSpec')))
        end
        function ret=or(i1,i2)
            if i1==i2
                ret=i1;
            else
                ret=TensorSpaceSpec([i1.base i2.base]);
            end
        end
        function ret=eq(i1,i2)
            ret = isequal(i1.base, i2.base);
        end
        function disp(obj)
            name=join(["" obj.base.name]);
            id=mod(keyHash(obj),1000);
            fprintf("space spec:%s, id:%d\n",name,id)
        end
    end
end