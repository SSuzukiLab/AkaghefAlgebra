classdef TypeParam<handle
    %TypeParam templateパラメータ用インターフェース
    %   詳細説明をここに記述
    properties
        dict dictionary =dictionary
        createFcn
    end
    methods
        function obj = TypeParam(createFcn)
            %TypeParam ファクトリーメソッドを登録する
            obj.createFcn=createFcn;
            obj.dict=dictionary;
        end
        function ret = get(obj,key)
            % get キャッシュにアクセス　なければ生成
            try
                ret=obj.dict(key);
            catch
                if isempty(obj.createFcn)
                    error("invalid access:"+string(key))
                end
                obj.dict=obj.dict.insert(key,{obj.createFcn(key)});
                ret=obj.dict(key);
            end
            ret=ret{1};
        end
        function insert(obj,key,val)
            obj.dict=obj.dict.insert(key,{val});
        end
    end
end