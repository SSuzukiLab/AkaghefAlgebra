classdef IAdditive
    %IADDITIVE このクラスの概要をここに記述
    %   詳細説明をここに記述

    methods(Abstract)
        o=plus(i1,i2)
        o=uminus(i1)
    end
    methods
        function ret=minus(i1,i2)
            ret=i1+(-i2);
        end

        function i1=uplus(i1)
        end
        function ret=sum(arr,dir)
            arguments
                arr 
                dir string {mustBeMember(dir,["1" "2" "all"])}="1"
            end
            if numel(arr)==0
                ret=0;
                return
            end
            if dir=="1"
                ret=arr(1,:);
                for i=2:size(arr,1)
                    ret=arrayfun(@plus,ret,arr(i,:));
                end
            elseif dir=="2"
                ret=arr(:,1);
                for i=2:size(arr,2)
                    ret=arrayfun(@plus,ret,arr(:,i));
                end
            elseif dir=="all"
                ret=arr(1);
                for i=2:numel(arr)
                    ret=ret+arr(i);
                end               
            end
        end
    end
end

