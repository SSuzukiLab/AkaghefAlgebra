classdef BaseConverter<handle
    %BASECONVERTER : Base conversion manager


    properties(Constant)
        H BaseConverter=BaseConverter();
    end
    properties
        enable (1,:) logical
        base0 (1,:) Bases
        base1 (1,:) Bases
        matrix (1,:) cell % cell array of conversion matrix from base0 to base1
    end

    methods(Access=private)
        function obj = BaseConverter()
        end
    end
    methods

        function  add(obj,base0,base1,matrix)
            %ADD add conversion rule to BaseConverter
            arguments
                obj
                base0 (1,:) Bases
                base1 (1,:) Bases
                matrix double
            end
            obj.enable(end+1)=true;
            obj.base0(end+1)=base0;
            obj.base1(end+1)=base1;
            obj.matrix{end+1}=matrix;
        end
        function setEn(obj,base,en)
            obj.enable(obj.base==base)=en;
        end

        function ret=getIdx(obj,base)
            [~,ret]=ismember(base,obj.base0);
            dis=find(~obj.enable);
            ret(ismember(ret,dis))=0;
        end
    end
end