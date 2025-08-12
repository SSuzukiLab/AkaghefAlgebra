classdef(InferiorClasses=?sym) Usl2Full<strAlg&UEAlg
    properties(Constant)
        PairClass='FSL2Full'
        B=Bases(3,["E" "F" "H"],"sl2")
    end
    methods(Static)
        function [O,E,F,H]=getGenerator
            O=Usl2Full();
            O.ctype="S";
            O=O.make(0,{[]});
            E=O.make(1,{1});
            F=O.make(1,{2});
            H=O.make(1,{3});
        end
    end

    methods

        function obj=make(obj,cf,pw,~)
            obj=obj.make@strAlg(cf,pw,Usl2Full.B);
            
        end
        function ret=getP(obj,p)
            
        end

        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel(1)=obj.make([1 -1 -2],{[3 1] [1 3] 1});
                S.rel(2)=obj.make([1 -1 2],{[3 2] [2 3] 1});
                S.rel(3)=obj.make([1 -1 1],{[2 1] [1 2] 1});
                S.comm=[nan;nan];
                S.inv=[nan;nan];
                S=obj.get2vRelation_(S);

                % disp("更新")
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        function ret=HPGen(i1,i2)
            ret=HPGen(i2,i1);
        end
        function ret=HP(i1,i2)
            % ret=HP@HopfAlg(i1,i2);
            ret=HP(i2,i1);
        end

        % 簡約化の無視
        % function arg=replace(arg,Ntimes)
        % 
        % end

        function ret=rep(arg)
            persistent mat 
            if isempty(mat)
                mat={[0,1;0,0],[0,0;1,0],[1,0;0,-1],};
            end
            C=arrayfun(@(c,p)c*fold(@mtimes,mat(p{1}),eye(2)),arg.cf,arg.pw,UniformOutput=false);
            ret=sum(cat(3,C{:}),3);
        end
    end

end