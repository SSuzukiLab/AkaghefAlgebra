classdef(InferiorClasses=?sym) FSL2Full<strAlg&HopfAlg
    %UFSL2FULL このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant,Hidden)
        B=Bases(6,["m"+[11,12,21,22] "d"+[1,2]],"FSL2")
        DeltaStorage=TypeParam(@createDelta)
        PairClass='Usl2Full'
    end
    %% generation
    methods(Static)

        function [O,M,D]=getGenerator()
            O=FSL2Full();
            O.ctype="S";
            O=O.make(0,{[]});
            M=repmat(O,2,2);
            M(1,1)=O.getGen(1);
            M(1,2)=O.getGen(2);
            M(2,1)=O.getGen(3);
            M(2,2)=O.getGen(4);
            D(1)=O.getGen(5);
            D(2)=O.getGen(6);
        end
    end
    methods

        function obj=make(obj,cf,pw,~)
            obj=obj.make@strAlg(cf,pw,FSL2Full.B);
            
        end
        function ret=getGen(O,n)
            ret=O.make(1,{n});
        end
        %% algebra
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel=FSL2Full.empty;
                T=combinations(1:6,1:6);
                S.comm=[T{:,:}'];
                S.inv=[nan;nan];
                S=obj.get2vRelation_(S);

                % disp("更新")
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        %% Hopf
        % comult
        function ret = Delta(obj)
            arr=obj.DeltaStorage.get(2);
            I=arr(end);
            ret=obj.algfun(@(p,~)arr(p),I);
        end

        % Counit
        function ret = counit(obj)
            % arr=[1,0,0,1,0,0];
            % I=obj.unit;
            % X=obj.algfun(@(p,~)I.set_cp(arr(p),I.pw,I.bs),I);
            % ret=sum(X.cf);
            idx=~cellfun(@(p)any(ismember(p,[2 3])),obj.pw);
            ret=sum(obj.cf(idx));
        end

        % Antipode
        function obj = antipode(obj)
            len=cellfun(@length,obj.pw);
            obj.cf=(-1).^len.*obj.cf;
        end

        %% pairing
        function ret=HPGen(i1,i2)
            persistent HPmat
            if isempty(HPmat)
                % ([a b c d a-1 d-1],[E F H]) pairng result
                HPmat=[0,0,1;1,0,0;0,1,0;0,0,-1;0,0,1;0,0,-1];
            end
            mustBeA(i1,i2.PairClass)
            if isa(i2,'FSL2Full')
                [i1,i2]=deal(i2,i1);
            end
            p1=i1.pw{1};
            p2=i2.pw{1};
            ret=HPmat(p1,p2);
        end
        function ret=HP(i1,i2)
            if i1.term==1&&isscalar(i1.pw{1})
                M=rep(i2).';
                ret=i1.cf*M(i1.pw{1});
            else
                ret=HP@HopfAlg(i1,i2);
            end
        end

        % 簡約化の無視
        % function arg=replace(arg,Ntimes)
        % 
        % end
    end

end
function ret=createDelta(N)
    assert(N==2)
    I=FSL2Full.getGenerator;
    I=I.unit;
    % 2×2行列の成分を配列 fun から取得
    m11 = I.getGen(1);
    m12 = I.getGen(2);
    m21 = I.getGen(3);
    m22 = I.getGen(4);

    d1=I.getGen(5);
    d2=I.getGen(6);

    % 2×2行列の成分ごとの余積を計算
    DeltaM11 = (m11 | m11) + (m12 | m21);
    DeltaM12 = (m11 | m12) + (m12 | m22);
    DeltaM21 = (m21 | m11) + (m22 | m21);
    DeltaM22 = (m21 | m12) + (m22 | m22);

    DeltaD=arrayfun(@(li)(li|I)+(I|li),[d1 d2]);
    ret= [DeltaM11, DeltaM12, DeltaM21, DeltaM22,DeltaD,I|I];
end
