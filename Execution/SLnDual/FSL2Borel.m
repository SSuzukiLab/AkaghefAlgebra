classdef FSL2Borel<strAlg&HopfAlg
    %UFSL2BOREL このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant,Hidden)
        B=Bases(5,["m"+[11,12,22] "l"+[1,2]],"FSL2")
        DeltaStorage=TypeParam(@createDelta)
        PairClass='Usl2Borel'
    end
    %% generation
    methods(Static)

        function [O,M,L]=getGenerator()
            O=FSL2Borel();
            O=O.make(0,{[]});
            M=repmat(O,2,2);
            M(1,1)=O.getGen(1);
            M(1,2)=O.getGen(2);
            M(2,2)=O.getGen(3);
            L(1)=O.getGen(4);
            L(2)=O.getGen(5);
        end
    end
    methods

        function obj=make(obj,cf,pw,~)
            obj=obj.make@strAlg(cf,pw,FSL2Borel.B);
            
        end
        function ret=getGen(O,n)
            ret=O.make(1,{n});
        end
        %% algebra
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel=FSL2Borel.empty;
                T=combinations(1:5,1:5);
                S.comm=[T{:,:}'];
                S.inv=[1;3];
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
            arr=[1,0,1,0,0];
            I=obj.unit;
            X=obj.algfun(@(p,~)I.set_cp(arr(p),I.pw),I);
            ret=X.cf;
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
                % ([a b d log(a) log(d)],[E H]) pairng result
                HPmat=[0,1;1,0;0,-1;0,1;0,-1];
            end
            mustBeA(i1,i2.PairClass)
            if isa(i2,'FSL2Borel')
                [i1,i2]=deal(i2,i1);
            end
            p1=i1.pw{1};
            p2=i2.pw{1};
            ret=HPmat(p1,p2);
        end
    end

end
function ret=createDelta(N)
    assert(N==2)
    I=FSL2Borel.getGenerator;
    I=I.unit;
    % 2×2行列の成分を配列 fun から取得
    m11 = I.getGen(1);
    m12 = I.getGen(2);
    m22 = I.getGen(3);

    l1=I.getGen(4);
    l2=I.getGen(5);

    % 2×2行列の成分ごとの余積を計算
    DeltaM11 = (m11|m11);
    DeltaM12 = (m11|m12) + (m12|m22);
    DeltaM22 = (m22|m22);
    DeltaL=arrayfun(@(li)(li|I)+(I|li),[l1 l2]);
    ret= [DeltaM11, DeltaM12, DeltaM22,DeltaL,I|I];
end
