classdef FSL3Borel<StrAlg&HopfAlg
    %UFSL2BOREL このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant,Hidden)
        B=Bases(9,["m"+[11,12,13,22,23,33] "l"+[1,2,3]],"FSL3")
        DeltaStorage=TypeParam(@createDelta)
        PairClass='Usl3Borel'
    end
    %% generation
    methods(Static)

        function [O,M,L]=getGenerator()
            O=FSL3Borel();
            O=O.make(0,{[]});
            M=repmat(O,3,3);
            M(1,1)=O.getGen(1);
            M(1,2)=O.getGen(2);
            M(1,3)=O.getGen(3);
            M(2,2)=O.getGen(4);
            M(2,3)=O.getGen(5);
            M(3,3)=O.getGen(6);
            L(1)=O.getGen(7);
            L(2)=O.getGen(8);
            L(2)=O.getGen(9);
        end
    end
    methods

        function obj=make(obj,cf,pw,~)
            obj=obj.make@StrAlg(cf,pw,FSL3Borel.B);

        end
        function ret=getGen(O,n)
            ret=O.make(1,{n});
        end
        %% algebra
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel=FSL3Borel.empty;
                T=combinations(1:9,1:9);
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
        function ret=repl(obj,N)
            % det=1に対応する関係式操作
        end
        %% Hopf
        % comult
        function ret = Delta(obj)
            arr=obj.DeltaStorage.get(3);
            ret=obj.algfun(@(p,~)arr(p),arr(end));
        end

        % Counit
        function ret = counit(obj)
            arr=[1,0,0,1,0,1,0,0,0];
            I=obj.unit;
            X=obj.algfun(@(p,~)I.set_cp(arr(p),I.pw,I.bs),I);
            ret=X.cf;
        end

        % Antipode
        function obj = antipode(obj)
            % ng
            len=cellfun(@length,obj.pw);
            obj.cf=(-1).^len.*obj.cf;
        end

        %% pairing
        function ret=HPGen(i1,i2)
            persistent HPmat
            if isempty(HPmat)
                % ([m11, m12, m13, m22, m23, m33, l1,l2,l3],[E1, E12, E2, H1, H2]) pairng result
                % --- 1. 表現行列の定義 ---
                M={
                    [0 1 0; 0 0 0; 0 0 0];
                    [0 0 1; 0 0 0; 0 0 0];
                    [0 0 0; 0 0 1; 0 0 0];
                    [1 0 0; 0 -1 0; 0 0 0];
                    [0 0 0; 0 1 0; 0 0 -1];
                    };
                indices = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
                HPmat=zeros(9,5);
                for i=1:6
                    for j=1:5
                        HPmat(i,j)=M{j}(indices(i,1),indices(i,2));
                    end
                end
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
    assert(N==3)
    I=FSL3Borel.getGenerator;
    I=I.unit;

    % 3×3行列の成分を配列 fun から取得
    m11 = I.getGen(1);
    m12 = I.getGen(2);
    m13 = I.getGen(3);
    m22 = I.getGen(4);
    m23 = I.getGen(5);
    m33 = I.getGen(6);

    l1=I.getGen(7);
    l2=I.getGen(8);
    l3=I.getGen(9);

    % 3×3行列の成分ごとの余積を計算
    DeltaM11 = (m11|m11);
    DeltaM12 = (m11|m12) + (m12|m22);
    DeltaM13 = (m11|m13) + (m12|m23) + (m13|m33);
    DeltaM22 = (m22|m22);
    DeltaM23 = (m22|m23) + (m23|m33);
    DeltaM33 = (m33|m33);

    % 余積を配列に格納
    DeltaM = [DeltaM11, DeltaM12, DeltaM13, DeltaM22, DeltaM23, DeltaM33];
    DeltaL=arrayfun(@(li)(li|I)+(I|li),[l1 l2 l3]);
    ret= [DeltaM,DeltaL,I|I];
end
