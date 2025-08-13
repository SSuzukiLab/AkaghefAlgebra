classdef(InferiorClasses=?sym) StrUqsl2<StrAlg&HopfAlg
    properties(Constant,Hidden)
        B=Bases(4,["E" "F" "K" "Ki"],"Uqsl2")
        q
    end
    
    %% generation
    methods(Static)
        function obj=make(cf,pw,~)
            obj=StrUqsl2().make@StrAlg(cf,pw,StrUqsl2.B);
            obj.ctype="S";
        end
        function [O,E,F,K,Ki]=getGenerator(q)
            O=StrUqsl2();
            O=O.make(0,{[]});
            E=O.make(1,{1});
            F=O.make(1,{2});
            K=O.make(1,{3});
            Ki=O.make(1,{4});
        end
    end
    methods
        
        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                q=sym('q');
                % KE=q^2EK
                S.rel(1)=StrUqsl2.make([1 -q^2],{[3 1] [1 3]});
                S.rel(2)=StrUqsl2.make([1 -q^-2],{[4 1] [1 4]});
                % KF=q^-2FK
                S.rel(3)=StrUqsl2.make([1 -q^-2],{[3 2] [2 3]});
                S.rel(4)=StrUqsl2.make([1 -q^2],{[4 2] [2 4]});
                % KK^-1=1
                S.rel(5)=StrUqsl2.make([1 -1],{[4 3] []});
                S.rel(6)=StrUqsl2.make([1 -1],{[3 4] []});
                S.rel(7)=StrUqsl2.make([1 -1 -[1 -1]/(q-q^-1)],{[1 2] [2 1] 3 4});
                S.comm=[nan;nan];
                S.inv=[nan;nan];
                S=obj.get2vRelation_(S);
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        %% representation
        function ret = repMono(obj)
            persistent arr I
            if isempty(arr)
                [O,x,qth]=StrWeylXQ.getGenerator(2);
                I=O.unit;
                q=sym('q');
                arr=[(x(1)/x(2))*(qth(2)-1/qth(2))*(q-q^-1)^-1, ...
                     (x(2)/x(1))*(qth(1)-1/qth(1))*(q-q^-1)^-1, ...
                     qth(1)/qth(2), qth(2)/qth(1)];

                % % 
            end
            converted=obj.algfun(@fun,I);
            ret=I.set_cp(converted.cf,converted.pw,converted.bs);
            ret.dimV=2;
            function ret=fun(p,b)
                % assert(b==StrUqsl2.B)
                ret=arr(p);
            end
        end

        %% Hopf
        % comult
        function ret = Delta(obj)
            persistent dict I
            if isempty(dict)
                [O,E,F,K,Ki]=obj.getGenerator;
                I=O.unit;
                dict=dictionary( ...
                    1,(E|I)+(K|E), ...
                    2,(F|Ki)+(I|F), ...
                    3,K|K, ...
                    4,Ki|Ki);
                % 
            end
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p,b)
                % assert(b==StrUqsl2.B)
                ret=dict(p);
            end
        end

        % Counit
        function ret = counit(obj)
            ret=obj.lfun(@fun);
            ret=sum(ret.cf);
            function [c,p,b]=fun(p,b)
                c=all(p{1}==0);
                p={[]};
                b={[]};
            end
        end

        % Antipode
        function ret = antipode(obj)
            % NG
            ret=obj.lfun(@fun);
            function [c,p,b]=fun(p,b)
                c=(-1)^length(p{1});
            end
        end

    end
end
      