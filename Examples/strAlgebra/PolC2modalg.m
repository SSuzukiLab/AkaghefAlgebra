classdef(InferiorClasses=?Pol) PolC2modalg<PolAlg
    %SYMC2 このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        l=2
        q=sym('q')
    end
    methods (Static)
        function v=getGenerator(~)
            v=rowfun(@(x)PolC2modalg(1,x),table(eye(4)));
            v=dictionary([-2:-1,1:2]',v.Var1);
        end
    end
    methods
        function obj = PolC2modalg(varargin)
            %SYMC2 このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj@PolAlg(varargin{:})
            assert(obj.dim==obj.l*2)
            B=Bases(4,"x_"+["m2","m1","1","2"],"symplecticSpace");
            B.ctype="S";
            % B.ptype="S";
            obj.base=B;
        end
        function ret=normaltimes(i1,i2)
            [i1,i2]=StrAlg.alignNum(i1,i2);
            ret=prodeval(i1|i2);
        end
        
        function ret=mtimes(i1,i2)
            % [i1,i2]=alignNum(i1,i2);
            ret=prodeval4(i1|i2)+0;
            ret.base=i1.base;
        end
        function ret=times(i1,i2)
            ret=prodeval2(i1|i2);
        end
        function ret=prodeval4(arg)
            persistent M oper
            if isempty(M)
                M=tril(ones(4),-1);
                M([4 7])=2
                I=StrEndV.makeV(4);
                q=arg.q;
                oper=I.make("xd",4,3,[0 -1 0 0],q)|I.make("xd",1,2,[0 0 -1 0],q);
                oper.cf=q-q^-1;
            end
            depth=max(min(arg.pw(:,[3 6])'));
            evaluated1=arg.lfun(@fun2);
            evaluated2=act(qNumS(arg.q).exp(oper,depth,true),evaluated1);
            ret=prodeval(evaluated2);

            function [c,p]=fun2(p)
                p1=p(:,1:4);
                p2=p(:,5:8);
                c=arg.q^(p1*M*p2.');
            end
        end
        function ret=prodeval3(arg)
            persistent M
            if isempty(M)
                M=tril(ones(4),-1);
                M([4 7])=2
            end
            q=arg.q;
            qN=qNumS(q);
            
            % arg2=arg.lfun(@fun);
            % ret=arg2.lfun(@fun2);
            arg2=arg.lfun(@fun2);
            ret=arg2.lfun(@fun);
            ret=prodeval(ret);
            function [c,p]=fun(p)
                % 6=4+2
                p0=p;
                N=min(p(3),p(6));
                v=zeros(1,8);
                v([3 6])=1;
                v([1 8])=-1;
                % v([4 5])=-1;
                p=p-(0:N)'*v;
                e1=arrayfun(@qN.n,p(:,[3 6]));
                if N>0
                    e1(2:end,:)=e1(1:end-1,:);
                end
                e1(1,:)=1;
                e2=arrayfun(@(e)qN.fac(e),(0:N)');
                e3=arrayfun(@(e)qN.T^e*q^(-e*(e-1)/2),(0:N)');
                e4=(p0(2)+p0(7));
                e5=q.^((0:N)'*e4);
                c=prod(cumprod(e1),2)./e2.*e3.*e5;

                
            end
            function [c,p]=fun2(p)
                p1=p(:,1:4);
                p2=p(:,5:8);
                c=q^(p1*M*p2.');
            end
        end
        function ret=prodeval2(arg)
            q=arg.q;
            while(1)
                arg=arg.lfun(@fun);
                if all(~arg.pw(:,5:8),"all")
                    break
                end
            end
            ret=~arg;
            function [c,p]=fun(p)
                idx=0;
                p1=p(1:4);
                p2=p(5:8);
                if p2(1)~=0
                    c=q^(p1(4)+sum(p1(2:4)));
                    idx=1;
                elseif p2(2)~=0
                    c=[q^(2*p1(3)+p1(4));
                       q^(sum(p1(2:4)))*(q^(2*p1(3))-1)];
                    p=[p1 p2]+[0 1 0 0,  0 -1 0 0;
                               1 0 -1 1, 0 -1 0 0];
                    return
                elseif p2(3)~=0
                    c=q^(p1(4));
                    idx=3;
                elseif p2(4)~=0
                    c=1;
                    idx=4;
                else
                    c=1;
                    return
                end
                    p1(idx)=p1(idx)+1;
                    p2(idx)=p2(idx)-1;
                    p=[p1,p2];
            end
        end
    end
end

function ret=getq()
    syms q
    ret=1;
end