classdef(InferiorClasses=?sym) symCnmodalg<symp
    %SYMC2 このクラスの概要をここに記述
    %   詳細説明をここに記述
properties(Constant,Hidden)
        B=TypeParam(@(l)Bases(2*l,["x_m"+(l:-1:1) "x_"+(1:l)],"symC"+l+"modalg"))
        
    end
    properties
        l
        q=sym('q')
    end
    properties(Dependent)
        idxs
    end
    methods (Static)
        function v=getGenerator(l)
            v=rowfun(@(x)symCnmodalg(1,x),table(eye(2*l)));
            v=dictionary([-l:-1,1:l]',v.Var1);
        end
    end
    methods
        function ret=get.idxs(obj)
            ret=[-obj.l:-1 1:obj.l];
        end
        function obj = symCnmodalg(varargin)
            %SYMC2 このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj@symp(varargin{:})
            obj.l=obj.dim/2;
            mustBeInteger(obj.l)
            B=obj.B.get(obj.l);
            B.ctype="S";
            B.ptype="S";
            % B.ptype="S";
            obj.base=B;
        end
        
        function ret=normaltimes(i1,i2)
            ret=prodeval(i1|i2);
        end
        
        function ret=mtimes(i1,i2)
            ret=prodeval5(i1|i2);
            % ret=prodeval(tmp);
            % ret=prodeval4(symCnmodalg(i1|i2))+0;
            ret.base=i1.base;
        end
        function ret=times(i1,i2)
            ret=prodeval(i1|i2);
        end

        function ret=prodeval4(arg)
            persistent MC operC
            if isempty(MC)
                MC=TypeParam(@getM);
                operC=TypeParam(@getOper);
            end
            l=arg.l;
            if l==2
                depth=max(min(arg.pw(:,[3 6])'));
            elseif l==3
                depth(1)=max(min(arg.pw(:,[5 8])'));
                depth(2)=max(min(arg.pw(:,[4 9])'));
            end
            M=MC.get(l);
            M=M{1};
            oper=operC.get(l);
            oper=oper{1};
            evaluated1=arg.lfun(@fun2);
            if l==2
            evaluated2=act(qNumS(arg.q).exp(oper,depth,true),evaluated1);
            elseif l==3
                evaluated2=act(qNumS(arg.q).exp(oper(1),depth(1),true),evaluated1);
                evaluated2=act(qNumS(arg.q).exp(oper(2),depth(2),true),evaluated2);
                evaluated2=act(qNumS(arg.q).exp(oper(3),depth(2),true),evaluated2);
            end
            ret=prodeval(evaluated2);

            function [c,p]=fun2(p)
                p1=p(:,1:2*l);
                p2=p(:,2*l+1:4*l);
                c=arg.q^(p1*M*p2.');
            end
            function ret=getM(ll)
                Mat=tril(ones(2*ll),-1);
                Mat(1+(2*ll-1)*(1:ll))=2
                ret={Mat};
                % if arg==2
                %     M=tril(ones(4),-1);
                %     M([4 7])=2;
                % elseif arg==3
                % 
                % end
            end
            function ret=getOper(ll)
                q=arg.q;
                tau=q-q^-1;
                if ll==2
                    I=strEndV.makeV(4);
                    op=tau*I.make("xd",4,3,[0 -1 0 0],q)|I.make("xd",1,2,[0 0 -1 0],q);
                    
                elseif ll==3
                    I=strEndV.makeV(2*ll);
                    op(1)=tau*I.make("xd",6,5,[0 1 0 0 0 0],q)|...
                          I.make("xd",1,2,[0 0 0 0 1 0],q);
                    op(2)=tau*q*I.make("xd",6,4,[0 1 1 0 1 0],q)|...
                          I.make("xd",1,3,[0 1 0 1 1 0],q);
                    op(3)=tau*I.make("xd",5,4,[0 0 1 0 0 0],q)|...
                          I.make("xd",2,3,[0 0 0 1 0 0],q);
                end
                ret={op};
                
            end
        end
        function ret=prodeval5(arg)
            persistent MC operC
            if isempty(MC)
                MC=TypeParam(@getM);
                operC=TypeParam(@getOper);
            end
            l=arg.l;
            M=MC.get(l);
            % M=M{1};
            oper=operC.get(l);
            expr1=arg.lfun(@fun2);
            expr2=act(oper,expr1);
            ret=prodeval(expr2);

            function [c,p]=fun2(p)
                p1=p(:,1:2*l);
                p2=p(:,2*l+1:4*l);
                c=arg.q^(p1*M*p2.');
            end
            function ret=getM(ll)
                ret=tril(ones(2*ll),-1);
                ret(1+(2*ll-1)*(1:ll))=2
            end
            function ret=getOper(l)
                q=arg.q;
                tau=q-q^-1;
                ret=1;
                I=strEndV.makeV(2*l);
                for i=1:l-1
                    op=1;
                    for j=i+1:l
                        v=zeros(2,2*l);
                        for k=i:j-1
                            v(1:2,[l+1-k,l+k])=1;
                        end
                        v(1,l+i)=v(1,l+i)-1;%+
                        v(2,l+1-i)=v(2,l+1-i)-1;%-
                        op=op+(tau*q^(l-j))*...
                        (I.make("xd",l+j,l+i,v(1,:),q)|I.make("xd",l+1-j,l+1-i,v(2,:),q));
                    end
                    
                    ret=ret*op;
                end
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