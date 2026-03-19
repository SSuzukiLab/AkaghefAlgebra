classdef(InferiorClasses=?sym) StrKPAlg<StrAlg&UEAlg
    % https://arxiv.org/pdf/2505.00645
    properties(Constant,Hidden)
        B=Bases(3,["x","y","z"],"KP8")
    end
    
    %% generation
    methods
        function obj=StrKPAlg()
        end
        function obj=make(obj,cf,pw)
            obj=obj.make@StrAlg(cf,pw,obj.B);
        end
    end
    methods(Static)
        function [O,X]=getGenerator(N)
            arguments
                N=2
            end
            O=StrKPAlg();
            O=O.make(0,{[]});
            X=arrayfun(@(i)O.make(1,{i}),1:3);
        end
    end
    methods
        
        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent RS
            if isempty(RS)
                RS=TypeParam(@createRel);
            end
            S=RS.get(2);
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
            function S_=createRel(~)
                O=StrKPAlg.getGenerator(2);
                S_=struct;
                % zx=yz, zy=xz, xy=yx, x^2=1, y^2=1, 2*z^2=1+x+y-xy
                S_.rel=[O.make([1,-1],{T('zx'),T('yz')}), ...
                    O.make([1,-1],{T('zy'),T('xz')}), ...
                    O.make([1,-1],{T('xy'),T('yx')}), ...
                    O.make([1,-1],{T('xx'),[]}), ...
                    O.make([1,-1],{T('yy'),[]}), ...
                    O.make([-2,1,1,1, -1],{T('zz'),[],T('x'),T('y'),T('xy')})];
                S_.comm=[1;2]; % x,y commute
                S_.inv=[nan;nan];
                S_=O.get2vRelation_(S_);
                
            end
        end

        %% Hopf
        % comult
        function ret = Delta(obj)
            I=obj.unit;
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p)
                if p==1||p==2
                    X=I.set_cp(1,{p});
                    ret=X|X; % x,y: grp-like
                elseif p==3
                    % % Δ(z)=1/2(1⊗1+x⊗1+1⊗y-x⊗y)(z⊗z)
                    X=@(i)I.set_cp(1,{T(i)});
                    ret=(1/2)*((I|I)+(X('x')|I)+(I|X('y'))-(X('x')|X('y')))*(X('z')|X('z')); 
                end
            end
        end

        % Counit
        function ret = counit(obj)
            ret=sum(obj.cf);
        end

        % Antipode
        function ret = antipode(obj)
            I=obj.unit;
            ret=obj.antialgfun(@fun,I);
            function ret=fun(p)
                ret=I.set_cp(1,{p});
            end
        end

    end
end
function b=MakeBase()
    dec2bin(0:7)=='1'
end
function ret=T(xyz)
    [~,ret]=ismember(xyz,'xyz');
end