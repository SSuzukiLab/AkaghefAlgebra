classdef GroupRing
    properties
        base
        M
    end
    methods
        function obj=GroupRing(M)
            obj.verifyM(M)
            obj.M=M;
            
        end
        function ret=mtimes(obj,i1,i2)

        end
    end
    methods(Static)
        function verifyM(M)
            assert(isa(M,"double"))
            N=size(M,1);
            assert(size(M,2)==N)
            assert(all(M(:) >= 0)&&all(M(:)<N))
            assert(all(M(1,:)==(0:N-1))&&all(M(:,1)==(0:N-1).'), ...
                "unit condition error")
            [Ml,Mr]=deal(zeros(N,N,N));
            for i = 1:N
                Ml(:,:,i) = reshape(M(M+1,i),[N,N,1]);
                Mr(i,:,:) = reshape(M(i,M+1),[1,N,N]);
            end
            assert(all(Ml==Mr,"all"),"associativity condition error")
        end
    end
end