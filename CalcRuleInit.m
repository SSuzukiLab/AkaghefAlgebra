ac0=AlgebraConfig.H;
%[text] Display mode for URD: 0=default, 1=short, 2=long
%[text] 
ac0.add( ...
    disp_StrAlg =2, ...
    disp_VectAlg =1, ...
    SP_elim_zero =1, ...
    disp_PolAlg=1, ...
    simplify_func_Sparse=@simplify_func_Sparse);

function arg=simplify_func_Sparse(arg)
    if isa(arg,'sym')
        arg=simplify(arg);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
