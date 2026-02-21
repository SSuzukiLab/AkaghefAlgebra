function ret = coeffsFixed(expr,maxDeg,var)
    %COEFFSFIXED 多項式exprをvarについてmaxDeg次までの係数配列を返す
    %   多項式でないと怒られる
    % varが複数変数，エラー処理は未対応
    %
    arguments
        expr (:,1)
        maxDeg
        var =symvar(expr)
    end
    inf_flag=isinf(maxDeg);
    maxDeg(inf_flag)=0;
    ret=zeros([maxDeg+1 numel(expr)],'sym');
    Nvars=numel(var);
    if isempty(Nvars)
        ret=zeros([maxDeg+1 1]);
        return
    end
    Ivars=1:Nvars;
    if Nvars==1, Ivars=2; end
        
    for ii=1:numel(expr)
        % try
        tmp=coeffs(expr(ii),var,"All");
        % catch ME
        %     % 多項式にsimplifiedされていないとエラーになるらしい
        %     if strcmp(ME.identifier,'symbolic:coeffs:NotAPolynomial')
        %         tmp=coeffs(simplify(expr(ii)),var,"All");
        %     end
        % end
        if isempty(tmp),continue; end
        cidx=arrayfun(@(x){x:-1:1},size(tmp,Ivars));
        ret(cidx{:},ii)=tmp;
    end
    if any((size(ret,1:Nvars)>maxDeg+1)&~inf_flag)
        if AlgebraConfig.H.coeffs_fidxed_warning
            warning("actual degree exceeds")
        end
        maxDeg(inf_flag)=size(ret,find(inf_flag))-1;
        cidx=arrayfun(@(x){1:x},maxDeg+1);
        ret=ret(cidx{:});
    end
    % ret=permute(ret,[2:Nvars+1 1]);

end

