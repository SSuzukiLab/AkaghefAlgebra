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
    ret=zeros([maxDeg numel(expr)],'sym');
    Nvar=1:numel(var);
    if isempty(Nvar)
        ret=zeros([maxDeg+1 1]);
        return
    end
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
        cidx=arrayfun(@(x){x:-1:1},size(tmp.',Nvar));
        ret(cidx{:},ii)=tmp;
    end
    ret=permute(ret,[1+Nvar 1]);

end

