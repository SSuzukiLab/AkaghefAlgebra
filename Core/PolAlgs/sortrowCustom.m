function [B, idx] = sortrowCustom(A)
% sortrowCustom: 行列Aを行ごとにソートする（sortrowsと同等）
% 高速な自前sortrows実装
try
    A=double(A);
catch
end
[n, m] = size(A);
if ~(isnumeric(A) && isreal(A))
    [B, idx] = sortrows(A);
    return;
end
idx = (1:n)';

if isnumeric(A) && isreal(A)
    % 数値型・実数の場合
    % 各行をユニークな値に変換してソート
    % 2進法で重み付け（桁落ち防止のためdoubleで処理）
    maxA = max(abs(A(:)));
    if maxA == 0
        scale = 1;
    else
        scale = 2^(ceil(log2(maxA+1)));
    end
    % 重みベクトル
    weights = scale.^(m-1:-1:0);
    keys = A * weights(:);
    [~, idx] = sort(keys);
    B = A(idx,:);
else
    % その他（cell, char, 複素数など）はMATLAB組み込みsortrowsを利用
    B = sortrows(A);
end
end
