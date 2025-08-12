

function [V, map] = allocateVariablesFromT(T)
    % T: table with columns Name, I, J, Expr
    % V: cell array where V{global_i, j} corresponds to T row
    % map: table of Name, LocalI, J, GlobalI, Expr

    % 1. Unique variable names
    names = unique(T.Name);
    copyCount = zeros(size(names));
    for k = 1:numel(names)
        copyCount(k) = max(T.I(T.Name == names(k)));
    end

    % 2. Compute offsets
    offsets = cumsum([0; copyCount(1:end-1)]);
    nameOffset = containers.Map(names, num2cell(offsets));

    % 3. Determine overall size
    totalCopies = sum(copyCount);
    rank = max(T.J);
    V = cell(totalCopies, rank);

    % 4. Fill V and mapping table
    G = zeros(height(T), 1);
    for i = 1:height(T)
        localI = T.I(i);
        j = T.J(i);
        name = T.Name(i);
        globalI = nameOffset(name) + localI;
        G(i) = globalI;
        V{globalI, j} = sprintf('%s{%d,%d}', name, localI, j);
    end

    map = T;
    map.GlobalI = G;
end
