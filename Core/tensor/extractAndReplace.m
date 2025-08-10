
function [outStr, T] = extractAndReplace(exprStr)
    pattern = '(?<name>[a-zA-Z_]\w*)\{(?<i>[^,}]+),(?<j>[^}]+)\}';
    tokens = regexp(exprStr, pattern, 'names');

    outStr = exprStr;
    nameList = strings(0);
    iList = [];
    jList = [];
    exprList = strings(0);

    for k = 1:length(tokens)
        name = tokens(k).name;
        iExpr = tokens(k).i;
        jExpr = tokens(k).j;

        try
            iVal = evalin('caller', iExpr);
            jVal = evalin('caller', jExpr);
        catch
            error('Failed to evaluate indices: %s or %s', iExpr, jExpr);
        end

        exprToken = sprintf('%s{%s,%s}', name, iExpr, jExpr);
        outStr = strrep(outStr, exprToken, 'tentative');

        nameList(end+1) = name;
        iList(end+1) = iVal;
        jList(end+1) = jVal;
        exprList(end+1) = exprToken;
    end

    T = table(nameList', iList', jList', exprList', ...
        'VariableNames', {'Name','I','J','Expr'});
end