('fa{1}*ev(g, Sn(fa{2}{3},-1)*?*fa{2}{1}) | fa{2}{2}*b') %[output:8efd55e3]
A=parse_mul_tensor('fa{1}*ev(gb{1}, ev(S^-1,fa{2}{3})*?*fa{2}{1}) | fa{2}{2}*gb{2}') %[output:097efb4b]
%%

function ast = parse_mul_tensor(s)
%PARSE_MUL_TENSOR  Build a binary AST using only '*' and '|' with precedence.
%   - '{}' are ignored for parsing (they don't affect depth).
%   - Parentheses '()' control depth; we never split inside them.
%   - Precedence: '*' > '|' .  Both associate to the left.
%   - Output node:
%       leaf: struct('type','leaf','value',<substring>)
%       op  : struct('type','op','oper','*' or '|','left',<node>,'right',<node>)
%
% Example:
%   s = 'fa{1}*ev(gb{1}, ev(S^-1,fa{2}{3})*?*fa{2}{1}) | fa{2}{2}*gb{2}';
%   ast = parse_mul_tensor(s);

    if ~(ischar(s) || isstring(s)), error('Input must be char or string.'); end
    s = char(s);
    s = strtrim(s);
    if isempty(s), error('Empty expression.'); end
    ast = parseExpr(s);
end

function node = parseExpr(s)
    s = strtrim(stripOuterParens(s));
    % 1) lowest precedence: split on top-level '|', rightmost for left-assoc
    k = findTopLevelOp(s, '|');
    if k > 0
        left  = strtrim(s(1:k-1));
        right = strtrim(s(k+1:end));
        assertNonEmpty(left,  s, k, '|');
        assertNonEmpty(right, s, k, '|');
        node = struct('type','op','oper','|', ...
                      'left', parseExpr(left), ...
                      'right',parseExpr(right));
        return;
    end
    % 2) next: split on top-level '*', rightmost for left-assoc
    k = findTopLevelOp(s, '*');
    if k > 0
        left  = strtrim(s(1:k-1));
        right = strtrim(s(k+1:end));
        assertNonEmpty(left,  s, k, '*');
        assertNonEmpty(right, s, k, '*');
        node = struct('type','op','oper','*', ...
                      'left', parseExpr(left), ...
                      'right',parseExpr(right));
        return;
    end
    % 3) leaf
    node = struct('type','leaf','value',s);
end

function k = findTopLevelOp(s, op)
% Return the index of the RIGHTMOST occurrence of op at depth==0, else 0.
    depth = 0; k = 0;
    for i = 1:numel(s)
        c = s(i);
        if c == '('
            depth = depth + 1;
        elseif c == ')'
            depth = depth - 1;
            if depth < 0, error('Mismatched ")": position %d in "%s"', i, s); end
        end
    end
    if depth ~= 0, error('Unbalanced parentheses in "%s"', s); end

    depth = 0;
    for i = numel(s):-1:1
        c = s(i);
        if c == ')'
            depth = depth + 1;
        elseif c == '('
            depth = depth - 1;
            if depth < 0, error('Mismatched "(" scanning "%s"', s); end
        elseif c == op && depth == 0
            k = i; return;
        end
    end
end

function out = stripOuterParens(s)
% Strip one pair of outer parentheses if they enclose the ENTIRE expression.
    s = strtrim(s);
    if numel(s) >= 2 && s(1) == '(' && s(end) == ')'
        depth = 0;
        for i = 1:numel(s)
            c = s(i);
            if c == '('
                depth = depth + 1;
            elseif c == ')'
                depth = depth - 1;
                if depth == 0 && i < numel(s)
                    out = s; return; % outer parens don't enclose all
                end
            end
        end
        if depth ~= 0, error('Unbalanced parentheses in "%s"', s); end
        out = stripOuterParens(s(2:end-1)); % recursively strip
    else
        out = s;
    end
end

function assertNonEmpty(part, s, k, op)
    if isempty(part)
        error('Empty operand around operator %s at position %d in "%s".', op, k, s);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:8efd55e3]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'fa{1}*ev(g, Sn(fa{2}{3},-1)*?*fa{2}{1}) | fa{2}{2}*b'"}}
%---
%[output:097efb4b]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"A","value":"     type: 'op'\n     oper: '|'\n     left: [1×1 struct]\n    right: [1×1 struct]"}}
%---
