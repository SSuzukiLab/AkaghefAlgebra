function ret=eqD(i1,i2,tol)
    % i1,i2の絶対誤差がtol以下かを返す
    ret=abs(i1-i2)<tol;
end