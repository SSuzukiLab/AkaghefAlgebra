function ret=eqR(i1,i2,tol)
    % i1,i2のi1に対する相対誤差がtol以下かを返す
    ret=abs(i1-i2)<i1*tol;
end