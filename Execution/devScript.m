M1=repmat(v(1)*0,4,4,4);
M2=repmat(v(1)*0,4,4,4);
index=[-2 -1 1 2];
for ii=1:4
    for jj=1:4
        for kk=1:4
            % (**)*
            M1(ii,jj,kk)=(v(index(ii))*v(index(jj)))*v(index(kk));
            % *(**)
            M2(ii,jj,kk)=v(index(ii))*(v(index(jj))*v(index(kk)));
        end
    end
end