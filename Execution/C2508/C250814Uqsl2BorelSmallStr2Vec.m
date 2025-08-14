%[text] # C250814$u\_q(sl\_2^+)\n$　Str→Vec
% Define
syms q
[O,E,F,K,Ki]=StrUqsl2().getGenerator(q);
assume(q>0)
%%
N=5;
Mu0=zeros(N*ones(1,6),'sym');
Mu1=zeros(N*ones(1,7));
T=combinations(1:N,1:N,1:N,1:N);
for i=1:height(T)
    j=T{i,:};
    tmp=E^(j(1)-1)*K^(j(2)-1)*E^(j(3)-1)*K^(j(4)-1);
    for k=1:tmp.term
        j(5)=sum(tmp.pw{1}==1)+1;
        if j(5)>N, continue; end
        j(6)=mod(sum(tmp.pw{1}==3),N)+1;
        [c,t]=coeffsL(tmp.cf,q,"int");
        t=mod(t/2,N)+1;
        for l=1:length(t)
            Mu1(j(1),j(2),j(3),j(4),j(5),j(6),t(l))=c(l);
        end
        Mu0(j(1),j(2),j(3),j(4),j(5),j(6))=tmp.cf;
    end
end
Mu1=double(reshape(Mu1,[N^2,N^2,N^2,N]));
Mu0=reshape(Mu0,[N^2,N^2,N^2]);
%%
Mu3=double(subs(Mu0,q,exp(2*pi*1i/N)));
sum(~eqD(Mu3,0,1e-100),"all") %[output:1771ebf2]
%%
size(TP(Mu3,Mu3,3,1)) %[output:3394ad55]
%%

Mu2=zeros(N,N,N);
for i=1:N
    for j=1:N
        Mu2(i,j,mod(i+j-2,N)+1)=1;
    end
end
%%
size(Mu1) %[output:44a02b22]
size(TP(Mu1,Mu1,3,1)) %[output:3ad06110]
A1=TP(TP(Mu1,Mu1,3,1),Mu2,[3,6],[1,2]);
size(A1) %[output:5923c5c8]
size(TP(Mu1,Mu1,2,3)) %[output:86f8a1a5]
A2=permute(TP(TP(Mu1,Mu1,2,3),Mu2,[3,6],[1,2]),[1,3,4,2,5]);
size(A2) %[output:0f1363d0]
%%
all(A1==A2,"all") %[output:604dbb99]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:1771ebf2]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"375"}}
%---
%[output:3394ad55]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"ans","rows":1,"type":"double","value":[["9","9","9","9"]]}}
%---
%[output:44a02b22]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"ans","rows":1,"type":"double","value":[["25","25","25","5"]]}}
%---
%[output:3ad06110]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["25","25","5","25","25","5"]]}}
%---
%[output:5923c5c8]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"ans","rows":1,"type":"double","value":[["25","25","25","25","5"]]}}
%---
%[output:86f8a1a5]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["25","25","5","25","25","5"]]}}
%---
%[output:0f1363d0]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"ans","rows":1,"type":"double","value":[["25","25","25","25","5"]]}}
%---
%[output:604dbb99]
%   data: {"dataType":"textualVariable","outputData":{"header":"logical","name":"ans","value":"   1"}}
%---
