function V0=correctCov(V0,area)


[U S V]=svd(V0);    
l1=sqrt(S(1,1));    
l2=sqrt(S(2,2));

k=area/(pi*l1*l2); 
S=k*S;            

V0=U*S*V;         