function Cov=enlargeCov(Cov,incre)

[U S V]=svd(Cov);

S(1,1)=(sqrt(S(1,1))+incre)^2;
S(2,2)=(sqrt(S(2,2))+incre)^2;

Cov=U*S*V;