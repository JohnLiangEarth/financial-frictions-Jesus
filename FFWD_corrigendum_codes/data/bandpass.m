function z=bandpass(fL,fH,K,y)
%Galo Nuño Bank of Spain
%May 2009
T=length(y);
z=zeros(T-2*K,1);
a=zeros(2*K+1,1);

alpha=zeros(K+1,1);
alpha(1)=2*(fH-fL);

for j=2:K+1;
    alpha(j)=(sin(2*pi*(j-1)*fH)-sin(2*pi*(j-1)*fL))/pi/(j-1);
end
theta=-(sum(alpha(2:K+1))*2+alpha(1))/(2*K+1);
a(K+1:2*K+1)=alpha+theta;

for i=1:K
    a(K+1-i)=alpha(i+1)+theta;
end

for i=1:T-2*K
    z(i)=sum(a.*y(i:(i+2*K)));
end







