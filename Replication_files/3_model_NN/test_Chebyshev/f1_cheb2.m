function PHI = f1_cheb2(x1,x2,a1,b1,a2,b2,n)
z1 = 2*x1/(b1-a1) - (a1+b1)/(b1-a1);
z2 = 2*x2/(b2-a2) - (a2+b2)/(b2-a2);
m = length(x1);
PHI = zeros(m,(n+1)^2);
T1  = zeros(m,(n+1));
T2  = zeros(m,(n+1));
T1(:,1) = 1;
T1(:,2) = z1;
T2(:,1) = 1;
T2(:,2) = z2;

for i =3:n+1
    T1(:,i) = 2*z1.*T1(:,i-1) - T1(:,i-2);
    T2(:,i) = 2*z2.*T2(:,i-1) - T2(:,i-2);
end

for j =1:n+1
    for i=1:n+1
        PHI(:,i+(j-1)*(n+1)) = T1(:,i).*T2(:,j);
    end
    
end



    
    