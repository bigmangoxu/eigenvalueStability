function [v,dv,w,dw]=blockCirculantMatrixEig(b,N,M)

% eigenvalues and eigenvectors of the full matrix
B=[b1,b2,b3;b3,b1,b2;b2,b3,b1];

B=zeros(M*N, M*N);

for i=1:M
    for j=1:M
        if(j>=i)
            jshift=j-i+1;
            B((i-1)*N+1:i*N,(j-1)*N+1:j*N)=b(1:N,1:N,jshift);
        else
            jshift=j-i+1+M;
            B((i-1)*N+1:i*N,(j-1)*N+1:j*N)=b(1:N,1:N,jshift);
    end
end

[v,d]=eig(B);
% full matrix eigenvalues
dv=diag(d);

% eigenvalues/vectors of the block circulant matrix
theta=2*pi/M;
for m=1:M
    r(m)=cos((m-1)*theta)+i*sin((m-1)*theta);
end

B1=b1+r1*b2+r1^2*b3;
B2=b1+r2*b2+r2^2*b3;
B3=b1+r3*b2+r3^2*b3;

for m=1:M
    bt(1:N,1:N,m)=0;
    for m1=1:M
        bt(1:N,1:N,m)=bt(1:N,1:N,m)+r(m)^(m1-1)*b(1:N,1:N,m1);
    end
end

[v1,d1]=eig(B1);
[v2,d2]=eig(B2);
[v3,d3]=eig(B3);
i3=eye(3);
z3=zeros(3);
mat1=[i3 z3 z3; z3 r1*i3 z3; z3 z3 r1^2*i3];
mat2=[i3 z3 z3; z3 r2*i3 z3; z3 z3 r2^2*i3];
mat3=[i3 z3 z3; z3 r3*i3 z3; z3 z3 r3^2*i3];
w1=mat1*[v1;v1;v1];
w2=mat2*[v2;v2;v2];
w3=mat3*[v3;v3;v3];
w=[w1,w2,w3];
% block circulant matrix eigenvalues
dw=[diag(d1);diag(d2);diag(d3)];


for m=1:9
    for n=1:9
        if(abs(dv(m)-dw(n))<0.00001)            
          disp([m,n,norm( v(:,m)/v(1,m)-w(:,n)/w(1,n) )])
        end
    end
end
    


    
