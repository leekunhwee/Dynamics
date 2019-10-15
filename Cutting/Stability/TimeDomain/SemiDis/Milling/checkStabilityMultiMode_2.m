function eigMax = checkStabilityMultiMode_2(Cq,Kq,invMq,UU,Delta,a_xx,a_xy,a_yx,a_yy,m,dt)
%Two modes per axis:
%Cq,Kq->[4x4]
%Px,Py, mode matrix, [2x2]

L=cell(1,m+1);
R=cell(1,m+1);

for ii=1:m+1
    A=[a_xx(ii),a_xy(ii);
        a_yx(ii),a_yy(ii);];
    UAU=UU'* A*UU;
    L{ii}=[
        zeros(4,4),eye(4);
        Delta*invMq*UAU-Kq,-Cq;
        ];
    R{ii}=[
        zeros(4,8);
        -Delta*invMq*UAU,zeros(4,4)
        ];
    
end

invL=cell(1,m+1);
for ii=1:m+1
    invL{ii}=inv(L{ii});
end
B=cell(1,m+1);
for ii=1:m+1
    B{ii}=[
        expm(L{ii}*dt), zeros(8,8*m-16) , 1/2*(expm(L{ii}*dt)-eye(8))*invL{ii}*R{ii},1/2*(expm(L{ii}*dt)-eye(8))*invL{ii}*R{ii};
        eye(8*m),zeros(8*m,8)
        ];
end

PHI=B{1};
for ii=2:m+1
    PHI=B{ii}*PHI;
end

eigValue=eig(PHI);
eigMax=max(abs(eigValue));