function [kd,km] = logsimilarity(interaction,nd,nm,logweight1,logweight2,logweight3)
gamad = nd/(norm(interaction,'fro')^2);
C=interaction;
kd=zeros(nd,nd);
D=C*C';
for i=1:nd
    for j=i:nd
        kd(i,j)=exp(-(gamad*(logweight1*D(i,j)+logweight2*(D(i,i)+D(j,j)-2*D(i,j))+logweight3*(nd-(D(i,i)+D(j,j)-D(i,j))))));
    end
end
kd = kd+kd'-diag(diag(kd));
kd=1./(1+kd);

gamam = nm/(norm(interaction,'fro')^2);
km=zeros(nm,nm);
E=C'*C;
for i=1:nm
    for j=i:nm
       km(i,j)=exp(-(gamam*(logweight1*E(i,j)+logweight2*(E(i,i)+E(j,j)-2*E(i,j))+logweight3*(nd-(E(i,i)+E(j,j)-E(i,j))))));
    end
end
km=km+km'-diag(diag(km));
km=1./(1+km);
end
