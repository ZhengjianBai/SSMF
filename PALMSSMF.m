function [W,H,Ptime,Terror,itera]=PALMSSMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2)
tstart=tic;
Ptime=zeros(10000,1);
Terror=Ptime;

for j=1:r
[sW0itemp,sW0indx]=maxk(W0(:,j),Ks1);
   sW0=zeros(1,m);
    sW0(sW0indx)=SimplexProj(sW0itemp');
    W0(:,j)=sW0';
end
for i=1:n
[sH0itemp,sH0indx]=maxk(H0(:,i),Ks2);
sH0=zeros(1,r);
sH0(sH0indx)=SimplexProj(sH0itemp');
H0(:,i)=sH0';
end

V0=W0*H0;
TV=V0.^(1/2)-V.^(1/2);
for l=1:n
Terror(1)=Terror(1)+(sum((TV(:,l)).^2)/2)^(1/2)/n;
end


itera=2;
while itera<maxit

%---------------------gradient descent-------------
gradW=-(V-W0*H0)*H0';
 stepw=1/((norm(H0,'fro'))^2+1e-6);
 W=W0-stepw*gradW;
  for i=1:r
    [Witemp,Windx]=maxk(W(:,i),Ks1);
    TWzero=zeros(1,m);
    TWzero(Windx)=SimplexProj(Witemp');
    W(:,i)=TWzero;
  end
    gradH=-W'*(V-W*H0);
%     beta=1/(norm(W*W')+1e-6);
    steph=1/((norm(W,'fro'))^2+1e-6);
    H=H0-steph*gradH;
    for j=1:n
   [Hitemp,Hindx]=maxk(H(:,j),Ks2);
   THzero=zeros(1,r);
    THzero(Hindx)=SimplexProj(Hitemp');
    H(:,j)=THzero;
    end

V0=W*H;
TV=V0.^(1/2)-V.^(1/2);
for l=1:n
Terror(itera)=Terror(itera)+(sum((TV(:,l)).^2)/2)^(1/2)/n;
end
tend=toc(tstart);
Ptime(itera)=tend;
W0=W;
H0=H;
itera=itera+1;


end

end

