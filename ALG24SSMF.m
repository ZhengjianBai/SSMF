function [W,H,Ptime,Terror,itera]=ALG24SSMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2)
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
W=W0;
H=H0;
while itera<maxit
    
    %---------------------gradient descent-------------
    
    %-------------------update W------------------------
    VW=V-W0*H0;
    
    % if norm(W0*H0
    for i=1:r
        VW=VW+W0(:,i)*H0(i,:);
        VWTh=VW*(H0(i,:))';
        tempnormh=(norm(H0(i,:),2))^2;
        if tempnormh>1e-14
            
            [Witemp,Windx]=maxk((VWTh)/tempnormh,Ks1);
            Witemp=SimplexProj(Witemp');
            barW=sparse(1,m);
            barW(Windx)=Witemp;
            f2=(norm(VW-barW'*H0(i,:),'fro'))^2/2;
            f1=(norm(VW-W0(:,i)*H0(i,:),'fro'))^2/2;
            if f1-f2>10^(-6)*(norm(W0(:,i)-barW', 'fro'))^2/2
                W(:,i)=barW';
            else
                stepw=1/((norm(H0(i,:),2))^2+10^(-6));
                gradWl=-(VW-W0(:,i)*H0(i,:))*(H0(i,:))';
                WHl=W0(:,i)-stepw*gradWl;
                [Witemp,Windx]=maxk(WHl',Ks1);
                Witemp=SimplexProj(Witemp');
                wanW=sparse(1,m);
                wanW(Windx)=Witemp;
                f2=(norm(VW-wanW'*H0(i,:),'fro'))^2/2;
                f1=(norm(VW-W0(:,i)*H0(i,:),'fro'))^2/2;
                W(:,i)=wanW';
                if f2>=f1
                    W(:,i)=W0(:,i);
                end
            end
        else
            W(:,i)=W0(:,i);
        end
        VW=VW-W(:,i)*H0(i,:);
    end    
    gradH=-W'*(V-W*H0);
    for j=1:n
        gradj=gradH(:,j);
        downHtemp=W*gradj;
        if (norm(gradj)>1e-14)&&(norm(downHtemp)>1e-14)
            downHtemp=W*gradj;
            steph=(gradj)'*gradj/(downHtemp'*downHtemp);
            HT=H0(:,j)-steph*gradj;
            [Hitemp,Hindx]=maxk(HT,Ks2);
            barH=zeros(1,r);
            barH(Hindx)=SimplexProj(Hitemp');
            HT=barH';
            f1=(norm(V(:,j)-W*H0(:,j),'fro'))^2/2-(norm(V(:,j)-W*HT,'fro'))^2/2;
            if (f1>1e-6*(norm(H0(:,j)-HT,'fro'))^2/2)
                H(:,j)=HT;
            else
                steph2=1/(norm(W'*W)+1e-6);
                H(:,j)=H0(:,j)-steph2*gradj;
                [Hitemp,Hindx]=maxk(H(:,j),Ks2);
                wanH=zeros(1,r);
                wanH(Hindx)=SimplexProj(Hitemp');
                H(:,j)=wanH';
                f1=(norm(V(:,j)-W*H0(:,j),'fro'))^2/2-(norm(V(:,j)-W*wanH','fro'))^2/2;
                if f1<=0
                    H(:,j)=H0(:,j);
                end
            end
        else
            H(:,j)=H0(:,j);
        end
    end
    H0=H;
    W0=W;
    V0=W*H;
    TV=V0.^(1/2)-V.^(1/2);
    for l=1:n
        Terror(itera)=Terror(itera)+(sum((TV(:,l)).^2)/2)^(1/2)/n;
    end
    
    tend=toc(tstart);
    Ptime(itera)=tend;
    
    itera=itera+1;
    
end
end

