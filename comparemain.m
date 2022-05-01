clear,clc;
% %-----------data 1----------------------------------------------------
% load('V.mat');
% 
% V=im2double(V);
% V=V';
% [m,n]=size(V);
% for i=1:n
%     V(:,i)=V(:,i)/sum(V(:,i));
% end
% 
% %-----------------------------------------------
% 
% r=30;
% 
% %--------------------------------------------------
% 
% maxit=350;
%  Ks1=4500;
%  Ks2=5;
%  
%  %-------------The initial point is determined in algorithm ARTMsparse-----
%  %---------------ARTMsparse----------------------------------------------
%  [AW,AH,Alls,Atime,Aerror,Aiter] = ARTMsparse(V,r,maxit);
%  save('AW.mat','AW');
%  save('AH.mat','AH');
%  save('Atime.mat','Atime');
%  save('Aerror.mat','Aerror');
%  
%  %-------------------sparse---------------------------
%  sparw4=sum(sum(AW~=0));
%  sparh4=sum(sum(AH~=0));
% 
%  
% load('W0.mat');
% load('H0.mat');
%  
% %-----------------------Alg. 2.2 for SSNMF-------------------------------
% [TW,TH,Ttime,Terror,Titera]=PALMSSNMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2);
%  save('TW.mat','TW');
%  save('TH.mat','TH');
%  save('Ttime.mat','Ttime');
%  save('Terror.mat','Terror');
%  sparw1=sum(sum(TW~=0));
%  sparh1=sum(sum(TH~=0));
% 
% 
% 
% %--------------------------Alg. 2.4 for SSNMF----------------------------
% [PW,PH,Ptime,Perror,Pitera]=ALG24SSNMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2);
% 
% 
%  save('PW.mat','PW');
%  save('PH.mat','PH');
%  save('Ptime.mat','Ptime');
%  save('Perror.mat','Perror');
% 
% sparw2=sum(sum(PW~=0));
%  sparh2=sum(sum(PH~=0));
% 
% 
%  
% figure(1) 
% semilogy((1:maxit),Terror(1:maxit),'r-',(1:maxit),Perror(1:maxit),'b-',(1:maxit),Aerror(1:maxit),'m-d','MarkerSize',6,'LineWidth',1.5)
% legend('Alg. 2.2','Alg. 2.4','ARTMsparse')
% xlabel('Total number of iterations')
% ylabel('$${H(V,W^kH^k)}$$','Fontsize',14,'FontWeight','bold','interpreter','latex') 
% 

%-------------------------------data2(sythen)---------------------------------
m=400;
n=200;
r=50;
s1=80;
s2=10;
W=rand(m,r);
for j=1:r
[sWitemp,sWindx]=maxk(W(:,j),s1);
   sW=zeros(1,m);
    sW(sWindx)=SimplexProj(sWitemp');
    W(:,j)=sW';
end
% 
H=rand(r,n);
for i=1:n
[sHitemp,sHindx]=maxk(H(:,i),s2);
sH=zeros(1,r);
sH(sHindx)=SimplexProj(sHitemp');
H(:,i)=sH';
end
 V=W*H;


%--------------------------------------------------
maxit=400;
 Ks1=s1;
 Ks2=s2;
 

 
 %---------------ARTMsparse----------------------------------------------
%  [AW,AH,Alls,Atime,Aerror,Aiter] = ARTMsparse(V,r,maxit);
 

%If the ARTMsparse algorithm is used, the initial point has been generated and only needs to be loaded.
%2.If the ARTMsparse algorithm is not used, the initial point needs to be generated randomly.

%-----1-----
% load('W0.mat');
% load('H0.mat');


%-----2-----
 W0=rand(m,r);
 for j=1:r
     W0(:,j)=W0(:,j)/sum(W0(:,j));
 end
 H0=rand(r,n);
 for i=1:n
 H0(:,i)=H0(:,i)/sum(H0(:,i));
 end

 
%-----------------------Alg. 2.2 for SSMF-------------------------------
[TW,TH,Ttime,Terror,Titera]=PALMSSMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2);
 sparw1=sum(sum(TW~=0));
 sparh1=sum(sum(TH~=0));



%--------------------------Alg. 2.4 for SSMF----------------------------
[PW,PH,Ptime,Perror,Pitera]=ALG24SSMF(W0,H0,V,maxit,m,n,r,Ks1,Ks2);

sparw2=sum(sum(PW~=0));
 sparh2=sum(sum(PH~=0));


%% Comparision of Alg. 2.2 and Alg. 2.4 
figure(1) 
semilogy((1:Titera),Terror(1:Titera),'r-',(1:Pitera),Perror(1:Pitera),'b-','MarkerSize',6,'LineWidth',1.5)
legend('Alg. 2.2','Alg. 2.4')
xlabel('Total number of iterations')
ylabel('$${H(V,W^kH^k)}$$','Fontsize',14,'FontWeight','bold','interpreter','latex') 
 
%%Comparision of Alg. 2.2, Alg. 2.4, and ARTMsparse
 
% figure (2)
% semilogy((1:Tiera),Terror(1:Titera),'r-',(1:Pitera),Perror(1:Pitera),'b-', (1:Aiter), Aerror(1:Aiter)'m-', 'MarkerSize',6,'LineWidth',1.5)
% legend('Alg. 2.2','Alg. 2.4','ARTMsparse')
% xlabel('Total number of iterations')
% ylabel('$${H(V,W^kH^k)}$$','Fontsize',14,'FontWeight','bold','interpreter','latex') 