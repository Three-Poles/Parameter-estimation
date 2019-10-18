function [ET,Es,T,Ebs]=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss)
% Ta   air temperature in oC
global LAIe
global fr
global rabs
%% calculated rsc
thetamax=.49;
% thetamin=.11;
k3=1;
k4=0;
% KA=0.45;


F1=Rn./(Rn+k1);
a=(Tmax-k2)/(k2-Tmin);
F2=((Ta-Tmin).*(Tmax-Ta).^a)./((k2-Tmin).*(Tmax-k2).^a);
F3=1./(1+(D/D50).^k3)*(1-k4)+k4;
F4=((theta-thetamin).*(thetamax-thetamin+k5))./((thetamax-thetamin).*(theta-thetamin+k5));
dmon=F1.*F2.*F3.*F4.*LAIe;
rsc=rSTmin./(dmon+eps);

%% Clumping model
Cp=1.013;  
Ra=(delta+gamma).*raa;
Rc=(delta+gamma).*rac+gamma.*rsc;
Rs=(delta+gamma).*ras+gamma.*rss;
Rbs=(delta+gamma).*rabs+gamma.*rss;

fenmu=(Rs.*Rc.*Rbs+(1-fr).*Rs.*Rc.*Ra+fr.*Rbs.*Rs.*Ra+fr.*Rbs.*Rc.*Ra);
Cc=Rbs.*Rs.*(Rc+Ra)./fenmu;
Cs=Rbs.*Rc.*(Rs+Ra)./fenmu;
Cbs=Rs.*Rc.*(Rbs+Ra)./fenmu;

Rns=Rn.*exp(-KA.*LAI);

A=Rn-G;
As=Rns-G;
Ac=Rn-Rns;
Abs=Rn-G;

ETc=(delta.*A+(rho.*Cp.*D-delta.*rac.*As)./(raa+rac))./(delta+gamma.*(1+rsc./(raa+rac)));

ETs=(delta.*A+(rho.*Cp.*D-delta.*ras.*Ac)./(raa+ras))./(delta+gamma.*(1+rss./(raa+ras)));

ETbs=(delta.*Abs+rho.*Cp.*D./(raa+rabs))./(delta+gamma.*(1+rss./(raa+rabs)));

Eps=fr.*Cs.*ETs;
Ebs=(1-fr).*Cbs.*ETbs;
Es=Eps+Ebs;
T=fr.*Cc.*ETc;
ET=Eps+Ebs+T;

