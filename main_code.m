clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             Bayesian estiamte S-W parameter for Maize in Linze Station
%
%                                 ^    ^
%                                   --
%                                 Take easy Baby
%
%              2012-10-16 at Lanzhou University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% WS	Rain	Ta	Rh	Vapor	SoilFLUX  D_SR	U_SR	D_LR	U_LR  Vwc  P_kPa	 LE	 H	  Fc	   LAI  f
% 1      2      3   4    5        6        7     8       9       10   11    12       13  14   15       16  17


load dataAug
dataApril=dataAug;

index=find(dataApril(:,7)~=-999 & dataApril(:,11)~=-999 & dataApril(:,12)~=-999 & dataApril(:,13)~=-999 & dataApril(:,1)>1);


data=dataApril(index,:);

global fr
%% Observed data
P=data(:,12);                                     % atmospheric pressure in kPa
u=data(:,1);                                      % windspeedat10m
Ta=data(:,3);                                     % air temperature in oC
theta=data(:,11);                                  % soil watetr content
RH=data(:,4);                                     %relative humidity
G=data(:,6);                                      %soil heat flux
Rn=data(:,7)+data(:,9)-data(:,8)-data(:,10);     %in w/m2
LEEC=data(:,13);
H=data(:,14);
LAI=data(:,16);
fr=data(:,17);
global LAIe

LAIe=zeros(size(data,1),1);
for i=1:size(data,1)
    if LAI(i)<=2
        LAIe(i,1)=LAI(i);
    elseif LAI(i)<=4
        LAIe(i,1)=2;
    else
        LAIe(i,1)=LAI(i)/2;
    end
end


%% Climatic variables
thetas=.49;
lambda=2500.78-2.3601*Ta;     

% saturaed vapour pressure in kPa
es=.6108*exp(17.27*Ta./(Ta+237.3));
% slope of pressure to temperature
delta=4098*es./(Ta+237.3).^2;

%air density
Ta=Ta+273.14;                   % convert oC to K
Rd=287/1000;                    % the gas constant in kJ/kg/K
ea=es.*RH/100;
D=es-ea;
rho=P./(Rd*Ta.*(1+.378*ea./P));

%
Cp=1.013;                %specific heat capacity of the dry air in kJ/kg/K;
epsilong=.622;           %the ratio between the mplecular weight of water vapor and air

gamma=Cp*P./(lambda*epsilong);


%% calculate raa
hc=.25;                                         % crop height in [m]
k=0.41;                                         % von K¨¢rm¨¢n constant
z0=0.13*hc;                                     % roughness length in [m];
d=0.63*hc;                                      % zero plane displacement in [m]
z0h=0.1*z0;                                     % roughness length to the heat flux in [m];
n=2.5;                                          % parameter in SW model
z=4;                                            % reference height

% LAI >4 
raa_inf=log((z-d)/z0)./(k*k*(u+eps))*(log((z-d)/(hc-d))+hc/(n*(hc-d))*(exp(n*(1-(d+z0)/hc))-1)); 

% for bare surface
z0s=0.01;
ras0=log(z/z0s)*log((d+z0)/z0s)./(k*k*(u+eps));
raa_bare=(log(z/z0s)*log(z/z0s))./(k*k*(u+eps))-ras0;     % bare surface

% 
raa=.25*LAI.*raa_inf+.25*(4-LAI).*raa_bare;

%% calculate ras

% LAI >4 
ras_inf=log((z-d)/z0)./(k*k*(u+eps))*hc/(n*(hc-d))*(exp(n)-exp(n*(1-(d+z0)/hc)));

% for bare surface
z0s=0.01;
ras_bare=log(z/z0s)*log((d+z0)/z0s)./(k*k*(u+eps));

ras=.25*LAI.*ras_inf+.25*(4-LAI).*ras_bare;

%% calculated rac
zw=10;                                                    % observation height 2m
z0w=.05;                                                 % roughness length
Fw=5000;                                                 % the fetch at weather station
zb=.33*Fw^.875*z0w.^.125;                                % the height of internal boundary layer
z0g=.01;                                                 % roughness length of ground
if hc<=1
    z0c=.13*hc;
elseif hc<10
    z0c=.139*hc-.009*hc*hc;
else
    z0c=.05*hc;
end

if hc==0
    Cd=1.4/1000;
else
    Cd=(-1+exp(.909-3.03*z0c/hc)).^4/4;
end

if LAI>=4
    d0=hc-z0c/.3;
else
    d0=1.1*hc*log(1+(Cd.*LAI).^(.25));
end

z0=min(.3*(hc-d0),z0g+.3*hc*sqrt(Cd.*LAI));

za=hc;
uh=u.*(log(zb./z0w)./log(zb./z0)).*(log((za-d0)./z0)./log(zw./z0w));

n=2.5;
w=.02;                                                          % leaf width
rb=100/n*sqrt(w./(uh+eps))./(1-exp(-n/2));
rac=rb./LAI/2;


%% rss
thetas=.49;
a=7.02;
b=2.27;
rss=exp(a-b*theta./thetas);

%% rabs
global rabs
zm=0.75*hc;
au=3;
um=uh*exp(au*(zm/hc-1));
rabs=log(zm/z0g)*log(zm/z0g)./(k*k*um);



%% 
% parameter KA,rSTmin,D50,Tmin,Tmax,thetamin,k1 k2  k3  k4 k5,sigma
%           1    2     3    4     5     6     7  8   9  10 11 12
cmin(1)=.2;cmin(2)=10;cmin(3)=0.5;cmin(4)=0;cmin(5)=35;
cmin(6)=0.05;cmin(7)=0;cmin(8)=20;cmin(9)=1;cmin(10)=0;cmin(11)=.02;cmin(12)=0;

cmax(1)=.9;cmax(2)=1000;cmax(3)=2.5;cmax(4)=10;cmax(5)=40;
cmax(6)=.25;cmax(7)=900;cmax(8)=30;cmax(9)=10;cmax(10)=1;cmax(11)=.20;cmax(12)=100;

KA=cmin(1)+(cmax(1)-cmin(1))*rand;
rSTmin=cmin(2)+(cmax(2)-cmin(2))*rand;
D50=cmin(3)+(cmax(3)-cmin(3))*rand;
Tmin=cmin(4)+(cmax(4)-cmin(4))*rand;
Tmax=cmin(5)+(cmax(5)-cmin(5))*rand;
thetamin=cmin(6)+(cmax(6)-cmin(6))*rand;
k1=cmin(7)+(cmax(7)-cmin(7))*rand;
k2=cmin(8)+(cmax(8)-cmin(8))*rand;
k3=cmin(9)+(cmax(9)-cmin(9))*rand;
k4=cmin(10)+(cmax(10)-cmin(10))*rand;
k5=cmin(11)+(cmax(11)-cmin(11))*rand;
sigma=cmin(12)+(cmax(12)-cmin(12))*rand;



%% Initialize the Metropolis-Hastings sampler
T      = 10000;                                  % Maximum number of iterations
burnin = 50;                                     % burnin period
thin   = 5;                                      % thinning parameter

%% Storage space for our samples

Ta=data(:,3);                                     % air temperature in oC
state=zeros(12,T);
accept=zeros(12,T);

s=1;

for t=1:T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  KA 1
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_KA=KA;
    while true
        new_KA=old_KA+(-.5+rand)*(cmax(1)-cmin(1))/5;
        if new_KA<cmax(1) && new_KA>cmin(1)
            break
        end
    end    
    Model_old=SW(old_KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);    
    Model_new=SW(new_KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
   
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        KA=new_KA;
        accept(1,t)=1;
    else
        KA=old_KA;
    end
    state(1,t)=KA;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  rSTmin 2
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_rSTmin=rSTmin;
    while true
         new_rSTmin=old_rSTmin+(-.5+rand)*(cmax(2)-cmin(2))/10;
        if new_rSTmin<cmax(2) & new_rSTmin>cmin(2)
            break
        end
    end
        
    Model_old=SW(KA,old_rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,new_rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));
    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        rSTmin=new_rSTmin;
        accept(2,t)=1;
    else
        rSTmin=old_rSTmin;
    end
    state(2,t)=rSTmin;
 
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  D50 3
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_D50=D50;
    while true
         new_D50=old_D50+(-.5+rand)*(cmax(3)-cmin(3))/10;
        if new_D50<cmax(3) & new_D50>cmin(3)
            break
        end
    end     
   
    Model_old=SW(KA,rSTmin,old_D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,new_D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        D50=new_D50;
        accept(3,t)=1;
    else
        D50=old_D50;
    end
    state(3,t)=D50;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  Tmin 4
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_Tmin=Tmin;
    while true
        new_Tmin=old_Tmin+(-.5+rand)*(cmax(4)-cmin(4))/5;
        if new_Tmin<cmax(4) & new_Tmin>cmin(4)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,old_Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,new_Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        Tmin=new_Tmin;
        accept(4,t)=1;
    else
        Tmin=old_Tmin;
    end
    state(4,t)=Tmin;
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  Tmax 5
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_Tmax=Tmax;
    while true
        new_Tmax=old_Tmax+(-.5+rand)*(cmax(5)-cmin(5))/5;
        if new_Tmax<cmax(5) & new_Tmax>cmin(5)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,old_Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,new_Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        Tmax=new_Tmax;
        accept(5,t)=1;
    else
        Tmax=old_Tmax;
    end
    state(5,t)=Tmax;
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  thetamin 6
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_thetamin=thetamin;
    while true
        new_thetamin=old_thetamin+(-.5+rand)*(cmax(6)-cmin(6))/5;
        if new_thetamin<cmax(6) & new_thetamin>cmin(6)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,old_thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,new_thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        thetamin=new_thetamin;
        accept(6,t)=1;
    else
        thetamin=old_thetamin;
    end
    state(6,t)=thetamin;
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  k1 7
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_k1=k1;
    while true
        new_k1=old_k1+(-.5+rand)*(cmax(7)-cmin(7))/5;
        if new_k1<cmax(7) & new_k1>cmin(7)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,old_k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,new_k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        k1=new_k1;
        accept(7,t)=1;
    else
        k1=old_k1;
    end
    state(7,t)=k1;
    
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  k2     8
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_k2=k2;
    while true
        new_k2=old_k2+(-.5+rand)*(cmax(8)-cmin(8))/5;
        if new_k2<cmax(8) & new_k2>cmin(8)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,old_k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,new_k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        k2=new_k2;
        accept(8,t)=1;
    else
        k2=old_k2;
    end
    state(8,t)=k2;
    
    
    
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  k3     9
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_k3=k3;
    while true
        new_k3=old_k3+(-.5+rand)*(cmax(9)-cmin(9))/5;
        if new_k3<cmax(9) & new_k3>cmin(9)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,old_k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,new_k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        k3=new_k3;
        accept(9,t)=1;
    else
        k3=old_k3;
    end
    state(9,t)=k3;
    
    
    
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  k4    10
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_k4=k4;
    while true
        new_k4=old_k4+(-.5+rand)*(cmax(10)-cmin(10))/5;
        if new_k4<cmax(10) & new_k4>cmin(10)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,old_k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,new_k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        k4=new_k4;
        accept(10,t)=1;
    else
        k4=old_k4;
    end
    state(10,t)=k4;
    
    
    
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  k5    11
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_k5=k5;
    while true
        new_k5=old_k5+(-.5+rand)*(cmax(11)-cmin(11))/5;
        if new_k5<cmax(11) & new_k5>cmin(11)
            break
        end
    end
    
    Model_old=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,old_k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    Model_new=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,new_k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    L_old_gmax=sum(log(normpdf(Model_old,LEEC,sigma)));    
    L_new_gmax=sum(log(normpdf(Model_new,LEEC,sigma)));    
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        k5=new_k5;
        accept(11,t)=1;
    else
        k5=old_k5;
    end
    state(11,t)=k5;
       
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Sampling  sigma 12
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    old_sigma=sigma;
    while true
        new_sigma=old_sigma+(-.5+rand)*100/5;
        if new_sigma<100 & new_sigma>0
            break
        end
    end
    
    Model=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,[],Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
    
    
    L_old_gmax=sum(log(normpdf(Model,LEEC,old_sigma)))-2*log(old_sigma);
    
    L_new_gmax=sum(log(normpdf(Model,LEEC,new_sigma)))-2*log(new_sigma); 
    if mod(t,100)==0
        [ET1(:,s),Es1(:,s),T1(:,s),Ebs1(:,s)]=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,[],Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);
        s=s+1;
    end
    
    ratio=(L_new_gmax)-L_old_gmax;
    % Accept or reject?
    r=log(unifrnd(0,1));
    if r<ratio
        sigma=new_sigma;
        accept(12,t)=1;
    else
        sigma=old_sigma;
    end
    state(12,t)=sigma;
    
end



%% estiamte parameter
parameter=[];
for i=1:12
    parameter(i,:)=prctile(state(i,5000:end),[2.5 50 97.5]);
end
save parameter

[ET,Es,T,Ebs]=SW(parameter(1,2),parameter(2,2),parameter(3,2),parameter(4,2),parameter(5,2),parameter(6,2),parameter(7,2),parameter(8,2),parameter(9,2),parameter(10,2),parameter(11,2),[],Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss);

plot(ET,'r')
hold on
plot(LEEC,'o')
plot(Es)
plot(T,'g')
plot(Ebs,'k')
figure
plot(LEEC,ET,'.')
hold on
plot(-100:800,-100:800)
xlim([-100 800])
ylim([-100 800])