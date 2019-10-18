利用Metropolis-Hasting (M-H)算法估计模型参数后验分布，这里以双源蒸散发发模型(S-W)为例进行算法应用演示。该程序包含2个主要的模块：

\1.    sw.m：双源蒸散发模型

[ET,Es,T,Ebs]=SW(KA,rSTmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,flag,Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss)

其中：KA,rstmin,D50,Tmin,Tmax,thetamin,k1,k2,k3,k4,k5,为模型参数；

KA:消光系数

rsTmin：最小气孔阻力

D50:水汽压有关参数

thetamin：土壤水分参数

Tmin,Tmax：温度控制参数

k1,k2,k3,k4,k5：为阻力参数

\-------------------------------------------------------------------

Ta,G,Rn,LAI,rho,D,theta,delta,gamma,raa,ras,rac,rss为输入的驱动数据。

Ta：空气温度

G:土壤热通量 [W m-2]

Rn：净辐射[W m-2]

LAI：叶面积指数

rho：空气密度

D：空气水汽压

theta：土壤水分

delta：饱和水汽压对温度的斜率

gamma：湿度计常数

raa,ras,rac,rss：空气动力阻力

\2.    main_code.m：主程序

在主程序中，输入驱动数据data，其行数为观测数据的个数，列为不同的观测项。从左到右分别为：风速、降雨量、空气温度、相对湿度、水汽气压、土塘热通量、向下短波辐射、向上短波辐射、向下长波辐射、向上长波辐射、大气压、蒸发、感热、摩擦风速、LAI、植被盖度。

计算相应的中间数据；最后进入M-H算法迭代过程，迭代次数为nsim=30000;

迭代完成后，统计参数的后验分布，完成模型参数估计

 

 