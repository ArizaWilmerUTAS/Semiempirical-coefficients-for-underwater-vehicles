%this script shows how everything works together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wilmer Ariza Ramirez
%Implementarion of calculation fo Preterus Thesis for Remus 100
%7/11/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read Stl file containt file to generate R(x) for integration center of figure should be CB
% we export a 3D figure and slide a layer to be used as 2d Template
% file has to be configure as ascii STL
clc; clear vars; close all;

triangles = read_ascii_stl('Remus2.stl',1);
original = triangles;
%triangles = orient_stl(triangles,'z');
%triangles = rotate_stl(triangles,'y',90);
slice_height = 4;
tic;[movelist, z_slices] = slice_stl_create_path(triangles, slice_height);toc;
Path=movelist{2};

plot_slices(movelist,z_slices, 0.0001)

daspect([1 1 1])
%plot capture data
figure(2)
Path(:,1)=-(Path(:,1)-611)/1000;
Path(:,2)=Path(:,2)/1000;
plot(Path(:,1),Path(:,2));
daspect([1 1 1])

%% Interpolation
Xinter = linspace(min(Path(:,1)),max(Path(:,1)),5000);

Path(isnan(Path(:,1)),:) = [];
index=ind2sub([233,1],find(Path(:,2)<0));
Path(index,:) = [];
while i<length(Path)
    if round(Path(i-1,1),2)==round(Path(i,1),2)
        Path(i,:)=[];
    else
        i=i+1;
    end
    
end

% indices to unique values in column 3
[~, ind] = unique(Path(:,1), 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(Path(:,1), 1), ind);

[x, index] = unique(Path(:,1)); 
y=Path(index,2);
I=find(y<0.0005);
x(I)=[];
y(I)=[];
plot(x,y,'*');
Yinter = interp1(x, y, Xinter); 

 
figure(2);
plot(Xinter,Yinter,'*');
daspect([1 1 1])
[m,n]=size(Path);
i=2;
%% REMUS parameters
d=1.91e-1;%vehicle diameter
l=1.33; %vehicle lenght
ro=1000; %Water density
afin=0.131;%ccorected from preterus
xt=-0.721;%aft end of tail section
xf=-0.685;%aft end of fin section
xf2=-0.611;%%fordward end of fin section
xb2=0.610;%forward end of bow section
Sfin=6.65e-3;%Fin Platform area
Xfin=6.38e-1;
%% Xudot
betha= 0.251892;%Blevins aspect ratio
Xudot=-4*betha*ro*pi/3*(d/2)^3;
%% Yvdot
I=find(Xinter>xt & Xinter<xf);

X=Xinter(I);
Y=Yinter(I);
ma_x=pi*ro*Y.^2;
Yvdot=-trapz(X,ma_x);

I=find(Xinter>xf & Xinter<xf2);
X=Xinter(I);
Y=Yinter(I);
maf_x=pi*ro*(afin^2-Y.^2+Y.^4./afin^2);
Yvdot=Yvdot-trapz(X,maf_x);

I=find(Xinter>xf2 & Xinter<xb2);
X=Xinter(I);
Y=Yinter(I);
I=find(isnan(Y) == 1);
X(I)=[];
Y(I)=[];
ma_x=pi*ro*Y.^2;
Yvdot=Yvdot-trapz(X,ma_x);

%% Zwdot
Zwdot=Yvdot;
%% Mwdot
I=find(Xinter>xt & Xinter<xf);

X=Xinter(I);
Y=Yinter(I);
ma_x=X.*pi*ro.*Y.^2;
Mwdot1=trapz(X,ma_x);

I=find(Xinter>xf & Xinter<xf2);
X=Xinter(I);
Y=Yinter(I);
maf_x=X.*pi*ro.*(afin^2-Y.^2+Y.^4./afin^2);
Mwdot2=trapz(X,maf_x);

I=find(Xinter>xf2 & Xinter<xb2);
X=Xinter(I);
Y=Yinter(I);
I=find(isnan(Y) == 1);
X(I)=[];
Y(I)=[];
ma_x=X.*pi*ro.*Y.^2;
Mwdot3=trapz(X,ma_x);
Mwdot=-Mwdot1+Mwdot2-Mwdot3;
%% Nvdot
Nvdot=-Mwdot;
%% Yrdot
Yrdot=Nvdot;
%% Zqdot
Zqdot=Mwdot;
%% Mqdot
I=find(Xinter>xt & Xinter<xf);

X=Xinter(I);
Y=Yinter(I);
ma_x=X.^2.*pi*ro.*Y.^2;
Mqdot1=trapz(X,ma_x);

I=find(Xinter>xf & Xinter<xf2);
X=Xinter(I);
Y=Yinter(I);
maf_x=X.^2.*pi*ro.*(afin^2-Y.^2+Y.^4./afin^2);
Mqdot2=trapz(X,maf_x);

I=find(Xinter>xf2 & Xinter<xb2);
X=Xinter(I);
Y=Yinter(I);
I=find(isnan(Y) == 1);
X(I)=[];
Y(I)=[];
ma_x=X.^2.*pi*ro.*Y.^2;
Mqdot3=trapz(X,ma_x);
Mqdot=-Mqdot1-Mqdot2-Mqdot3;

%% Nrdot
Nrdot=Mqdot;
%% Kpdot

Kpdot=2/pi*ro*0.1172^4*(xf2-xf);
%% added mass Cross terms

Xwq=Zwdot;
Yura=Xudot;
Zuqa=-Xudot;
Muwa=-(Zwdot-Xudot);
Nuva=-(Xudot-Yvdot);
Xqq=Zqdot;
Ywp=-Zwdot;
Zvp=Yvdot;
Mvp=-Yrdot;
Nwp=Zqdot;
Xvr=-Yvdot;
Ypq=-Zqdot;
Zrp=Yrdot;
Mrp=(Kpdot-Nrdot);
Npq=-(Kpdot-Mqdot);
Xrr=-Yrdot;
Muqa=-Zqdot;
Nura=Yrdot;

%% Yuvl
Cydbetha=0.003*l/d*180/pi;


Yuvl=-1/2*ro*d^2*Cydbetha;
%% Zuwl
Zuwl=Yuvl;
%% Muwl
Xcp=-0.65*l+0.61;
Muwl=-0.5*ro*d^2*Cydbetha*Xcp;
%% Nuv
Nuvl=Muwl;
%% Fin lift
Clalpha=(1/(2*0.9*pi)+1/(pi*2.21))^-1;
ARe=2.21;
Yuudr=ro*Clalpha*6.65e-3/2;
Yuvf=-Yuudr;
Zuuds=-ro*Clalpha*6.65e-3/2;
Zuwf=Zuuds;
Yurf=-ro*Clalpha*6.65e-3/2*-0.638;
Zuqf=-Yurf;
%% Moment lift
Muuds=ro*Clalpha*6.65e-3/2*-0.638;
Muwf=Muuds;
Nuudr=ro*Clalpha*6.65e-3/2*-0.638;
Nuvf=Nuudr;
Muqf=-ro*Clalpha*6.65e-3/2*-0.638^2;
Nurf=Muqf;
%% combined terms
Yuv=Yuvl+Yuvf;
Yur=Yura+Yurf;
Zuw=Zuwl+Zuwf;
Zuq=Zuqa+Zuqf;
Muw=Muwa+Muwl+Muwf;
Muq=Muqa+Muqf;
Nuv=Nuva+Nuvl+Nuvf;
Nur=Nura+Nurf;
%% Axial Drag
Af=2.85e-2;%Vehicle fron area
Ap=l*d;
Css=3.397e-3;%% from principles of Naval architecture
Cd=Css*pi*Ap/Af*[1+60*(d/l)^3+0.0025*l/d];%aproximation of axial drag
Xuu=-0.5*ro*Cd*Af;

%% Crossflow drag
I=find(Xinter>xt & Xinter<xb2);

X=Xinter(I);
Y=Yinter(I);
I=find(isnan(Y) == 1);
X(I)=[];
Y(I)=[];
Cdc=1.1;%%Hoerner drag coefficient of a cylinder
Cdf=0.56;%Crossflow drag coefficient
Zww=-0.5*ro*Cdc*(trapz(X,2*Y))-2.*(0.5*ro*Sfin*Cdf);
Yvv=Zww;
Mww=0.5*ro*Cdc*(trapz(X,2*X.*Y))-2*Xfin*(0.5*ro*Sfin*Cdf);
Nvv=-Mww;
Yrr=-0.5*ro*Cdc*(trapz(X,2*X.*abs(X).*Y))-2*Xfin*abs(Xfin)*(0.5*ro*Sfin*Cdf);
Zqq=-Yrr;
Mqq=-0.5*ro*(trapz(X,2*Cdc*l*Y))-2*Xfin*(0.5*ro*Sfin*Cdf);
Nrr=Mqq;
%% Rolling Drag
yvvf=2.*(0.5*ro*Sfin*Cdf);
Kpp=yvvf*afin^3;
