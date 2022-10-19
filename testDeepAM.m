clc
clear all
close all

warning off
addpath('./utilis/','./goalplus/');

Path            = './models/sigma_0_IPADAtomNum1_35/';% path of the learned model

if_show         = 1;% if show dictionaries
if_SR           = 5;% 5 or 14 for Set5 and Set14
if_rescale      = 0;% if rescale the thresholds (0 or 1), default 0
sigma           = 0;% testing noise level

%% Dictionary Setting
d_sz            = [36,256,256,256];% dictionary size
sf              = 2;% up-sampling factor
sd              = 1;% stride

psz             = sqrt(d_sz(1));
%% Load Dictionary
load_D = 1;
if load_D == 1% Load Dictionaries
    OmegaNm   = sprintf([Path,'Omega_sf_%d_ps_%d_%d_%d_%d_%d.mat'],sf,d_sz);
    TNm       = sprintf([Path,'T_sf_%d_ps_%d_%d_%d_%d_%d.mat'],sf,d_sz);
    InfNm       = sprintf([Path,'Inf_sf_%d_ps_%d_%d_%d_%d_%d.mat'],sf,d_sz);

    load(OmegaNm,'Omega');    
    load(TNm,'T'); 
    load(InfNm,'Inform'); 
end
%% Rescale the thresholds
if if_rescale
    p = sigma/Inform{2};
    T{1}(1:Inform{1})       = T{1}(1:Inform{1})*p^2;
    T{1}(Inform{1}+1:end)   = T{1}(Inform{1}+1:end)*p;
end
%% Show Dictionaries
if if_show
    showDict(Omega)
end

%% Perform Super-Resolution
if if_SR==5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Super-Resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [imS1,imB1,psnr1]=DeepAM_SR(psz,sf,'./data/testing/Set5/baby.bmp',Omega,T,sd,sigma);
    [~,~,psnr2]=DeepAM_SR(psz,sf,'./data/testing/Set5/bird.bmp',Omega,T,sd,sigma);
    [~,~,psnr3]=DeepAM_SR(psz,sf,'./data/testing/Set5/butterfly.bmp',Omega,T,sd,sigma);
    [~,~,psnr4]=DeepAM_SR(psz,sf,'./data/testing/Set5/head.bmp',Omega,T,sd,sigma);
    [~,~,psnr5]=DeepAM_SR(psz,sf,'./data/testing/Set5/woman.bmp',Omega,T,sd,sigma);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PSNR = [psnr1,psnr2,psnr3,psnr4,psnr5]'
    mPSNR = mean(PSNR)
elseif if_SR == 14
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Super-Resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,~,psnr1]=DeepAM_SR(psz,sf,'./data/testing/Set14/baboon[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr2]=DeepAM_SR(psz,sf,'./data/testing/Set14/barbara[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr3]=DeepAM_SR(psz,sf,'./data/testing/Set14/bridge[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr4]=DeepAM_SR(psz,sf,'./data/testing/Set14/coastguard[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr5]=DeepAM_SR(psz,sf,'./data/testing/Set14/comic[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr6]=DeepAM_SR(psz,sf,'./data/testing/Set14/face[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr7]=DeepAM_SR(psz,sf,'./data/testing/Set14/flowers[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr8]=DeepAM_SR(psz,sf,'./data/testing/Set14/foreman[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr9]=DeepAM_SR(psz,sf,'./data/testing/Set14/lenna[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr10]=DeepAM_SR(psz,sf,'./data/testing/Set14/man[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr11]=DeepAM_SR(psz,sf,'./data/testing/Set14/monarch[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr12]=DeepAM_SR(psz,sf,'./data/testing/Set14/pepper[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr13]=DeepAM_SR(psz,sf,'./data/testing/Set14/ppt3[1-Original].bmp',Omega,T,sd,sigma);
    [~,~,psnr14]=DeepAM_SR(psz,sf,'./data/testing/Set14/zebra[1-Original].bmp',Omega,T,sd,sigma);
    toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PSNR = [psnr1,psnr2,psnr3,psnr4,psnr5,psnr6,psnr7,psnr8,psnr9,psnr10,...
        psnr11,psnr12,psnr13,psnr14]'
    mPSNR = mean(PSNR)
end