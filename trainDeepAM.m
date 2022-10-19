clc
clear all
close all

addpath('./utilis/','./goalplus/');

rng(12345,'twister');
%% Parameters for Saving Model
Path        = './models/sigma_0_IPADAtomsNum1_35/';% path of the learned model         
mkdir(Path);
%% Extract Patches for Training
d_sz        = [36,256,256,256];% defines the number of atoms in each layer
layers      = size(d_sz,2) - 1;% number of layers

sf          = 2;% super-resolution factor
sigma       = 0;% noise level
stride      = 2;% stride for patch extraction
Images      = './data/training/ScsR/';% training data path

[XH, XL]    = trainData(Images, sqrt(d_sz(1)), sf, stride, sigma);
XL          = dataDeMean(XL);
XH          = dataDeMean(XH);

IPADAtomNum1= 0;
% set the number of IPAD atoms in 1st layer, if set to 0, IPADAtomNum1 =
% ranhk of the input data

%% Testing and Validation Data
num = 1e5;
Index1 = randi([1 size(XL,2)],num,1);
Index2 = randi([1 size(XL,2)],num,1);

SlT = XL(:,Index1);      ShT = XH(:,Index1);% Testing Set
SlV = XL(:,Index2);      ShV = XH(:,Index2);% Validation Set

%% Parameters for Save Batch Data
PATH         = './batchdata/';          
mkdir(PATH);
param.PATH   = PATH;                    
param.bsize  = 4e4;% batch size
param.bnum   = 15;% batch number                       
Operator_type   = {'RAND'};

%% Start Learning
Xi              = XL;
OmegaEff        = 1;
for l  = 1:layers
    fprintf('\nLearning Layer-%d Analysis Dictionary\n', l);
    saveBatchData(Xi, XH, param);
    
    Pos(l)          = rank(single(XL*XL'));
    K(l)            = d_sz(l+1) - Pos(l);
    
    if l==1 && IPADAtomNum1>Pos(l)
        K(l) = K(l) - (IPADAtomNum1 - Pos(l));
    end
    
    %% Init IPAD
    OmegaIPAD{l}    = initDict(Operator_type, Xi, sqrt(d_sz(l)),(d_sz(l+1)- K(l))/d_sz(l));
    paramIPAD       = setParamIPAD(OmegaIPAD{l},sqrt(d_sz(l)),param);
    
    %% Init CAD
    OmegaCAD{l}     = initDict(Operator_type, XH, sqrt(size(XH,1)),K(l)/(size(XH,1)));
    paramCAD        = setParamCAD(OmegaCAD{l},sqrt(d_sz(l)),param);
    
    %% Learning IPAD
    fprintf('\nLearning IPAD\n');
    tmp = Xi;
    if l > 1
        tmp(size(OmegaIPAD{l-1},1)+1:end,:) = 0;
    end
    [Ui,~,~]        = svd(tmp*tmp');
    
    paramIPAD       = goalPlusIPAD(Pos(l), Ui, paramIPAD);
    OmegaIPAD{l}    = paramIPAD.Omega;
    ThresIPAD{l}    = zeros(size(OmegaIPAD{l},1),1);
    
    if d_sz(l+1) - K(l) > Pos(l)
        ThresIPAD{l}    = findThreRatio(OmegaIPAD{l}, Xi);
        [ThresIPAD{l},~]= threshScaleIPAD(OmegaIPAD{l}, SlT, ShT, SlV, ShV, ThresIPAD{l});
    end
    
    %% Learning CAD
    fprintf('\nLearning CAD\n');
    if K(l) > 0
        Dtmp        = XH*Xi'*pinv(Xi*Xi' + 1e0*eye(size(Xi*Xi')));
        [Uh,~,~]    = svd(XH*XH');
        PosCAD(l)   = rank(single(XH*XH'));
        paramCAD    = goalPlusCAD(PosCAD(l), Uh, Dtmp, paramCAD);
        
        OmegaCAD{l} = paramCAD.Omega*Dtmp;
        OmegaCAD{l} = bsxfun(@times,OmegaCAD{l},1./sqrt(sum(OmegaCAD{l}.^2,2)));
        ThresCAD{l} = findThreRatio(OmegaCAD{l}, Xi);
    else
        OmegaCAD{l} = [];
        ThresCAD{l} = [];
    end
    
    Omega{l} = [OmegaIPAD{l}; OmegaCAD{l}];
    Thres{l} = [ThresIPAD{l}; ThresCAD{l}];
    
    if K(l) > 0
        [Thres{l},S] = threshScaleCAD(Omega{l}, SlT, ShT, SlV, ShV, Thres{l}, 1./ThresCAD{l});
    end
    
    %% Prepare Data for the next layer learning
    Xi      = wthreshVec(Omega{l}*Xi,'s',Thres{l});
    SlV     = wthreshVec(Omega{l}*SlV,'s',Thres{l});
    SlT     = wthreshVec(Omega{l}*SlT,'s',Thres{l});
end

if layers == 0
    l=0;
    Thres = {};
end
%%
fprintf('\nLearning Synthesis Dictionary\n', l)
D           = XH*Xi'*pinv(Xi*Xi' + 0*eye(size(Xi*Xi')));
Omega{l+1}  = D;
T           = Thres;

Inform{1}   = size(OmegaIPAD{1},1);
Inform{2}   = sigma;

OmegaNm     = sprintf([Path,'Omega_sf_%d_ps_%d_%d_%d_%d.mat'],sf,d_sz);
TNm         = sprintf([Path,'T_sf_%d_ps_%d_%d_%d_%d.mat'],sf,d_sz);
InfNm       = sprintf([Path,'Inf_sf_%d_ps_%d_%d_%d_%d.mat'],sf,d_sz);

save(OmegaNm,'Omega');
save(TNm,'T');
save(InfNm,'Inform')

showDict(Omega)

