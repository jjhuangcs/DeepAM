%% set the parameters for IPAD learning

function Learn_para = setParamIPAD(Omega,p_sz,param)

% LogAbs, LogSquare, PNormAbs, PNormSquare, AtanAbs, AtanSquare
Learn_para.Sp_type  = 'LogSquare';%'PNormAbs';%LogSquare LogAbs
Learn_para.verbose  = 1;
Learn_para.logger   = [];
Learn_para.Omega    = Omega;
Learn_para.p_sz     = p_sz;
Learn_para.max_iter = 100;

n = size(Omega,2);
k = size(Omega,1);

Learn_para.p        = 1;
Learn_para.q        = 2;%0 %% q = 2 to minimize the variance

Learn_para.kappa    = k*1*1e0;
Learn_para.mu       = 0;

Learn_para.nu       = 1e2*n;
Learn_para.trainSize= 1;

Learn_para.path  	= param.PATH;
Learn_para.tBatchNum= param.bnum;
end