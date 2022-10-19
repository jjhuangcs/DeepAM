% Function finds the best analysis operator Omega that lies on the oblique
% manifold via a riemannian conjugate gradient method. The cost function to
% be minimized is f(Omega) = ||Omega*Y||_p - kappa*log(det(Omega'*Omega))
% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de

%% Goal+ learning algorithm for CAD
function para = goalPlusCAD(Pos, Uy, D, para)

fprintf('Rank %d\n',Pos);
if Pos < size(Uy,1)
    E = Uy(:,Pos+1:end)';% E is the nullspace vectors of the input data
    Eall = Uy(:,Pos+1:end)';
else
    E = [];
    Eall = [];
end

U = Uy(:,1:Pos);
if isempty(E)
    U = eye(size(para.Omega,2));
end
%%
% Normalize the initial Operator
Omega       = bsxfun(@times,para.Omega,1./sqrt(sum(para.Omega.^2,2)));

if size(Omega,1) < size(Omega,2)
    U = Uy(:,1:size(Omega,1));
    Eall = Uy(:,size(Omega,1)+1:end)';
end

I           = eye(size(Omega,1)) + ones(size(Omega,1));
I_ext_grad  = I == 2;
%%
nu 			= para.nu;
q 			= [para.q,para.p];
Sp_type     = para.Sp_type;
kappa       = para.kappa*1/(log(size(para.Omega*U,2))*size(para.Omega*U,2));
mu          = para.mu;

% Linesearch parameters
alpha 		= 1e-1;
beta  		=   .9;
Qk          =    0;
Ck          =    0;
nuk         =   .0;

SUB_LOG_DET = size(Omega*U,2)*log(size(Omega,1));%klog(m)

numAtoms    = size(Omega,1);

numStop     = 0;
for iterBatch = 1:para.tBatchNum
    fprintf('Batch %d\n',iterBatch);
    %%
    miniY = [];
    miniX = [];
    for ld = 1:para.trainSize
        bNum = randi([1 para.tBatchNum],1,1);
        xh = load(sprintf([para.path,'XH_%d'],bNum));
        miniY = [miniY,xh.xh];
        
        xl = load(sprintf([para.path,'XL_%d'],bNum));
        miniX = [miniX,xl.xl];
        clear xh xl;
    end
    
    M           = numel(miniY);
    miniX       = D*miniX;
    miniR       = miniY-miniX;
    %%
    OY  = Omega*miniY;
    
    [f0_sp0,q_w0] = Sparsifying_functions_norm(Sp_type, 'Evaluate', OY, q, nu);
    f0             = 1/M*f0_sp0;
    
    [f0_sp2,q_w2] = Sparsifying_functions_norm_joint(Sp_type, 'Evaluate', Omega, miniR, miniX, q, nu);
    f0             = f0 + mu*1/M*f0_sp2;
        
    f0_sp = f0;
    if para.verbose
        Sp_input = f0;
    end
    
    if kappa ~= 0
        [~,Si,Vi]    = svd(Omega*U);
        LogDetTerm   = log(prod(diag(Si).^2))-SUB_LOG_DET;
        g_term       = -kappa*(LogDetTerm);
        f0           = f0 + 1*g_term;
    end
    
    for k = 1:para.max_iter
        
        if iterBatch==1 && k==2
            Sp_input = f0_sp;
        end
        
        %% Computing the gradient
        Omega_grad = zeros(size(Omega));
        
        Omega_grad0 = Sparsifying_functions_norm(Sp_type, 'Derivative', OY, q, nu, [], q_w0);
        Omega_grad  = Omega_grad + 1/M*Omega_grad0*miniY';
        
        Omega_grad2 = Sparsifying_functions_norm_joint(Sp_type, 'Derivative', Omega, miniR, miniX, q, nu, [], q_w2);
        Omega_grad2 = 2*(Omega_grad2.*(Omega*miniX)*miniX' - Omega_grad2.*(Omega*(miniR))*(miniR)');
        Omega_grad  = Omega_grad + mu*1/M*Omega_grad2;
        
        %         Omega_grad = Omega_grad - (2*knorm1*Omega*M1'*M1./(vecnorm(Omega*M1')) + 2*knorm2*Omega*M2'*M2./(vecnorm(Omega*M2')));
        
        if kappa ~= 0
            Omega_grad  = Omega_grad - 2*1*kappa*Omega*U*Vi*diag(1./diag(Si).^2)*Vi'*U';
        end
        
        %% Transpose for update step on the Oblique manifold
        Omega_grad = Omega_grad';
        Omega      = Omega';
        
        szO1 = size(Omega,1);
        for pp = 1:size(Omega,2)
            %% Old Equations before 2018-01-29
            A = [2*Omega(:,pp) Eall'];
            P = eye(szO1) - A/(A'*A)*A';
            % orthogonal projection onto A. This is much better! 2018-07-04
            Omega_grad(:,pp) = P*Omega_grad(:,pp);
        end
        
        if numel(Omega_grad(isnan(Omega_grad)))
            Omega_grad(isnan(Omega_grad)) = 1e10;
            return
        end
        
        %% Computation of Conjugate Direction
        if k == 1 || ~mod(k,numel(Omega))
            dx          = -Omega_grad;
            t           = 2*pi/max(sqrt(sum(dx.^2)));
            t_initial   = t;
        else
            tau_g       =  parallel_transport(g0, Omega_old, dx, t_prev, Norm_dx );
            tau_dx      =  parallel_transport([], Omega_old, dx, t_prev, Norm_dx);
            
            yk          =  Omega_grad - tau_g;
            denom       =  tau_dx(:)'*yk(:);
            
            cg_beta     = (Omega_grad(:)'*yk(:))/denom;
            cg_beta_dy  = norm(Omega_grad(:))^2/denom;
            cg_beta     = max(0,min(cg_beta_dy, cg_beta));
            
            dx          = -Omega_grad + cg_beta*tau_dx;
        end
        
        val      = Omega_grad(:)'*dx(:);
        Norm_dx  = sqrt(sum(dx.^2));
        sel      = Norm_dx > 0;
        ls_iter  = 0;
        
        Ck = (nuk*Qk*Ck+f0)/(nuk*Qk+1);
        Qk =  nuk*Qk+1;
        
        %% Backtracking Linesearch
        while (ls_iter == 0) || ... % Quasi Do while loop
                (f0 > Ck + alpha * t * val) && ... % Check Wolfe Condition
                (ls_iter < 100) % Check Maximum number of Iterations
            
            % Update the step size
            t  = t*beta;
            % Store the step length as it is required for the parallel
            % transport and the retraction
            t_prev    = t;
            
            % Exponential mapping and transpose back to get standing operator
            Omega_c = exp_mapping(Omega, dx, t, Norm_dx, sel)';
            % Evalute the costfunction
            OY = Omega_c*miniY;
            [f0_sp0,q_w0] = Sparsifying_functions_norm(Sp_type, 'Evaluate', OY, q, nu);
            f0            = 1/M*f0_sp0;
            
            [f0_sp2,q_w2] = Sparsifying_functions_norm_joint(Sp_type, 'Evaluate', Omega_c, miniR, miniX, q, nu);
            f0            = f0 + mu*1/M*f0_sp2;
            
            f0_sp         = f0;
            
            if kappa ~= 0
                [~,Si,Vi]  = svd(Omega_c*U);
                g_term     =  -kappa*(log(prod(diag(Si).^2))-SUB_LOG_DET);
                f0         = f0 + 1*g_term;
            end
            
            % Increase the number of linesearch iterates as we only
            % allow a certain amount of steps
            ls_iter = ls_iter + 1;
        end
        para.logger = [para.logger,f0];
        
        Omega = Omega';
        
        if  k==para.max_iter
            fprintf('Change of the operator: %f\n',norm(Omega-Omega_c,'fro'))
            fprintf('Current Sparsity %e ~ Input Sparsity %e\n', f0_sp, Sp_input)
            OO      = Omega_c*(eye(size(Omega_c,2))-1/size(Omega_c,2)*ones(size(Omega_c,2)));
            Sings   = svd(Omega_c);
            Cohe    = sort(abs(OO*OO'));
            fprintf('Condition: %f ~ Mutual Coherence: %f\n',Sings(1)/Sings(end-size(E,1)),max(Cohe(end-1,:)))
        end
        
        if  norm(Omega-Omega_c,'fro')/size(Omega,1)/size(Omega,2) < 1e-16
            fprintf('************ ITERATION %d **********\n',k);
            fprintf('************ NO PROGRESS **********\n');
            numStop     = numStop + 1;
            break;
        end
        
        para.Omega = Omega_c;
        
        if ls_iter == 100
            fprintf('************ ITERATION %d **********\n',k);
            fprintf('************ FINISHED **********\n');
            numStop     = numStop + 1;
            break;
        else
            Omega_old   = Omega';
            g0          = Omega_grad;
            Omega       = Omega_c;
            t           = t/beta^2;
        end
    end
    if numStop >= 3% default value 3
        break;
    end
end
end
