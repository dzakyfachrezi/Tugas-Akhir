function [mECG_ab,best_MSE,opt_mu,opt_nC] = LMS_extraction(mECG,aECG,ext_row)

% Parameter adaptive filter
grid_mu = [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]; % range of step size
grid_nC = [1,5,10,15,20,25,30,40,50,60]; % range of filter length
opt_mu = 0; % the optimum step size after grid search
opt_nC = 0; % the optimum filter length after grid search
best_MSE = Inf; % the highest SNR

% Normalize signal (scaling berdasarkan segmen 10 detik pertama)
mECG_norm = normecg(mECG,1,(ext_row*2));
aECG_norm = normecg(aECG,1,(ext_row*2));

% Find the optimum mu and nC using grid search
for mu = grid_mu
    for nC = grid_nC
        [~,e] = LMS_adaptivefiltering(mECG_norm,aECG_norm,mu,nC);
        e = e(ext_row+1:end);
        MSE = mean(abs(e).^2);
        if MSE < best_MSE
            best_MSE = MSE;
            opt_mu = mu;
            opt_nC = nC;
        end
    end
end

% Compute residual using the optimum mu and nC
[y_opt,~] = LMS_adaptivefiltering(mECG_norm,aECG_norm,opt_mu,opt_nC);
mECG_ab = y_opt(ext_row+1:end);

end