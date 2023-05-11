function [opt_kernel, opt_basis,min_score] = optimizer(x,y)

kernel_options = {'Matern32', 'Matern52','SquaredExponential','ardmatern32','ardmatern52','ardsquaredexponential'};
basis_options = {'None', 'Constant', 'Linear', 'PureQuadratic'};

% Define cross-validation parameters
K = 9; % number of folds
rng(1); % set random seed for reproducibility

% Initialize variables to store results
BIC_scores = zeros(length(kernel_options), length(basis_options));

% Loop over all combinations of kernel and basis functions
for i = 1:length(kernel_options)
    for j = 1:length(basis_options)
        % Define Gaussian Process model with current kernel and basis functions
        kernel = kernel_options{i};
        basis = basis_options{j};
        X = x';
        Y = y;
        gprMdl = fitrgp(X, Y, 'KernelFunction', kernel, 'BasisFunction', basis);
        
        
        % Perform K-fold cross-validation
        cv = cvpartition(size(X,1), 'KFold', K);
        cv_BIC_scores = zeros(K, 1);
        for k = 1:K
            trIdx = cv.training(k);
            teIdx = cv.test(k);
            gprMdl_k = fitrgp(X(trIdx,:), Y(:,trIdx), 'KernelFunction', kernel, 'BasisFunction', basis);
            yhat = predict(gprMdl_k,X(teIdx,:));
            r = yhat - Y;
            n = size(X,1);
            [~,p] = postFitStatistics(gprMdl_k);
            RSS = sum(r.^2, "all");
            BIC = n*log(RSS/n) + p*log(n);
            cv_BIC_scores(k) = BIC;
        end
        
        % Compute mean BIC score across folds
        mean_BIC_score = mean(cv_BIC_scores);
        
        % Store result in matrix
        BIC_scores(i,j) = mean_BIC_score;
        
        BIC_scores(i,j) = BIC;
    end
end

% Find combination with lowest BIC score
[min_BIC_score, min_idx] = min(BIC_scores(:));
[min_kernel_idx, min_basis_idx] = ind2sub(size(BIC_scores), min_idx);
best_kernel = kernel_options{min_kernel_idx};
best_basis = basis_options{min_basis_idx};

opt_kernel = best_kernel;
opt_basis = best_basis;
min_score = min_BIC_score;

end