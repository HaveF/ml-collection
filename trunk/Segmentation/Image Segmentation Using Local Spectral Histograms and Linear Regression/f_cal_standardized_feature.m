function X_standardized = f_cal_standardized_feature(X)

% X: p x n
m_val = mean(X, 2);
std_val = std(X, 0, 2);

X_standardized = X - repmat(m_val, 1, size(X, 2));
X_standardized = X_standardized ./ repmat(std_val, 1, size(X, 2));
