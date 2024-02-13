function [y,e] = LMS_adaptivefiltering(ref_input,org_input,mu,nC)

d = org_input;              % original input signal
N = length(ref_input);      % jumlah sampel data
w = zeros(nC,1);            % vektor koefisien filter 
e = zeros(size(ref_input)); % vektor error
y = zeros(size(ref_input));
for n = nC:N
    x = ref_input(n:-1:n-nC+1); % reference input signal
    y(n) = w'*x;                   % filter prediction : mECG component of aECG signal
    e(n) = d(n) - y(n);            % evaluate error
    w = w + mu*e(n)*x;          % update filter weight
end

end