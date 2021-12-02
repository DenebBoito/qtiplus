function out = create_new_example_dataset(bten, SNR_1, SNR_2)
% this function allows the creation of a new example dataset with the 
% specified level of noise
sigma1 = 1 / SNR_1;
sigma2 = 1 / SNR_2;
n = size(bten,1);

% create the final data container
out = zeros(64,64,6,n);

%% define parameters and signal for a non-central Wishart distribution for an isotropic mean diffusion tensor
% parameters
p = 2;   
D = eye(3) * 0.7e-9;
Sigma = D ./ (5*p);
Omega = 4 * p * Sigma;

% generate signal
signal = zeros(n,1);
for i = 1:n
    signal(i) = det(eye(3) + convert_1x6_to_3x3(bten(i,:)) * Sigma)^(-p) * (exp(-trace(convert_1x6_to_3x3(bten(i,:))*inv(eye(3) + Sigma * convert_1x6_to_3x3(bten(i,:))) * Omega)));
end

% insert value in out container
for i = 1:64
    for j = 1:64
        out(i,j,1,:) = signal;
        out(i,j,2,:) = add_rician_noise(signal,sigma1,sigma1);
        out(i,j,3,:) = add_rician_noise(signal,sigma2,sigma2);
    end
end

%% define parameters and signal for a non-central Wishart distribution for a anisotropic mean diffusion tensor
% parameters
p = 4;
D = eye(3) * 0.2e-9; D(9) = 1.3e-9; D(5) = 0.6e-9; 
Sigma = D ./ (5*p);
Omega = 4 * p * Sigma;

% generate signal
signal = zeros(n,1);
for i = 1:n
    signal(i) = det(eye(3) + convert_1x6_to_3x3(bten(i,:)) * Sigma)^(-p) * (exp(-trace(convert_1x6_to_3x3(bten(i,:))*inv(eye(3) + Sigma * convert_1x6_to_3x3(bten(i,:))) * Omega)));
end

% insert value in out container
for i = 1:64
    for j = 1:64
        out(i,j,4,:) = signal;
        out(i,j,5,:) = add_rician_noise(signal,sigma1,sigma1);
        out(i,j,6,:) = add_rician_noise(signal,sigma2,sigma2);
    end
end

end

function s_n = add_rician_noise(s,sigma1,sigma2)
% helper function to create a noisy signal
g_1 = normrnd(0,sigma1,size(s));
g_2 = normrnd(0,sigma2,size(s));
s_n = sqrt((s+g_1).^2 + g_2.^2);
end