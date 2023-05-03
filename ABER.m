%%%%%%% ABER-RIS-SSK %%%%%%%%%
close all;
clc; clear; 
%%%
M=2;
N=128;
n_t=2;
n_r=1;

p_db = -10:5:30;  % db scale
p = 10.^(p_db/10); % linear scale

aber_denom = 1/(n_t*M*log2(n_t*M));  % constant term
aber = zeros(1, length(p_db));

e=eye(n_t);
%theta=(2*pi(m-1)/M)*eye(M);


for val = 1:length(p_db)  % index value 
D_ham = 0;
for m=1:M
    for i=1:n_t
        for m_hat=1:M 
            for i_hat=1:n_t 
                if i_hat~=i
                    D1=pdist([e(:,i)';e(:,i_hat)'],'hamming');
                else
                    D1=0;
                end
                theta=(2*pi*(m-1)/M)*eye(M);
                if m_hat~=m
                    D2=pdist([theta(:,m)';theta(:,m_hat)'],'hamming');
                else
                    D2=0;
                end
                D_ham=D1+D2+D_ham;
            end
        end
    end
end

gamma_tilde = (p(val)*N)/2;
mu_tilde=0.5*(1-sqrt(gamma_tilde/(1+gamma_tilde)));

%pep = 0;

for k=0:n_r-1
    sum=nchoosek(n_r-1-k,k)*(1-mu_tilde)^k; 
end
pep = (mu_tilde^n_r)*sum;

aber(val) = aber_denom*D_ham*pep;

end

semilogy(p_db, aber);
title('RIS-SSK-RPM ')
xlabel('SNR')
ylabel('ABER')
