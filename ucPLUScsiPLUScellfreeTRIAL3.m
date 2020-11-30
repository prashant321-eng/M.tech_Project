%All the AP's are located within a 1KM^2 area randomly with coordinates x,y
clear;
M=100; K=10;
rng(1);
x=rand(1,M)*1000;y=rand(1,M)*1000;
subplot(3,1,1)
scatter(x,y)
hold on;
%All the MS's are located within a 1KM^2 area randomly with coordinates x,y
rng(1);
a=rand(1,K)*1000;b=rand(1,K)*1000;
scatter(a,b)
hold off;
m_final = [x;y]%co-ordinates of AP's
k_final = [a;b]%co-ordinates of MS's
for p=1:M
    for q=1:K
        dd_km(p,q)=sqrt((x(p)-a(q))^2 + (y(p)-b(q))^2);
    end
end
d_km = dd_km';%distance between MS's and AP's
%for calculating g_km ,pathloss and showding correlation cofficient needed
%path loss
d1=50; d0=10;f=1900;h_AP=15;h_MS=1.65;
L=46.3+33.9*log10(f)-13.82*log10(h_AP)-(1.11*log10(f)-0.7)*h_MS+1.56*log10(f)-0.8;
for p=1:K
    for q=1:M
        if d_km(p,q)>d1
            PL_km(p,q) = -L-35*log10(0.001*d_km(p,q))
        elseif d_km(p,q)>d0&d_km(p,q)<=d1
            PL_km(p,q) = -L-10*log10(((0.001*d1)^1.5)*((0.001*d_km(p,q))^2))
        elseif d_km(p,q)<d0
            PL_km(p,q) = -L-10*log10(((0.001*d1)^1.5)*((0.001*d0)^2))
        end
    end
end
%shadowing correlation coefficient calculation - doubt how to calculate?
% for p=1:M
%     for  q=1:M
%         if p==q
%             d_ap_mM(p,q)=0;
%         else 
%                 d_ap_mM(p,q)= sqrt((x(p)-x(q))^2 + (y(p)-y(q))^2);
%         end
%     end
% end
% d_decorr = 100;
% for p=1:M
%     for  q=1:M
%         if p==q
%             real_a_m(p,q)=1;%this is actually a_m but have doubt
%         else 
%                 real_a_m(p,q) = (1/(2^(d_ap_mM(p,q)/d_decorr)));
%         end
%     end
% end
% for p=1:K
%     for  q=1:K
%         if p==q
%             d_ms_kK(p,q)=0;
%         else 
%                 d_ms_kK(p,q)= sqrt((x(p)-x(q))^2 + (y(p)-y(q))^2);
%         end
%     end
% end
% for p=1:K
%     for  q=1:K
%         if p==q
%             real_b_k(p,q)=1;%this is actually b_k but have doubt
%         else 
%                 real_b_k(p,q)= (1/(2^(d_ms_kK(p,q)/d_decorr)));
%         end
%     end
% end
%aa_m = randn(1,M);a_m=repelem(aa_m,K,1);bb_m=randn(K,1);b_m = repelem(bb_m,1,M);
%shadowing correlation coefficient
z_km = sqrt(0.5).*randn(K,M)+sqrt(0.5).*randn(K,M);sigma_sh=8;
A_km = (10.^(PL_km/10)); B_km = (10.^((sigma_sh.*z_km)/10))
beta_km = A_km .* B_km ;%beta for the calculation of g_km
h_km = abs(sqrt(0.5)*randn(K,M)+i*sqrt(0.5)*randn(K,M));
g_km = sqrt(beta_km).*h_km;%channel gains
t_p = K;%pilot 
for p=1:t_p
    for q=1:K
        if p==q
            pilot_matrix(p,q)=1;
        else
            pilot_matrix(p,q)=0;
        end
    end
end
%roh_k = repelem(0.1,1,K);
%w_m  = zeros(100,10);
%w1 = sqrt(2*(10^-13.4))*randn(1,M); w2 = sqrt(2*(10^-13.4))*randn(1,M);
%ww_m = abs(w1+i*w2);
%for p=1:M
    %for q=1:K
        %if p==q
            %w_m(p,q) = ww_m(p);
        %end
    %end
%end
noise_variance = 2*(10^-13.4);%how to take into account the noise figure?
w_m  = abs(sqrt(noise_variance/2)*(randn(K,M)+i*randn(K,M)));
roh_k = 0.1;%uplink power for channel estimation
y_m = (K*sqrt(roh_k)).*(pilot_matrix*g_km) + w_m;%output matrix
estimate_pm = (1/sqrt(roh_k))*(pilot_matrix'*y_m);%pilot matched channel estimate
g_hat_pm = ((abs(estimate_pm)).^2);%square of the channel pilot matched estimate
sorted_g_hat_pm = fliplr(sort(g_hat_pm',2));%sorting the pm matrix
%estimate_LMMSE = ((sqrt(0.1)*(beta_km.*(pilot_matrix'*y_m)))/noise_variance);%channel estimate LMMSE
%g_hat_LMMSE = (abs(estimate_LMMSE)).^2;%square of channel estimate LMMSE
%sorted_g_hat_LMMSE = fliplr(sort(g_hat_LMMSE,2));%sorted LMMSE
n=2;%no of MS's to be served as strongest channels
%just do for pilot matched from here
AAA = sorted_g_hat_pm; BBB = g_hat_pm';
C = zeros(M*n,1);D = zeros(M*n,1);t =1;u=1;
%for every AP's(in C Matrix), the user with strongest channel (in D Matrix) 
for p=1:M
    for q = 1:n
        for r=1:M
            for s=1:K
                if AAA(p,q)==BBB(r,s)
                    C(t)=r;D(u)=s;
                    t=t+1;u=u+1;
                end
            end
        end
    end
end

for p=1:M
    sum=0;
    for q=0:n-1
        deno_n_km(p)=sum+(abs(estimate_pm(D(p+(p-1)+q),p))^2);
    end
end
roh_m_dl = 0.2 
for p=1:K
    for q =1:M
        n_km(p,q) = (roh_m_dl/deno_n_km(q));% n_km doesnt depend on k rather it depends on the set of k for every AP's
    end
end
E = [C D];
%this below is for the set M(k)where 1st colum is MS's and second is AP's
F = fliplr(sortrows(E,2));
q=1;
for r=1:K
    sum=0;
    while r==F(q,1) && q<M*n
        NUM(r)=sum+(abs(sqrt(n_km(r,F(q,2)))*g_km(r,F(q,2))*(estimate_pm(r,F(q,2))))^2);
        q=q+1;
    end
end
sum=0;
%calculating denominator for rate
for r=1:K
    sum_deno=sum+NUM(r);
end
sigma_zsquare = 2*(10^-13.4);W=20*(10^6);
for r = 1:K
    for s = 1:K
        if s~=r
            sum_deno(r) = sum+NUM(r);
        end
    end
    SINR_DL(r) = (NUM(r)/(sigma_zsquare + sum_deno(r)));
    Rate(r)= W*log2(1+SINR_DL(r));            
end
p = min(Rate):0.1*(10^7):max(Rate);
q = expcdf(p,mean(Rate));
subplot(3,1,3);
p1 = plot(p,q,'b');
hold on;
%xticks(0:20*(10^6):140*(10^6));
xlim([0,140*10^6]);
set(gca,'XTick',[0 : 20*10^6 : 140*10^6]);
%xlim([0 140*(10^6)]);


%USER CENTRIC PILOT MATCHED DONE ABOVE
%CELL FREE PILOT MTACHED BELOW




n=K;%no of MS's to be served as strongest channels
%just do for pilot matched from here
AAA = sorted_g_hat_pm; BBB = g_hat_pm';
C = zeros(M*n,1);D = zeros(M*n,1);t =1;u=1;
%for every AP's(in C Matrix), the user with strongest channel (in D Matrix) 
for p=1:M
    for q = 1:n
        for r=1:M
            for s=1:K
                if AAA(p,q)==BBB(r,s)
                    C(t)=r;D(u)=s;
                    t=t+1;u=u+1;
                end
            end
        end
    end
end

for p=1:M
    sum=0;
    for q=0:n-1
        deno_n_km(p)=sum+(abs(estimate_pm(D(p+(p-1)+q),p))^2);
    end
end
roh_m_dl = 0.2 
for p=1:K
    for q =1:M
        n_km(p,q) = (roh_m_dl/deno_n_km(q));% n_km doesnt depend on k rather it depends on the set of k for every AP's
    end
end
E = [C D];
%this below is for the set M(k)where 1st colum is MS's and second is AP's
F = fliplr(sortrows(E,2));
q=1;
for r=1:K
    sum=0;
    while r==F(q,1) && q<M*n
        NUM(r)=sum+(abs(sqrt(n_km(r,F(q,2)))*g_km(r,F(q,2))*(estimate_pm(r,F(q,2))))^2);
        q=q+1;
    end
end
sum=0;
%calculating denominator for rate
for r=1:K
    sum_deno=sum+NUM(r);
end
sigma_zsquare = 2*(10^-13.4);W=20*(10^6);
for r = 1:K
    for s = 1:K
        if s~=r
            sum_deno(r) = sum+NUM(r);
        end
    end
    SINR_DL(r) = (NUM(r)/(sigma_zsquare + sum_deno(r)));
    Rate(r)= W*log2(1+SINR_DL(r));            
end
p = min(Rate):0.1*(10^7):max(Rate);
q = expcdf(p,mean(Rate));
p2 = plot(p,q,'r');
hold on;
%xticks(0:20*(10^6):140*(10^6));
xlim([0,140*10^6]);
set(gca,'XTick',[0 : 20*10^6 : 140*10^6]);
%xlim([0 140*(10^6)]);
yticks(0:0.1:1);axis([0 4e7 0 1]);



%USER CENTRIC PERFECT CSI




g_km_square = ((abs(g_km)).^2);%square of the channel pilot matched estimate
%sorted_g_hat_pm = fliplr(sort(g_hat_pm',2));%sorting the pm matrix
sorted_g_km_square = fliplr(sort(g_km_square',2));
%estimate_LMMSE = ((sqrt(0.1)*(beta_km.*(pilot_matrix'*y_m)))/noise_variance);%channel estimate LMMSE
%g_hat_LMMSE = (abs(estimate_LMMSE)).^2;%square of channel estimate LMMSE
%sorted_g_hat_LMMSE = fliplr(sort(g_hat_LMMSE,2));%sorted LMMSE
n=2;%no of MS's to be served as strongest channels
%just do for pilot matched from here
AAA = sorted_g_km_square; BBB = g_km_square';
C = zeros(M*n,1);D = zeros(M*n,1);t =1;u=1;
%for every AP's(in C Matrix), the user with strongest channel (in D Matrix) 
for p=1:M
    for q = 1:n
        for r=1:M
            for s=1:K
                if AAA(p,q)==BBB(r,s)
                    C(t)=r;D(u)=s;
                    t=t+1;u=u+1;
                end
            end
        end
    end
end

for p=1:M
    sum=0;
    for q=0:n-1
        deno_n_km(p)=sum+(abs(g_km(D(p+(p-1)+q),p))^2);
    end
end
roh_m_dl = 0.2 
for p=1:K
    for q =1:M
        n_km(p,q) = (roh_m_dl/deno_n_km(q));% n_km doesnt depend on k rather it depends on the set of k for every AP's
    end
end
E = [C D];
%this below is for the set M(k)where 1st colum is MS's and second is AP's
F = fliplr(sortrows(E,2));
q=1;
for r=1:K
    sum=0;
    while r==F(q,1) && q<M*n
        NUM(r)=sum+(abs(sqrt(n_km(r,F(q,2)))*g_km(r,F(q,2))*(g_km(r,F(q,2))))^2);
        q=q+1;
    end
end
sum=0;
%calculating denominator for rate
for r=1:K
    sum_deno=sum+NUM(r);
end
sigma_zsquare = 2*(10^-13.4);W=20*(10^6);
for r = 1:K
    for s = 1:K
        if s~=r
            sum_deno(r) = sum+NUM(r);
        end
    end
    SINR_DL(r) = (NUM(r)/(sigma_zsquare + sum_deno(r)));
    Rate(r)= W*log2(1+SINR_DL(r));            
end
p = min(Rate):0.1*(10^7):max(Rate);
q = expcdf(p,mean(Rate));
p3 = plot(p,q,'g');
hold on;
%xticks(0:20*(10^6):140*(10^6));
xlim([0,140*10^6]);
set(gca,'XTick',[0 : 20*10^6 : 140*10^6]);
%xlim([0 140*(10^6)]);
yticks(0:0.1:1);axis([0 4e7 0 1]);



%CELL FREE PERFECT CSI




g_km_square = ((abs(g_km)).^2);%square of the channel pilot matched estimate
%sorted_g_hat_pm = fliplr(sort(g_hat_pm',2));%sorting the pm matrix
sorted_g_km_square = fliplr(sort(g_km_square',2));
%estimate_LMMSE = ((sqrt(0.1)*(beta_km.*(pilot_matrix'*y_m)))/noise_variance);%channel estimate LMMSE
%g_hat_LMMSE = (abs(estimate_LMMSE)).^2;%square of channel estimate LMMSE
%sorted_g_hat_LMMSE = fliplr(sort(g_hat_LMMSE,2));%sorted LMMSE
n=K;%no of MS's to be served as strongest channels
%just do for pilot matched from here
AAA = sorted_g_km_square; BBB = g_km_square';
C = zeros(M*n,1);D = zeros(M*n,1);t =1;u=1;
%for every AP's(in C Matrix), the user with strongest channel (in D Matrix) 
for p=1:M
    for q = 1:n
        for r=1:M
            for s=1:K
                if AAA(p,q)==BBB(r,s)
                    C(t)=r;D(u)=s;
                    t=t+1;u=u+1;
                end
            end
        end
    end
end

for p=1:M
    sum=0;
    for q=0:n-1
        deno_n_km(p)=sum+(abs(g_km(D(p+(p-1)+q),p))^2);
    end
end
roh_m_dl = 0.2 
for p=1:K
    for q =1:M
        n_km(p,q) = (roh_m_dl/deno_n_km(q));% n_km doesnt depend on k rather it depends on the set of k for every AP's
    end
end
E = [C D];
%this below is for the set M(k)where 1st colum is MS's and second is AP's
F = fliplr(sortrows(E,2));
q=1;
for r=1:K
    sum=0;
    while r==F(q,1) && q<M*n
        NUM(r)=sum+(abs(sqrt(n_km(r,F(q,2)))*g_km(r,F(q,2))*(g_km(r,F(q,2))))^2);
        q=q+1;
    end
end
sum=0;
%calculating denominator for rate
for r=1:K
    sum_deno=sum+NUM(r);
end
sigma_zsquare = 2*(10^-13.4);W=20*(10^6);
for r = 1:K
    for s = 1:K
        if s~=r
            sum_deno(r) = sum+NUM(r);
        end
    end
    SINR_DL(r) = (NUM(r)/(sigma_zsquare + sum_deno(r)));
    Rate(r)= W*log2(1+SINR_DL(r));            
end
p = min(Rate):0.1*(10^7):max(Rate);
q = expcdf(p,mean(Rate));
p4 = plot(p,q,'c');
hold on;
%xticks(0:20*(10^6):140*(10^6));
xlim([0,140*10^6]);
set(gca,'XTick',[0 : 20*10^6 : 140*10^6]);
%xlim([0 140*(10^6)]);
yticks(0:0.1:1);axis([0 4e7 0 1]);
h = [p1;p2;p3;p4];
legend(h,'uc pm','cf pm','uc perfect csi','cf perfect csi');