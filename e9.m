clc;
close all;
N=20000;
SNR_limit=35;
SNR_db=-5:0.5:SNR_limit;
SNR=10.^(SNR_db/10);
u=rand(1,N);
m=floor(2*rand(1,N));
var=1; % sigma^2
nstd=sqrt(var);
y=(zeros(1,N));
Pe_BPSK_sim=(zeros(1,length(SNR)));
Pe_BFSK_sim=(zeros(1,length(SNR)));
Pe_DPSK_sim=(zeros(1,length(SNR)));
%4(a)generate rayleighranom variable
r=sqrt(-(2*var*log(u)));
figure(1);
hist(r,100);
title('rayleigh random variable histogram plot');
xlabel('random variable R');
ylabel('frequency')
a=[0:0.01:10];
R=(a/var).*exp(-(a.*a)/(2*var));
figure(2);
plot(a,R);
title('rayleigh PDF');
xlabel('random variable');
ylabel('probability');
legend('variance =1')
%BPSK simulation
Pe_BPSK_id=0.5*(1-sqrt((var*SNR)./(1+var*SNR)));
%BFSK simulation (coharent)
BFSK_id = 0.5*(1-sqrt(var*SNR./(2+(var*SNR))));
%DPSK simulation
Pe_DPSK_id=0.5./(1+var*SNR);

%Comparison of Error performance for AWGN and rayleigh fading channels
Pe_BPSK_NF=0.5*(erfc(sqrt(SNR)));
Pe_BFSK_NF=0.5*(erfc(sqrt(SNR/2)));
Pe_DPSK_NF=0.5*exp(-SNR); %non-coharent
figure(3);
semilogy(SNR_db,Pe_BPSK_id,'r.-', SNR_db,BFSK_id,'r*-',SNR_db,Pe_DPSK_id,'r--', SNR_db,Pe_BPSK_NF,'b.-',SNR_db,Pe_BFSK_NF,'b*-',SNR_db,Pe_DPSK_NF,'b--');
axis([-5 SNR_limit 0.000001 1]);
title('performance of BPSK,BFSK,DPSK')
xlabel('SNR(db)');
ylabel('probability of error');
legend('Pe of BPSK with fading','Pe of BFSK with fading','Pe of DPSK with fading','Pe of BPSK without fading','Pe of BFSK without fading','Pe of DPSK without fading');

%BFSK
clc;
%close all;
N=20000;
SNR_limit=35;
SNR_db=-5:0.5:SNR_limit;
SNR=10.^(SNR_db/10);
u=rand(1,N);
m=floor(2*rand(1,N));
var=1; % sigma^2
nstd=sqrt(var);
y=(zeros(1,N));
Pe_BPSK_sim=(zeros(1,length(SNR)));
Pe_BFSK_sim=(zeros(1,length(SNR)));
Pe_DPSK_sim=(zeros(1,length(SNR)));
%4(a)generate rayleighranom variable
r=sqrt(-(2*var*log(u)));
figure(4);
hist(r,100);
title('rayleigh random variable histogram plot');
xlabel('random variable R');
ylabel('frequency')
a=[0:0.01:10];
R=(a/var).*exp(-(a.*a)/(2*var));
figure(5);
plot(a,R);
title('rayleigh PDF');
xlabel('random variable');
ylabel('probability');
legend('variance =1')
%BPSK simulation
Pe_BPSK_id=0.5*(1-sqrt((var*SNR)./(1+var*SNR)));
%BFSK simulation (coharent)
BFSK_id = 0.5*(1-sqrt(var*SNR./(2+(var*SNR))));
%DPSK simulation
Pe_DPSK_id=0.5./(1+var*SNR);
%Comparison of Error performance for AWGN and rayleigh fading channels
Pe_BPSK_NF=0.5*(erfc(sqrt(SNR)));
Pe_BFSK_NF=0.5*(erfc(sqrt(SNR/2)));
Pe_DPSK_NF=0.5*exp(-SNR); %non-coharent
figure(6);
semilogy(SNR_db,Pe_BPSK_id,'r.-', SNR_db,BFSK_id,'r*-',SNR_db,Pe_DPSK_id,'r--', SNR_db,Pe_BPSK_NF,'b.-',SNR_db,Pe_BFSK_NF,'b*-',SNR_db,Pe_DPSK_NF,'b--');
axis([-5 SNR_limit 0.000001 1]);
title('performance of BPSK,BFSK,DPSK')
xlabel('SNR(db)');
ylabel('probability of error');
legend('Pe of BPSK with fading','Pe of BFSK with fading','Pe of DPSK with fading','Pe of BPSK without fading','Pe of BFSK without fading','Pe of DPSK without fading');

%BPSK
%%Performance of Binary Modulation in Rayleigh Fading Channel
%%transmission through a frequency nonselective channel
clc;
Eb= 1; % Energ per bit
EbNo_dB= 0:5:35; % vary the average SNR
No_over_2= Eb*10.^(-EbNo_dB/10);% Noise power
sigma= 1 ; % Rayleigh parameter
var=sigma^2;
BER= zeros(1,length(EbNo_dB));

% Calculation of error probabilit using Monte Carlo simulation:
for i = 1:length(EbNo_dB)
no_errors = 0;
no_bits = 0;
% Assumption: m = 0 (All zero codeword is transmitted):
while no_errors <= 10
u = rand;
% rand returns a single uniformly distributed
%random number in the interval (0,1)
alpha = sigma*sqrt(- 2*log(u)) ; % alpha is non-negative
%Rayleigh distrbuted with variance selected to be unity.
noise = sqrt(No_over_2(i))*randn; %randn gives Gaussian, with zero mean and unit variance
y = alpha*sqrt(Eb) + noise; %simulate the input to the detector
if y <= 0
y_d = 1;
else
y_d = 0;
end
no_bits = no_bits + 1 ;
no_errors = no_errors + y_d;
end
BER(i) = no_errors/no_bits ;%estimated error probability
end
% Calculation of error probabilit using the theoretical formula:
rho_b = Eb./No_over_2*var;
P2 = 1/2*(1-sqrt(rho_b./(1+rho_b))); %the theoretical value
% Plot the results:
semilogy(EbNo_dB,BER, '-* ' ,EbNo_dB, P2,'-o')
title('Montecarlosimualtion for Performance of BPSK signal');
xlabel('Average SNR/bit (dB)')
ylabel('Error Probability')
legend('Monte Carlo simulation','Theoretical value')