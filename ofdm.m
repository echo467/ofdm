clc;
clear;
%% Part 1
% Calculate  parameters
% Estimating the RMS delay spread of the channel
pvec = [-2,-4,-9,-18]; % multi-path channel attenuation coefficient dB
tvec = [100,300,500,1300]; %multi-path channel delay 1e-9s
pvec_W = 10.^(pvec/10); % convert dB to mW

% Calculate the mean excess delay
t_delay = (sum(pvec_W.*tvec)/sum(pvec_W)); % unit: ns
t_delay_squared = (sum(pvec_W.*(tvec.^2))/sum(pvec_W)); % delay squared
%  Calcualte the RMS delay spread
delay_spread  = sqrt(t_delay_squared-t_delay^2)*10^-9; % unit: s

%  Estimating the coherence bandwidth of the channel 
% Assuming 0.5 frequency correlation
Bc = 1/(5*delay_spread);

%  Designing an OFDM System
% Initialise the design requirements
Rb = 60*10^6; % Bit rate (60Mbps)
BW = 20*10^6; % Channel bandwidth (20MHz)

% Assumptions
dirac = 0.1; % proportion of gaurd interval
Rc = 1; % code rate
beta = 0.25; 

% Using the modulating index (u) to find the no. of sub carriers (N)
% NOTE: u = log2(M)
u = (Rb/BW)*(1+beta); % raised-cosine pulse shaping
u =ceil(u); % add one for u has to be greater than (Rb/BW)*(1+beta)
N = 2^8; % N was chosen as 2^8 = 256 as it was closest base2 value

% Using delta_f to find the symbol duration
delta_f = BW/N;
Ts = 1/delta_f; % symbol period

% Guard interval Tg > 10*delay_spread hence set Tg to 10*delay_spread
% rounded up to the nearest 7th decimal place since dealing with micro
% seconds
Tg = round(3*delay_spread,7);
T_total = Ts + Tg; % symbol period + guard 

% Evaluating the data rate of the OFDM system designed 
R = (N*u)/(T_total); 

% Evaluating the BW of the system
BW_system = R/u; % 17.181 MHz (< 20MHz channel requirement)

% Bandwidth of Single Carrier system
T_single = 1/BW;
BW_single = BW/u; % BW_single = 2.8571MHz 

M = 2^(u); % alphabet size
k = u;  % bits/Symbol
Nsyms = N; % number of modulated symbols
maxNumBits = 1e5;

% Calculate OFDM Parameters
ofdmBw = BW;
deltaF = ofdmBw/Nsyms; % OFDM sub carrier separation
Tsamp = Ts/Nsyms; % sample time
Ncp = floor(Tg/Tsamp)+1;  % number of samples for Guard interval
EbNoVec = (0:30)'; % EbNo vector (0 - 30 dB)
% Set the errorRate array and BER array
berVec_flat = zeros(length(EbNoVec),1);
theoryBER_flat = berfading(EbNoVec, 'qam',M,1);

% calculates the modulator scale factor for M-QAM
input = 0:M-1;
const = qammod(input,M);
scale = modnorm(const,'avpow',1);

% Calculate multi-path channel parameters
fad_pdB = pvec; % multipath delay vector in dB
fad_p = 10.^(fad_pdB/10);   % multipath delay vector in W
alpha = 1/sum(fad_p);   % estimated normalisation factor 
sigma_tau = alpha.*fad_p;   % variance of each path 
delayT = tvec*10^-9;    % Time delay with correct unit (nano-seconds)
delaySample = floor(delayT/Tsamp)+1;    % number of samples for each delay 

%% Part 2 
%design OFDM
% Use M-QAM modulation. 
hMod = comm.RectangularQAMModulator('ModulationOrder',M);
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',M);
EbNoVec = (0:30)'; % EbNo vector (0 - 30 dB)
snrVec = EbNoVec + 10*log10(k); % Find SNR from given EbNo

% Set the errorRate array and BER array
errorRate = comm.ErrorRate('ResetInputPort',true);
berVec_multifade = zeros(length(EbNoVec),1);
berVec_multifade1 = zeros(length(EbNoVec),1);
errorStats = zeros(1,3);

% iterate through each EbNo vector and calculate average BER.
% the actual BER will be calculated using biterr() function but errorStats
% will be used for controlling the while loop to ensure the data captures
% sufficient rounds of run.
for n=1:length(EbNoVec)
    % These two variables are used for averaging the BER over 
    % number of test runs in the while loop
    berSUM = 0;
    berSUM1 = 0;
    count = 0;
    
    % The while loop ensures we capture good BER result with high enough
    % max number of bits
    while errorStats(3) <= maxNumBits 
        dataIn = randi([0 M-1],Nsyms,1);  % Message signal in symbol forms
        dataIn_bin = reshape(de2bi(dataIn),k*Nsyms,1);  % message signal in binary form
        y = scale*step(hMod,dataIn);  % modulate the data using hMod object
        OFDM = ifft(y)*sqrt(Nsyms); % Generate the OFDM symbol by taking IFFT of the modulated symbols
        
        OFDM_CP =  [OFDM(Nsyms-Ncp+1:Nsyms); OFDM];  % Append cyclic prefix to the beginning of the symbols
        OFDM_CP1 =  [zeros(Ncp,1); OFDM]; % Zero filling prefix
        
        % Genrate a mutli-path fading channel
        fad_multipath = zeros(Nsyms,1); % create empty vector with length of total symbols
        fad_multipath(delaySample) = sqrt(sigma_tau/2)'.*(randn(length(fad_pdB),1)+1i*randn(length(fad_pdB),1)); % Generate Rayleigh fading channel   
        
        % Transmit signal through a fading + AWGN channel.
        OFDM_fad = filter(fad_multipath,1,OFDM_CP); % apply multi-fading
        ynoisy = awgn(OFDM_fad(1:length(OFDM_CP)),snrVec(n));   % add AWGN 
        OFDM_fad1 = filter(fad_multipath,1,OFDM_CP1); % apply multi-fading
        ynoisy1 = awgn(OFDM_fad1(1:length(OFDM_CP1)),snrVec(n));   % add AWGN 
        
        OFDM_Rx=ynoisy(Ncp+1:(length(OFDM_CP))); % remove cyclic prefix
        OFDM_Rx1=ynoisy1(Ncp+1:(length(OFDM_CP1))); % remove cyclic prefix
        
        OFDM_decode=fft(OFDM_Rx,Nsyms)/sqrt(Nsyms); % decode the message by fft 
        OFDM_decode1=fft(OFDM_Rx1,Nsyms)/sqrt(Nsyms); % decode the message by fft
        
        fadF = fft(fad_multipath,Nsyms);    % Single tap euqalisation to remove the channel effect.
        yeq = OFDM_decode./fadF;
        yeq1 = OFDM_decode1./fadF;
         
        dataOut = step(hDemod,yeq/scale);   % Demodulate to recover the message.  
        dataOut1 = step(hDemod,yeq1/scale);   % Demodulate to recover the message.  
        
        dataOut_bin = reshape(de2bi(dataOut),k*Nsyms,1);    % received message in binary form 
        dataOut_bin1 = reshape(de2bi(dataOut1),k*Nsyms,1);    % received message in binary form 
         
        [error1, ratio1] = biterr(dataIn_bin,dataOut_bin);  % calculate current run's BER
        [error2, ratio2] = biterr(dataIn_bin,dataOut_bin1);  % calculate current run's BER
        
        berSUM = berSUM + ratio1;     % sum of the BER for each iteration. will be used to calculate average BER.
        berSUM1 = berSUM1 + ratio2; 
        
        count = count + 1;  % counter for keeping track of total number of iteration
        errorStats = errorRate(dataIn,dataOut,0);     % check bits used doesnt exceed maxNumBits. If yes, exit while loop
    end % end while loop 
    
    errorStats = errorRate(dataIn,dataOut,1);         % Reset the error rate calculator
    berVec_multifade(n) = berSUM/count;   % calculate the average BER for each EbNo
    berVec_multifade1(n) = berSUM1/count;   % calculate the average BER for each EbNo
end % end for loop

% plot flat fading channel and multifading channel and observe the difference 
figure(1)
semilogy(EbNoVec,berVec_multifade,'b-s');
hold on;
semilogy(EbNoVec,berVec_multifade1,'r-o');
hold on;

legend('循环前缀','补零前缀');
hold on;
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
title('多(4)径瑞利衰落信道下误比特率曲线');
grid on
hold off


