function PlotFrequency(data,timeStep,totalTime)

Fs = (1/timeStep)*0.1;   
t = timeStep:timeStep:totalTime;
x = data;
X = fft(x);
X=X(1:length(X)/2+1); %one-sided DFT
P = (abs(X)/length(x)).^2;     % Compute the mean-square power
P(2:end-1)=2*P(2:end-1); % Factor of two for one-sided estimate
% at all frequencies except zero and the Nyquist
Hmss=dspdata.msspectrum(P,'Fs',Fs,'spectrumtype','onesided');
freq = Hmss.frequencies;
power = Hmss.Data;
power = power/max(power);
figure(2)
plot(freq,power,'k','LineWidth',2); 

xlabel('Frequency');
ylabel('Normalized Power');

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',20); 

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',20);

set(gca,'XMinorTick','on');
set(gca,'FontSize',16);

% Fs = 1000; % Sampling frequency
% T = 1/Fs; % Sample time
% L = totalTime; % Length of signal
% % t = (timeStep:timeStep:L)*T; % Time vector
% y = data; % Sinusoids plus noise
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% % Plot single-sided amplitude spectrum.
% figure(2)
% plot(f,2*abs(Y(1:NFFT/2+1)));
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

