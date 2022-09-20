clc
clearvars
close all
set(0,'defaultTextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%
fs = 781e3;         % [Hz] sampling frequency for frequency response plot
fo = 35e3;          % [Hz] cutoff frequency
L  = 10e-6;         % [H]  inductance
C  = 2e-6;          % [F]  capacitance


% Noise transfer function
% First-order  (1-z^-1)
NTF1    = filt([1 -1],1)
[h1,w1] = freqz([1 -1],1,2048,fs);     
figure(1);
zplane([1 -1],1)

% Second-order (1-z^-1)^2
NTF2    = filt([1 -2 1],1); 
[h2,w2] = freqz([1 -2 1],1,2048,fs);
figure(2);
zplane([1 -2 1],1) 

% Third-order  (1-z^-1)^3
NTF3    = filt([1 -3 3 -1],1); 
[h3,w3] = freqz([1 -3 3 -1],1,2048,fs);
figure(3);
zplane([1 -3 3 -1],1) 

% Third-order modified (1-z^-1)(1- K*z^-1 + z^-2)
K = 2*cos(2*pi*fo/fs);
b1 = [1 -1];
b2 = [1 -K 1];
b = conv(b1,b2);
[hm,wm] = freqz(b,1,2048,fs);
zplane(b,1)

wb = linspace(1,fs*pi, 2048)';
hb = freqs(1,[L*C C 1], wb);

figure(7);
hold on
plot(w1,20*log10(abs(h1)),w2,20*log10(abs(h2)),w3,20*log10(abs(h3)), ...
     wm,20*log10(abs(hm)),wb/(2*pi),20*log10(abs(hb)), 'Linewidth', 1.5)
xlim([0 fs/2])
ylim([-80 30])

legend('$1^{st}$ Order','$2^{nd}$ Order','$3^{rd}$ Order',...
'$3^{rd}$ Order Mod','$2^{nd}$ order Buck Output Filter', 'location', 'southeast')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency response noise shaping transfer functions')
grid on

%% WITHOUT SHAPING
fs = 1000e3;        % simulation sampling frequency
Tsim = 5e-3;        % [s]  simulation time
t = (0:1/fs:Tsim-1/fs)';
d_sine_sim = 0.8*sin(2*pi*300*t);
n_bit = 4;
Q = 2/2^n_bit;      % quantization step
carrier = sawtooth(2*pi*50000*t,0.5);
car_q = quant(carrier, Q);

figure()
subplot(2,1,1)
hold on
stairs(t, car_q)
subplot(2,1,2)
hold on
plot(t,d_sine_sim)
stairs(t,car_q)
xlabel('Time (s)')

d_sine_q = quant(d_sine_sim, Q);

pwm0 = zeros(size(t));
for i = 1:1:length(t)
    if car_q(i) <= d_sine_q(i)  % if modulating greater than carrier pwm=1
        pwm0(i)=1;
    else
        pwm0(i)=0;
    end
end

figure();
hold on
plot(t,car_q,'k',t, d_sine_q,'r', t, pwm0);
ylim([-1.5 1.5])
title("$\Sigma\Delta$ DPWM simulation")
legend("Carrier","PWM input", "output $g$")
xlabel("Time (s)")

%% 1 order
x1 = zeros(size(t));
d1 = zeros(size(t));
q_e_sim1 = zeros(size(t));

for i = 1:length(t)
    if i == 1
        x1(i) = d_sine_sim(i);
        d1(i) = quant(x1(i), Q);
        q_e_sim1(i) = d1(i) - x1(i);
    else
        x1(i) = d_sine_sim(i) - q_e_sim1(i-1);
        d1(i) = quant(x1(i), Q);
        q_e_sim1(i) = d1(i) - x1(i);
    end
end

pwm1 = zeros(size(t));
for i = 1:1:length(t)
    if car_q(i) <= d1(i)
        pwm1(i)=1;
    else
        pwm1(i)=0;
    end
end
figure();
plot(t,car_q,'k', t, d1,'r', t, pwm1);
ylim([-1.5 1.5])
%xlim([0 1e-3])
title("$\Sigma\Delta$ DPWM simulation")
legend("Carrier","PWM input", "output $g$")
xlabel("Time (s)")
%% 2 order
x2 = zeros(size(t));
d2 = zeros(size(t));
q_e_sim2 = zeros(size(t));

for i = 1:length(t)
    if i == 1
        x2(i) = d_sine_sim(i);
        d2(i) = quant(x2(i), Q);
        q_e_sim2(i) = d2(i) - x2(i);
    elseif i == 2
        x2(i) = d_sine_sim(i) - 2*q_e_sim2(i-1);
        d2(i) = quant(x2(i), Q);
        q_e_sim2(i) = d2(i) - x2(i);    
    else
        x2(i) = d_sine_sim(i) - 2*q_e_sim2(i-1) + q_e_sim2(i-2);
        d2(i) = quant(x2(i), Q);
        q_e_sim2(i) = d2(i) - x2(i);
    end
end
% d is the modulating signal for the dpwm
pwm2 = zeros(size(t));
for i = 1:1:length(t)
    if car_q(i) <= d2(i)  % if modulating greater than carrier pwm=1
        pwm2(i)=1;
    else
        pwm2(i)=0;
    end
end

figure();
plot(t,car_q,'k', t, d2,'r', t, pwm2);
ylim([-1.5 1.5])
%xlim([0 1e-3])
title("$\Sigma\Delta$ DPWM simulation")
legend("Carrier","PWM input", "output $g$")
xlabel("Time (s)")
%% 3 order
x3 = zeros(size(t));
d3 = zeros(size(t));
q_e_sim3 = zeros(size(t));

for i = 1:length(t)
    if i == 1
        x3(i) = d_sine_sim(i);
        d3(i) = quant(x3(i), Q);
        q_e_sim3(i) = d3(i) - x3(i);
    elseif i == 2
        x3(i) = d_sine_sim(i) - 3*q_e_sim3(i-1);
        d3(i) = quant(x3(i), Q);
        q_e_sim3(i) = d3(i) - x3(i);
    elseif i == 3
        x3(i) = d_sine_sim(i) - 3*q_e_sim3(i-1) + 3*q_e_sim3(i-2);
        d3(i) = quant(x3(i), Q);
        q_e_sim3(i) = d3(i) - x3(i);
    else 
        x3(i) = d_sine_sim(i) - 3*q_e_sim3(i-1) + 3*q_e_sim3(i-2) - q_e_sim3(i-3);
        d3(i) = quant(x3(i), Q);
        q_e_sim3(i) = d3(i) - x3(i);        
    end    
end
% d is the modulating signal for the dpwm
pwm3 = zeros(size(t));
for i = 1:1:length(t)
    if car_q(i) <= d3(i)  % if modulating greater than carrier pwm=1
        pwm3(i)=1;
    else
        pwm3(i)=0;
    end
end

figure();
plot(t,car_q,'k', t, d3,'r', t, pwm3);
ylim([-1.5 1.5])
title("$\Sigma\Delta$ DPWM simulation")
legend("Carrier","PWM input", "output $g$")
xlabel("Time (s)")
%xlim([0 1e-3])

% low pass filtering example of sigma delta output
% -> the third order signal is the one with less noise 
figure();
plot(t, lowpass(d3,fo,fs),'r'); % d_sine_q   d1  d2
ylim([-1.5 1.5])
title("$\Sigma\Delta$ DPWM simulation")
xlabel("Time (s)")
%xlim([0 1e-3])

%% QUANTIZATION NOISE
figure()
hold on
histogram(q_e_sim3)  % uniform distribution q_e_sim1,q_e_sim2

% linear model
d1_ntf = filter([1 -1],1, q_e_sim1);
d2_ntf = filter([1 -2 1],1, q_e_sim2);
d3_ntf = filter([1 -3 3 -1],1, q_e_sim3);
[pxx1, f1] = periodogram(d1_ntf+d_sine_sim,rectwin(length(d1_ntf+d_sine_sim)),[],fs);
[pxx2, f2] = periodogram(d2_ntf+d_sine_sim,rectwin(length(d2_ntf+d_sine_sim)),[],fs);
[pxx3, f3] = periodogram(d3_ntf+d_sine_sim,rectwin(length(d3_ntf+d_sine_sim)),[],fs);

% simulation
[pxxd1, fd1] = periodogram(d1,rectwin(length(d1)),[],fs);
[pxxd2, fd2] = periodogram(d2,rectwin(length(d2)),[],fs);
[pxxd3, fd3] = periodogram(d3,rectwin(length(d3)),[],fs);

figure()
hold on
plot(f1,10*log10(pxx1),f2,10*log10(pxx2),f3,10*log10(pxx3))
grid on
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('$1^{st}$ Order','$2^{nd}$ Order','$3^{rd}$ Order', 'location', 'best')
title('PSD estimate of the linear model $d = d_h + NTF(z)q_e$')
xlim([0 fs/2])
ylim([-140 -50])

figure()
hold on
plot(fd1,10*log10(pxxd1),fd2,10*log10(pxxd2),fd3,10*log10(pxxd3))
grid on
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('$1^{st}$ Order','$2^{nd}$ Order','$3^{rd}$ Order', 'location', 'best')
title('PSD estimate of $\Sigma\Delta$ simulation output')
xlim([0 fs/2])
ylim([-140 -50])

%% PSPECTRUM
[pspect0, f0] = pspectrum(pwm0,fs); %d_sine_q
[pspect1, f1] = pspectrum(pwm1,fs); %d1
[pspect2, f2] = pspectrum(pwm2,fs); %d2
[pspect3, f3] = pspectrum(pwm3,fs); %d3

figure()
hold on
plot(f1, 10*log10(pspect1))
plot(f2, 10*log10(pspect2))
plot(f3, 10*log10(pspect3))
plot(f0, 10*log10(pspect0))
legend('$1^{nd}$ Order','$2^{rd}$ Order',...
'$3^{rd}$ Order','Without shaping', 'location', 'best')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Frequency content of the output signal $g$')
xlim([0 fs/2])

%%
% [pxx, f] = periodogram(car_q,hanning(length(car_q)),length(car_q)/5,fs,'psd');
% sinad(pxx,f,'psd')
% thd(pxx,f,'psd')
% 
% [pxx, f] = pwelch(car_q,[],[],[],fs);
% thd(pxx,f,'psd')
% sinad(pxx,f,'psd')

%% SINAD and THD dpwm output
[pxx,f] = periodogram(pwm3,hanning(length(pwm0)),length(pwm0)/5,'power');
rbw = enbw(hanning(length(pwm0)/5)); % pwm1 pwm2 pwm3
sinad(pxx,f,rbw,'power')

[pxx,f] = periodogram(pwm2,hanning(length(pwm0)),length(pwm0)/10,fs,'power'); 
rbw = enbw(hanning(length(pwm0)/10),fs); % pwm1 pwm2 pwm3
thd(pxx,f,rbw,'power')

