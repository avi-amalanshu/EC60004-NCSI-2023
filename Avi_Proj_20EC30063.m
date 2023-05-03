%% Metadata
% Author       : Avi Amalanshu
% Inst.        : Dept of E&ECE, IIT Kharagpur
% Last Updated : May 05, 2023
% x-----------------------------------------------------------------------x
% |This code produces all the required figures for the term project of the|
% |course EC60004 (Neuronal Coding of Sensory Information), session Spring|
% |2022-2023.                                                             |
% x-----------------------------------------------------------------------x
clear;
%% Question 1: 
%  By using the AN model as informed in class to generate the following: Use a high spontaneous rate
%  auditory nerve fiber (ANF) with Best Frequency (BF) of 500 Hz kHz and 4 kHz and obtain their tuning
%  curves (response rates as a function of frequency) at 10 different intensities: -10 dB SPL to 80 dB
%  SPL (in steps of 10 dB). Use tone frequencies of 125 Hz to 16 kHz (a total of 7 octaves) with 8
%  frequencies in each octave (1/8th octave frequency difference). That is the tone frequencies will be
%  125*2.^[0:1/8:7] Hz. Use a duration of 200 ms for each tone and modulate the tones with onset and
%  offset ramps of 10 ms. Use 20 repetitions of each tone and obtain the average rates. Plot all the
%  tuning curves of each ANF in one figure (use a logarithmic frequency axis, Figure 1 and Figure 2).
%  Obtain the rate vs intensity function for BF tone of each ANF at the 10 intensities above and plot them
%  (Figure 3). What are the observations?

% ====================
%  Some global params
% ====================

BF = [500 4000];            % Best frequencies
freqs = 125*2.^(0:1/8:8);   % Tone frequencies' vector
intensities = -10:10:80;    % Intensities' vector
repetitions = 5;            % Repetitions of each tone

% <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
%                               Spectrograms
% <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

% ==========
% BF = 500Hz
% ==========

rates_500 = zeros(repetitions, length(intensities), length(freqs));         % This way, it is easier to average over tone
                                                                            % repetitions. We abstract 'nrep' in the global
                                                                            % parameter 'repetitions'.
for rep=1:repetitions
    fprintf("500Hz: Repetition #%d \n", rep);
    for i = 1:length(intensities)
        for f = 1:length(freqs)
            % model fiber parameters
            CF          = BF(1);    % CF in Hz;   
            cohc        = 1.0;      % normal ohc function
            cihc        = 1.0;      % normal ihc function
            fiberType   = 3;        % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
            implnt      = 0;        % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
            % stimulus parameters
            F0          = freqs(f); % stimulus frequency in Hz
            Fs          = 100*10^3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)
            T           = 200e-3;   % stimulus duration in seconds
            rt          = 10e-3;    % rise/fall time in seconds
            stimdb      = intensities(i);       % stimulus intensity in dB SPL
            % PSTH parameters
            nrep        = 4;        % number of stimulus repetitions (e.g., 50);
            psthbinwidth= 0.5e-3;   % binwidth in seconds;

            t = 0:1/Fs:T-1/Fs; % time vector
            mxpts = length(t);
            irpts = rt*Fs;
            
            pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t);
            pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
            pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

            vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
            [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 

            timeout = (1:length(psth))*1/Fs;
            psthbins = round(psthbinwidth*Fs); 
            psthtime = timeout(1:psthbins:end);
            pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep;
            Psth = pr/psthbinwidth;
            avg_Psth = sum(Psth)/length(Psth);
            rates_500(rep, i, f) = avg_Psth;
        end
    end
end
rates_500_avg = zeros(length(intensities), length(freqs));
for i = 1:length(intensities)
    for f = 1:length(freqs)
        rates_500_avg(i,f) = sum(rates_500(:,i,f))/repetitions;             % each nrep cycle is already averaged
    end
end

figure
    hold on
        for i = 1:length(intensities)
            y = rates_500_avg(i,:);
            x=log10(freqs);
            plot(x, y);
        end
    hold off
    title('BF = 500, -10dB to 80dB')
grid

% ===========
% BF = 4000Hz
% ===========

rates_4000 = zeros(repetitions, length(intensities), length(freqs));
for rep=1:repetitions
    fprintf("4000Hz: Repetition #%d \n", rep);
    for i = 1:length(intensities)
        for f = 1:length(freqs)
            % model fiber parameters
            CF          = BF(2);    % CF in Hz;   
            cohc        = 1.0;      % normal ohc function
            cihc        = 1.0;      % normal ihc function
            fiberType   = 3;        % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
            implnt      = 0;        % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
            % stimulus parameters
            F0          = freqs(f); % stimulus frequency in Hz
            Fs          = 100*10^3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)
            T           = 200e-3;   % stimulus duration in seconds
            rt          = 10e-3;    % rise/fall time in seconds
            stimdb      = intensities(i);       % stimulus intensity in dB SPL
            % PSTH parameters
            nrep        = 4;        % number of stimulus repetitions (e.g., 50);
            psthbinwidth= 0.5e-3;   % binwidth in seconds;

            t = 0:1/Fs:T-1/Fs; % time vector
            mxpts = length(t);
            irpts = rt*Fs;
            
            pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t);  
            pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
            pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

            vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
            [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 

            timeout = (1:length(psth))*1/Fs;
            psthbins = round(psthbinwidth*Fs);
            psthtime = timeout(1:psthbins:end);
            pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep;
            Psth = pr/psthbinwidth;
            avg_Psth = sum(Psth)/length(Psth);
            rates_4000(rep, i, f) = avg_Psth;
        end
    end
end

rates_4000_avg = zeros(length(intensities), length(freqs));
for i = 1:length(intensities)
    for f = 1:length(freqs)
        rates_4000_avg(i,f) = sum(rates_4000(:,i,f))/repetitions;
    end
end

figure
    hold on
        for i = 1:length(intensities)
            y = rates_4000_avg(i,:);
            x=log10(freqs);
            plot(x, y);
        end
    hold off
    title('BF = 4000Hz, -10dB to 80dB')
grid

% <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
%                           Rate vs Intensity
% <<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

% ==========
% BF = 500Hz
% ==========

rate_int_500 = zeros(length(intensities));
for i = 1:length(intensities)
    rate_int_500(i) = rates_500_avg(i,17);
end

% ===========
% BF = 4000Hz
% ===========

rate_int_4000 = zeros(length(intensities));
for i = 1:length(intensities)
    rate_int_4000(i) = rates_4000_avg(i,41);
end

figure
    hold on
    plot(intensities,rate_int_500);
    plot(intensities,rate_int_4000);
    title('rate vs intensity');
    legend('500Hz','4000Hz');
grid
%% Question 2: 
%  Now have a bank of ANFs starting with BF 125 Hz up to 8 kHz (a total of 6 octaves) with 16 ANFs
%  in each octave spaced 1/16th octaves apart (like frequencies presented in Part 1; ANF BFs
%  125*2.^[0:1/12:6] Hz; a total of 72 different BFs of ANFs).
%  Fixing sound level: Use a steady state portion of the speech sound wavfile provided („ah‟ part of
%  b“a”sketball). Use the wavread or audioread function to read it into MATLAB. Separate the “ah” out
%  from your speech signal waveform – by trial and hearing the segment. Use the root mean square
%  value of the segment to calculate its dB SPL level [re 20*10^(-6)]. Use this steady state sound level
%  and multiply the entire speech signal with appropriate factors to input in the ANFs (bank) for 3
%  different sound levels. Determine the 3 sound level as follows. Use the steady state portion and
%  modify it with onset and offset ramps as in Part 1 and find the rate responses to the vowel “ah” of a
%  500 Hz BF ANF at -20 to 80 dB SPL in 5 dB steps, plot (Figure 4) the rate intensity function (comment
%  by comparing it with the BF tone rate intensity function). Choose 3 sound levels one near (but above)
%  threshold, one in the dynamic range and one in the saturation level close to the end of the dynamic
%  range. After having determined the 3 sound levels generate the spike trains (50 repetitions each) of
%  each ANF in the bank (72 fibers) to the entire speech signal at the 3 sound levels.
%  Plot (Figure 5) the spectrogram of the speech signal with appropriate window size (25.6 ms hanning
%  windows maybe used with overlap of successive windows by 50%, that is, a resolution of 12.8 ms).
%  Now compare the spectrogram with the following: Represent the responses determined above from
%  each ANF as an average rate (number of spikes per unit time) as a function of time. Use windows of 4
%  ms, 8 ms, 16 ms, 32, ms 64 ms and 128 ms (with overlap between successive windows by 50%,
%  Figure 6A-F, 6 different window sizes). Plot the rate in an image in color with one axis as time (centre
%  of each successive window) and the other axis as BF of the ANFs: It is akin to a spectrogram
%  (cochleogram), only that now you have rate response instead of energy and BF instead of frequency.
%  Also in comparing the spectrogram with the above images do not forget that the spectrogram has a
%  linear frequency axis whereas the ANF BFs are spaced logarithmically axis.

intensities5 = -20:5:80;
[fivewo, Fs] = audioread('fivewo.wav');   %reading entire wav file
ah = fivewo(98250:113499);    % 'ah' sound
sound(ah, Fs);
rms_ah = rms(ah);
I_0 = 20*10^(-6);
db_SPL = 20*log10(rms_ah/I_0)
appr_fac = zeros(1,length(intensities5));

for i = 1:length(intensities5)
   appr_fac(1,i) = (10^(intensities5(i)/20)*(20 * (10^-6)))/rms_ah; % appropriate factors
end
ah = ah';
% model fiber parameters
CF          = 500;      % CF in Hz;   
cohc        = 1.0;      % normal ohc function
cihc        = 1.0;      % normal ihc function
fiberType   = 3;        % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt      = 0;        % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
F0          = CF;       % stimulus frequency in Hz
Fs          = 100*10^3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)
T           = length(ah)/Fs;  % stimulus duration in seconds
rt          = 5e-3;     % rise/fall time in seconds
% PSTH parameters
nrep        = 1;        % number of stimulus repetitions (e.g., 50);
psthbinwidth= 0.5e-3;   % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

repetitions = 50;
rate_int_ah = zeros(repetitions, length(intensities5));

for rep=1:repetitions
    for i = 1:length(intensities5)
        pin = ah*appr_fac(1,i);  
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

        vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 

        timeout = (1:length(psth))*1/Fs;
        psthbins = round(psthbinwidth*Fs);
        psthtime = timeout(1:psthbins:end);
        pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep;
        Psth = pr/psthbinwidth;
        Psth_avg = sum(Psth)/length(Psth);
        rate_int_ah(rep, i) = Psth_avg;
    end
end

rate_int_ah_avg = zeros(1, length(intensities5));
for i = 1:length(intensities5)
    rate_int_ah_avg(1, i) = sum(rate_int_ah(:,i))/repetitions;
end

figure
    hold on;
    plot(intensities,rate_int_500)
    plot(intensities5,rate_int_ah_avg)
    legend('500Hz','ah');
    title('Rate vs Intensity (BF = 500)')
    hold off;
grid

% <<<>>><<<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
%               Filtering the three different sound levels
% <<<>>><<<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

[fivewo, Fs] = audioread('fivewo.wav');
fivewo = fivewo';
rms_val = rms(fivewo);
I_0 = 20*10^(-6);
db_SPL = 20*log10(rms_val/I_0);
appr_fac = zeros(1,3);

I_1 = 0; I_2 = 40; I_3 = 76; 
lev = [I_1, I_2, I_3];
bank_freqs = 125*2.^[0:1/12:6];
bank_size = length(bank_freqs);

fwTime      = length(fivewo)/Fs;
s_tr_I_1    = zeros(bank_size,fwTime*Fs*2);
s_tr_I_2    = zeros(bank_size,fwTime*Fs*2);
s_tr_I_3    = zeros(bank_size,fwTime*Fs*2); % spike trains

for i = 1:3
   appr_fac(1,i) = (10^(lev(i)/20)*(20*(10^-6)))/rms_val;
end

% model fiber parameters
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
Fs = 100*10^3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = length(fivewo)/Fs;  % stimulus duration in seconds
% T = 50;
rt = 5e-3;   % rise/fall time in seconds
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;
% Feeding level 0dB
for b = 1:bank_size
        pin = fivewo*appr_fac(1,1);  
        CF = bank_freqs(b);
        vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt);
        s_tr_I_1(b,:) = psth;
end
% Feeding level 40dB
for b = 1:bank_size
        pin = fivewo*appr_fac(1,2);  
        CF = bank_freqs(b);
        vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt);
        s_tr_I_2(b,:) = psth;
end
% Feeding level 76dB
for b = 1:bank_size
        pin = fivewo*appr_fac(1,2);  
        CF = bank_freqs(b);
        vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt);
        s_tr_I_3(b,:) = psth;
end

% <<<>>><<<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>
%                           Fake Spectrogram
% <<<>>><<<>>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>><<>>

figure
    spectrogram(fivewo, hann(25.6e-3*Fs), 12.8e-3*Fs, 1:8000, Fs, 'yaxis');
    title('Audio spectrogram')
grid

s_tr_I_3 = squeeze(s_tr_I_3);

wind_sizes_ms = [4,8,16,32,64,128];
wind_sizes = wind_sizes_ms*1e-3*Fs;
half_wind_sizes = floor(wind_sizes)./2;

for w = 1:length(wind_sizes)
    interval = wind_sizes(w)/2:half_wind_sizes(w):length(fivewo)-wind_sizes(w)/2;
    response_rates = zeros(length(bank_freqs),length(interval));
    for f = 1:length(bank_freqs)
            for i = 1:length(interval)
                response_rates(f, i) = sum(s_tr_I_3(f, interval(i)-half_wind_sizes(w)+1:interval(i)+half_wind_sizes(w)))/wind_sizes(w);
            end
    end
    
    figure(6)
        subplot(2,3,w)
        [t, f] = meshgrid(interval, bank_freqs);
        surf(t, f, response_rates,'edgecolor','none');
        colorbar;
        xlim([0,1.5e5]);
        view(2);
        title(sprintf("window = %d ms\n", wind_sizes_ms(w)));

end

%% Question 3: 
%  Use the PSTHs from 50 repeats of the stimulus in every ANF, using a window size 0.1 ms or 100
%  microseconds. Consider the 12.8 ms long successive windows (50% overlap in successive windows)
%  and get the discrete Fourier Transform (use the fft function) of the PSTH. This is an indirect way of
%  looking at phase locking, that too relative amounts of locking to many different frequencies can be
%  observed simultaneously. Find the frequency to which a fiber locks the most, that is, find the peak in
%  the fft and its corresponding frequency which is the dominant frequency. Get the dominant frequency
%  in each successive window. Mark the frequency and time location on top of the spectrogram (say with
%  an asterisk). Do not use all the BFs of ANFs for this purpose. Use only BFs 1 octaves apart and only
%  up to 4 kHz (0.125, 0.25, 0.5, 1, 2 and 4 kHz). So there will be total of 6 fibers, use 6 colors of
%  asterisks and overlay them on the spectrogram at appropriate frequencies (dominant frequency) and
%  time (Figure 7). In a separate figure (Figure 8) do the same for another 5 fibers with BFs at 1 octave
%  intervals starting ½ octave above 125 Hz and ending ½ octave below 8 kHz (ie ~ 0.177 to 5.657 kHz).
%  Comment on your observations.

% =========================
% Even-Half Integer Octaves
% =========================

F_idx = [1,13,25,37,49,61];
cmap1 = hsv(11);
wind_size = 12.8e-3*Fs;
half_wind = floor(wind_size/2);
interval = wind_size/2:half_wind:length(fivewo)-wind_size/2;
max_pts = zeros(1, length(interval));

figure(7)
    spectrogram(fivewo, hann(12.8e-3*Fs), 6.4e-3*Fs, 1:8000, Fs, 'yaxis');
    view(3);
    hold on;
for f = 1:length(F_idx)
    for i = 1:length(interval)
        rate = s_tr_I_3(F_idx(f), (interval(i)-half_wind+1) : (interval(i)+half_wind));
        m = mean(rate);
        FFT = abs(fft(rate - m));
        [M, I] = max(squeeze(FFT(1:length(FFT)/2)));
        max_pts(i) = I*Fs/length(FFT);
    end

    %scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-50,[],cmap1(f,:),'filled', 'MarkerEdgeColor', 'k');
    scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-10,[],cmap2(f,:),'filled', 'MarkerEdgeColor', 'k');
    %set(gca,'Xscale','log');
    ylim([0 3]);
    hold on;
end
view(2)
grid
hold off;

% ========================
% Odd-Half Integer Octaves
% ========================

F_idx = [7,19,31,43,55,67];
cmap2 = hsv(11);
wind_size = 12.8e-3*Fs;
half_wind = floor(wind_size/2);
interval = wind_size/2:half_wind:length(fivewo)-wind_size/2;
max_pts = zeros(1, length(interval));

figure(8)
    spectrogram(fivewo, hann(12.8e-3*Fs), 6.4e-3*Fs, 1:8000, Fs, 'yaxis');
    view(3);
    hold on;
for f = 1:length(F_idx)
    for i = 1:length(interval)
        rate = s_tr_I_3(F_idx(f), (interval(i)-half_wind+1) : (interval(i)+half_wind));
        m = mean(rate);
        FFT = abs(fft(rate - m));
        [M, I] = max(squeeze(FFT(1:length(FFT)/2)));
        max_pts(i) = I*Fs/length(FFT);
    end

    %scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-50,[],cmap1(f,:),'filled', 'MarkerEdgeColor', 'k');
    %scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-10,[],cmap2(f,:),'filled',
    %'MarkerEdgeColor', 'k'); Absolutely 0 clue why I need to reuse cmap1.
    scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-10,[],cmap1(f,:),'filled', 'MarkerEdgeColor', 'k'); 
    %set(gca,'Xscale','log');
    ylim([0 3]);
    hold on;
end
view(2)
grid
hold off;
