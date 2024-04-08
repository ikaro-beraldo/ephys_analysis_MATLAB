function Delta_band = wavelet(data, fs)
% Function to obtain the delta frequency band applied wavelet filter bank
%
% input: data
%        fs sample frequency of data
% output: delta frequency
%
% Natalia Espinosa 21/07/2023
% ------------------------------------------------------------------------

% Calculate the order of wavelet filter bank based on fs
counter = 1;
delta_freq = fs/2;
while (abs(delta_freq - 4) > 1)
   counter = counter + 1;
   delta_freq = delta_freq/2;  % ~ 4Hz
 
end

wavlevel = counter; % level of wavelet 1 to 24 for 'db' Daubechies wavelet family

wname = ['db' num2str(wavlevel)]; % wavelet type
        [LoD,HiD,LoR,HiR] = wfilters(wname); % Low and High pass filters Decomposition and Reconstruction
        [C,L]             = wavedec(data, wavlevel, LoD, HiD); % wavelet decomposition n levels
        Delta_band        = wrcoef('a', C, L, wname, wavlevel);  % DELTA

end