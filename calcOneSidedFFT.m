% Function for calculating the frequency vector f and the single sided FFT
% Author Pieter Try, 27.03.2023
% Input: data, N, Fs, boolLog

function [f, P] = calcOneSidedFFT(data, N, Fs, boolLog)
%Aref = 1;
if isempty(N)
    N = size(data,1);
end
if N

    P_orig = fft(data, N,1);  %calculate FFT
    P = abs(P_orig) ./ N;  %Get Abs of complex value and divide by number of datapoints
    P = P(1:floor(N/2+1)); %get one sided Spectrum
    P(2:end-1) = 2* P(2:end-1);


    if boolLog
        P = 10 * log10(abs(P/max(P)));
    end
    f = Fs/2 .* linspace(0,1,length(P))';
end