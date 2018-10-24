%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Short-Time Fourier Transform            %
%               with MATLAB Implementation             %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov        12/21/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stft, f, t] = stft(x, wlen, hop, nfft, fs)

% function: [stft, f, t] = stft(x, wlen, hop, nfft, fs)
% x - signal in the time domain
% wlen - length of the analysis Hamming window
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
% stft - STFT matrix (only unique points, time across columns, freq across rows)
% f - frequency vector, Hz
% t - time vector, s

% represent x as column-vector
x = x(:);

% length of the signal
xlen = length(x);

% form a periodic hamming window
win = hamming(wlen, 'periodic');

% stft matrix estimation and preallocation
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/hop);      % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix

% initialize the signal time segment index
indx = 0;

% perform STFT
for col = 1:coln
    % windowing
    xw = x(indx+1:indx+wlen).*win;
    
    % FFT
    X = fft(xw, nfft);
    
    % update the stft matrix
    stft(:, col) = X(1:rown);
    
    % update the index
    indx = indx + hop;
end

% calculate the time and frequency vectors
t = (wlen/2:hop:wlen/2+(coln-1)*hop)/fs;
f = (0:rown-1)*fs/nfft;

end