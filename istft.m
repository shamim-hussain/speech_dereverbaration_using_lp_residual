%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Inverse Short-Time Fourier Transform        %
%               with MATLAB Implementation             %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov        12/26/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, t] = istft(stft, wlen, hop, nfft, fs)

% function: [x, t] = istft(stft, wlen, hop, nfft, fs)
% stft - STFT matrix (only unique points, time across columns, freq across rows)
% wlen - length of the sinthesis Hamming window
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
% x - signal in the time domain
% t - time vector, s

% signal length estimation and preallocation
coln = size(stft, 2);
xlen = wlen + (coln-1)*hop;
x = zeros(1, xlen);

% form a periodic hamming window
win = hamming(wlen, 'periodic');

% initialize the signal time segment index
indx = 0;

% perform ISTFT (via IFFT and Weighted-OLA)
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    for col = 1:coln
        % extract FFT points
        X = stft(:, col);
        X = [X; conj(X(end:-1:2))];
        
        % IFFT
        xprim = real(ifft(X));
        xprim = xprim(1:wlen);
        
        % weighted-OLA
        x((indx+1):(indx+wlen)) = x((indx+1):(indx+wlen)) + (xprim.*win)';
        
        % update the index
        indx = indx + hop;
    end
else                                % even nfft includes Nyquist point
    for col = 1:coln
        % extract FFT points
        X = stft(:, col);
        X = [X; conj(X(end-1:-1:2))];
        
        % IFFT
        xprim = real(ifft(X));
        xprim = xprim(1:wlen);
        
        % weighted-OLA
        x((indx+1):(indx+wlen)) = x((indx+1):(indx+wlen)) + (xprim.*win)';
        
        % update the index
        indx = indx + hop;
    end
end

% scale the signal
W0 = sum(win.^2);                  
x = x.*hop/W0;                      

% generate time vector
t = (0:xlen-1)/fs;                 

end