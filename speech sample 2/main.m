clc; clear all; close all

fs=8e3;
Lsig=20;
Limp=256e-3;
Nsig=Lsig*fs;
Nimp=Limp*fs;
nsig=1:Nsig;
tsig=nsig/fs;
vwind=1:fs*5;

% load data.mat x
[x,fs]=audioread('male.wav');
load rir.mat h

h=h(1:Nimp).';
x=x(1:Nsig);
% x=awgn(x, 20, 'measured');

x=x/max(abs(x));
y=filter(h, 1, x);
y=y/max(abs(y));
% y=awgn(y, 20, 'measured');

xr=lpres(x);
yr=lpres(y);

figure
subplot(221);plot(tsig(vwind),x(vwind))
xlabel('time(s)');title('Clean speech, x(t)')
subplot(222);plot(tsig(vwind),y(vwind))
xlabel('time(s)');title('Reverberant speech, y(t)')
subplot(223);plot(tsig(vwind),xr(vwind))
xlabel('time(s)');title('LP residula of clean speech, x_r(t)')
subplot(224);plot(tsig(vwind),yr(vwind))
xlabel('time(s)');title('LP residual of reverberant speech, y_r(t)')

figure
subplot(211);plot(tsig(vwind),xr(vwind))
xlabel('time(s)');title('LP residula of clean speech, x_r(t)')
subplot(212);plot(tsig(vwind),yr(vwind))
xlabel('time(s)');title('LP residual of reverberant speech, y_r(t)')

% figure
% subplot(311);plot(xr)
% subplot(312);plot(yr)
% subplot(313);plot(filter(h,1,xr))

lseg=400;
theta=1e-6;

xkurt=kurt(xr,lseg,theta);
ykurt=kurt(yr,lseg,theta);


figure
subplot(221);plot(tsig(vwind),xr(vwind))
xlabel('time(s)');title('LP residual of clean speech, x_r(t)')
subplot(222);plot(tsig(vwind),yr(vwind))
xlabel('time(s)');title('LP residual of reverberant speech, y_r(t)')
subplot(223);plot(tsig(vwind),xkurt(vwind))
xlabel('time(s)');title('Kurosis of LP residual of clean speech, kurt(x_r)')
subplot(224);plot(tsig(vwind),ykurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of reverberant speech, kurt(y_r)')

figure
h1=subplot(211);plot(tsig(vwind),xkurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of clean speech, kurt(x)')
h2=subplot(212);plot(tsig(vwind),ykurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of reverberant speech, kurt(y_r)')
linkaxes([h1,h2]);

figure
subplot(211);plot(tsig(vwind),xr(vwind))
xlabel('time(s)');title('x_r(t)')
subplot(212);plot(tsig(vwind),xkurt(vwind))
xlabel('time(s)');title('kurt(x_r)')
figure
subplot(211);plot(tsig(vwind),yr(vwind))
xlabel('time(s)');title('y_r(t)')
subplot(212);plot(tsig(vwind),ykurt(vwind))
xlabel('time(s)');title('kurt(y_r)')


Lf=250;
niter=150;
mu=3e-6;
p=Lf;                       %filter order
ss=Lf;                      %stepsize
bs=ss+p-1;                  %blocksize
% gh=[1;zeros(p-1,1)];
gh=(1:p).';
gh=gh./sqrt(sum(abs(gh).^2));
Gh=fft([gh; zeros(ss-1,1)]);
zkurt=zeros(niter,1);

zr2=zeros(bs,1);
zr3=zeros(bs,1);
zr4=zeros(bs,1);

tic
for m=1:niter
    yrn=zeros(bs,1);
    for k=1:ss:Nsig
        yrn(1:p-1)=yrn(end-p+2:end);
        yrn(p:end)=yr(k:k+ss-1);
        
        Yrn=fft(yrn);
        cYrn=conj(Yrn);
        zrn=ifft(Gh.*Yrn);
        
        zrn(1:p-1)=0;
        zr2(p:end)=zrn(p:end).^2;
        zr3(p:end)=zrn(p:end).^3;
        zr4(p:end)=zrn(p:end).^4;
        
        Z2=sum(zr2(p:end));
        Z4=sum(zr4(p:end));
        
        zkurt(m)=max(zkurt(m),Z4/(Z2^2+1e-15)*ss);
        
%         z3y=fft(zr3).*cYrn;
%         zy=fft(zrn).*cYrn;
% 
%         gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
%         Gh=Gh+mu*gJ;
        
                
        z3y=ifft(fft(zr3).*cYrn);
        z3y(p+1:end)=0;
        zy=ifft(fft(zrn).*cYrn);
        zy(p+1:end)=0;

        gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
        Gh=Gh+mu*fft(gJ);
        
        
        Gh=Gh./sqrt(sum(abs(Gh).^2)/bs);
    end
end
toc

gh=ifft(Gh);
gh=gh(1:p);
save ifilt.mat gh
figure
plot(zkurt);
xlabel('Iterations'); ylabel('Maximum kurtosis');

load ifilt.mat gh
figure;plot((1:length(gh))/fs,gh);
xlabel('time(s)');title('Inverse filter, \bf{g}');

figure
plot((1:Nimp)/fs,h); hold on
plot((1:Nimp)/fs,filter(gh,1,h), 'LineWidth',1)
xlabel('time(s)');title('Original vs Inverse filtered RIR');
legend('Original RIR', 'Inverse Filtered RIR');

figure
plot((1:Nimp)/fs,filter(gh,1,h), 'LineWidth', 1)
xlabel('time(s)');title('Inverse filtered RIR');

z=filter(gh, 1, y);
% z=z/max(abs(z));
zr=lpres(z);
zkurt=kurt(z,lseg,theta);

figure
subplot(311);plot(tsig(vwind),xr(vwind))
xlabel('time(s)');title('LP residual of clean speech, x_r(t)')
subplot(312);plot(tsig(vwind),yr(vwind))
xlabel('time(s)');title('LP residual of reverberant speech, y_r(t)')
subplot(313);plot(tsig(vwind),zr(vwind))
xlabel('time(s)');title('LP residual of inverse filtered speech, z_r(t)')

figure
h1=subplot(311);plot(tsig(vwind),xkurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of clean speech, kurt(x)')
h2=subplot(312);plot(tsig(vwind),ykurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of reverbarant speech, kurt(y)')
h3=subplot(313);plot(tsig(vwind),zkurt(vwind))
xlabel('time(s)');title('Kurtosis of LP residual of inverse filtered speech, kurt(z)')
linkaxes([h1,h2,h3])

figure; 
subplot(311); plot(tsig(1:fs*5),x(1:fs*5));
xlabel('time(s)');title('Clean speech, x(t)')
subplot(312); plot(tsig(1:fs*5),y(1:fs*5));
xlabel('time(s)');title('Reverbarant speech, y(t)')
subplot(313); plot(tsig(1:fs*5),z(1:fs*5));
xlabel('time(s)');title('Inverse filterd speech, z(t)')

figure; 
subplot(311); spectrogram(x(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Clean speech, x(t)')
subplot(312); spectrogram(y(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Reverbarant speech, y(t)')
subplot(313); spectrogram(z(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Inverse filtered speech, z(t)')

wlen = 64e-3*fs;
hop = 8e-3*fs;
nfft = 1024;
[Sz, f_stft, t_stft] = stft(z, wlen, hop, nfft, fs);
Pz=abs(Sz).^2;
phz=angle(Sz);

% figure
% [T,F]=meshgrid(t_stft,f_stft);
% subplot(211)
% hfig=pcolor(T,F,log10(Pz));
% set(hfig, 'EdgeColor', 'None')
% subplot(212);spectrogram(x(1:fs*5),hamming(400),200,'yaxis');

ro_w=7;
a_w=5;
i_w=-ro_w:15;
wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
wS(i_w<-a_w)=0;
figure;plot(i_w,wS)

gS=.25;
epS=1e-3;
Pl=gS*filter(wS,1,Pz.').';

P_att=(Pz-Pl)./Pz;
P_att(P_att<epS)=epS;

Pxh=Pz.*P_att;
Sxh=sqrt(Pxh).*exp(1i*phz);

nu1=.01;
nu2=7;
Ez=sum(abs(Sz).^2, 1)/size(Sz,1);
Exh=sum(abs(Sxh).^2, 1)/size(Sxh,1);
P_att=ones(size(Sz,2),1);
P_att(Ez<nu1 & Ez./Exh>nu2)=1e-3;
Sxh=Sxh*spdiags(P_att,0,length(P_att),length(P_att));

[xh, t_istft] = istft(Sxh, wlen, hop, nfft, fs);
xh=xh/max(abs(xh));

figure; 
tsig=(1:Nsig)/fs;
subplot(221); plot(tsig(1:fs*5),x(1:fs*5));
xlabel('time(s)');title('Clean speech, x(t)')
subplot(223); plot(tsig(1:fs*5),y(1:fs*5));
xlabel('time(s)');title('Reverbarant speech, y(t)')
subplot(222); plot(tsig(1:fs*5),z(1:fs*5));
xlabel('time(s)');title('Inverse filtered speech, z(t)')
subplot(224); plot(tsig(1:fs*5),xh(1:fs*5));
xlabel('time(s)');title('Final processed speech, xh(t)')

figure; 
subplot(221); spectrogram(x(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Clean speech, x(t)')
subplot(223); spectrogram(y(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Reverbarant speech, y(t)')
subplot(222); spectrogram(z(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Inverse filtered speech, z(t)')
subplot(224); spectrogram(xh(1:fs*5),hamming(400),200,'yaxis');colorbar off
title('Final processed speech, xh(t)')


audiowrite('speech_x.wav',x,fs)
audiowrite('speech_y.wav',y,fs)
audiowrite('speech_z.wav',z,fs)
audiowrite('speech_xh.wav',xh,fs)