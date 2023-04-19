

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY GIANMARCO PINTON 
% FIRST CREATED: 2017-05-01
% LAST MODIFIED: 2020-04-05
% NONLINEAR PROPAGATION IN HOMOGENEOUS MEDIUM
% MODIFIED ANGULAR SPECTRUM 
% RUSANOV
% FREQUENCY DOMAIN ATTENUATION AND DISPERSION
% KRAMERS KRONIG ATTENUATION AND DISPERSION
% ABSORBING BOUNDARY LAYERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load hydrophone measurements of imaging transducer taken near transducer face
load source_pressure_4C1_30V_M1_9_100MHz_single
p0=maxmax(pressure);

dTT=1/100e6;    % this is the pulse duration (PD)
dXX=(mean(diff(elPos)))*1e-3;

a=32e-3; % source radius
f0=3e6; % center frequency of pulse
omega0 = f0*2*pi;
c0=1500; % speed of sound
lambda=c0/f0; % wavelength
dX=lambda/3; dY=dX; % grid spacing

dZ=dX*6;
dT=dX/5/c0; % 2*pi/omega0/20;
k=omega0/c0;
beta=3.5; % nonlinear coefficient
rho0=1000; % equilibrium density
N=beta/2/c0^3/rho0;
%p0=1e6; % transmitter pressure
prop_dist=4.5e-2; % propagation distance in cm

%% INTERPOLATE PRESSURE %%
tmp=interp1easy(pressure(:,round(end/2),round(end/2)),dTT/dT);
p2=zeros(length(tmp),size(pressure,2),size(pressure,3));
for j=1:size(pressure,2)
    for k=1:size(pressure,3)
        p2(:,j,k)=interp1easy(pressure(:,j,k),dTT/dT);
    end
end
p2=single(p2);

outimg = squeeze(p2(round(end/2),:,:));
outimg = interp2easy(outimg,dXX/dX,dXX/dY);
p3 = zeros(size(p2,1),size(outimg,1),size(outimg,2),'single');

for i=1:size(p2,1)
  p3(i,:,:) = interp2easy(squeeze(p2(i,:,:)),dXX/dX,dXX/dY);
end
clear p2;
%p3=p3(100:end,:,:);
%p3=flipdim(p3,1);
p4=permute(p3,[2 3 1]);

% GRID SIZE
nX=size(p4,1); nY=size(p4,2); nT=size(p4,3)% grid size

if(mod(nX,2)-1) % keep the grid odd
    nX=nX+1;
    p4(nX-1:nX,:,:)=0;
end
if(mod(nY,2)-1)
    nY=nY+1;
    p4(:,nY-1:nY,:)=0;
end
if(mod(nT,2)-1)
    nT=nT+1;
    p4(:,:,nT-1:nT,:,:)=0;
end

xaxis=(0:nX-1)*dX;xaxis=xaxis-mean(xaxis);
yaxis=(0:nY-1)*dY;yaxis=yaxis-mean(yaxis);
taxis = (0:nT-1)*dT;

apa=p4;
% figure(1), imagesc(t,yaxis,(squeeze(apa(:,round(nY/2),:)))), 
% title('initial condition (Pa)'), xlabel('t (s)'), ylabel('y (m)'), colorbar

%% PRECALCULATE MODIFIED ANGULAR SPECTRUM 
[HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
%% ABSORBING BOUNDARY LAYER
[abl] = precalculate_abl(nX,nY,nT);
%% PRECALCULATE ABSORPTION/DISPERSION FILTER
%  alpha0=2.17e-3; % dB/MHz^2/cm for water, note attenuation law proportional to f^2 here
%  [afilt3d] = precalculate_ad(alpha0,nX,nY,nT,dZ,dT);
alpha0=0.0; pow=1; % dB/MHz^pow/cm power law attenuation here
[afilt3d] = precalculate_ad_pow2(alpha0,nX,nY,nT,dZ,dT,c0,omega0/2/pi,pow); 

% 
% figure(1), plot(real(squeeze(afilt3d(round(end/2),round(end/2),:))))
% figure(2), plot(imag(squeeze(afilt3d(round(end/2),round(end/2),:))))
% figure(3), imagesc(real(abl(:,:,round(end/2))))
%

%% MARCH ALONG PROPAGATION AXIS %%
apaz=apa;   % apaz is pressure in Z
zvec=dZ;    % zvec is currently just the step in Z
dZa=dZ;
cc=1;
isppa = zeros(round(prop_dist/dZ) + 5,1);
fs = 1/dT;
nfft = nT;
fft_freq = fs/2*[0:2/nfft:1];
while(sum(zvec)<prop_dist)
    disp(sum(zvec))
    disp(prop_dist)
    % N*dZ/dT*max(max(max(apaz)))
    if(N*dZ/dT*max(max(max(apaz)))>0.1) % stability criterion
        disp('Stability criterion violated, retrying with smaller step size')
        dZ=0.075*dT/max(max(max(apaz)))/N % adjusting dZ by .75
        cc=1;
        apaz=apa;
        zvec=dZ;
        
        isppa = zeros(round(prop_dist/dZ) + 5,1);
        fs = 1/dT;
        nfft = nT;
        fft_freq = fs/2*[0:2/nfft:1-2/nfft];
        
        [HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
        [abl] = precalculate_abl(nX,nY,nT);
        [afilt3d] = precalculate_ad_pow2(alpha0,nX,nY,nT,dZ,dT,c0,omega0/2/pi,pow); 
    end
    zvec(cc)=dZ; sum(zvec);
    disp(['Propagation distance = ' num2str(sum(zvec)) ' m'])

    apaz=march_asr2(apaz,dZ,dT,N,HH,abl,afilt3d);   % this is some application (like fullwave) that does the thing

    % store all pulse waveforms to process later
    % insert code here to calc time averaged intensity integral, and
    % compare it to the last one
    % remember that we are marching through space, not time

     p_pulse = squeeze(apaz(round(nX/2),round(nY/2),:));
     pfft = abs(fftshift(fft(p_pulse)));
     pfft = pfft(ceil(nfft/2):end);
     pfft = pfft/(max(pfft));
%     p_derated = derate(p_pulse,alpha0,f0,sum(zvec));
%     isppa_z = compute_isppa(p_derated,taxis, c0, rho);
%     isppa(cc) = isppa_z;

    figure(1), imagesc(t,yaxis,(squeeze(apaz(round(nX/2),:,:)))), xlabel('t (s)'), ylabel('m in y'), title([num2str(sum(zvec)) 'm']), colorbar
    %str = ['print -djpeg movie/yimg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    %figure(2), imagesc(t,xaxis,(squeeze(apaz(:,round(nY/2),:)))), xlabel('t (s)'), ylabel('m in x'), title([num2str(sum(zvec)) 'm']), colorbar
    %str = ['print -djpeg movie/ximg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    figure(3), plot(taxis,p_pulse), xlabel('t (s)'), ylabel('Pa'), title( [num2str(sum(zvec)) 'm'])
    %str = ['print -djpeg movie/timg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    figure(4), plot(fft_freq, pfft),xlabel('Frequency (Hz)'), ylabel('Normalized FFT'), title('FFT of Waveform')
    %drawnow
    cc=cc+1;
end
save(['pressure_focal_dist_alpha_' num2str(alpha0) '_f0_' num2str(f0 * 1e-6) 'MHz_'  '.mat'], 'alpha0', 'f0', 'apaz', 'dZ', 'prop_dist', 'taxis', 'isppa');