
% for Ace2N-2AA-mNeon and other GEVIs with a overshoot depolarizing
% (failing) edge

%% Intialize parameters
clear all; clc;
% Movie loading path
basepath = 'X:\91 Data and analysis\YJunqi\Screening method\Stable HEK293T cell line_ArcLight\20210708 Resting potential and sensitivity test\Dish1\cell1';
subfolder = '\154622_HASAP1_step';

pathname = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 176;         % x
nrow = 96;          % y
camera_bias = 395.8;          % background due to camera bias (100 for bin 1x1)
dt_mov = 0.9452;    % exposure time in millisecond (1058Hz)
DAQname = '\movie_DAQ.txt';
dnsamp = 20;        % downsampling rate = DAQ rate/camera rate
DAQStart = 16001;    % DAQ steps start at this datapoint (included)
frmStart = 802;     % Movie steps start at this frame (included)
frmPer = 800;       % frames period
frmCyc = 14;        % number of cycles
frmWindow_up = 140;    % number of frames to calculate time constants
frmWindow_dn = 140;
frmWindow = 140;

% load bin movie
fname = [pathname movname];
[mov, nframes] = readBinMov(fname, nrow, ncol);
mov = double(mov);
mov = mov-camera_bias;      % remove camera bias

% Set up a time axis for movie
t_mov = [0:size(mov,3)-1]*dt_mov/10^3;      % camera time axis in second

% Plot the raw whole-field intensity of the movie.
intens = squeeze(mean(mean(mov)));
figure(1)
plot(intens)
xlabel('Frame number')
ylabel('Whole-field Intensity')

%% Remove baseline drifts.  rem_pbleach calculates the minimum values in a
% sliding window, and then performs linear interpolation to make a smooth
% photobleaching trace.  The raw intensity is divided element-wise by the photobleaching trace.
% This function removes any baseline drift, not just photobleaching.
% The smoothing window should exceed the interval between stimuli.
%[intensN, pbleach] = rem_pbleach(intens, frmPer*400+1);
[intensN, pbleach] = rem_pbleach(intens, frmPer*15+1);
% figure(2)
% plot(t_mov, intensN)
% xlabel('Time')
% ylabel('Normalized intensity')
pbleach = pbleach/max(pbleach); % Normalize the pbleach trace to the maximum value

% Correct the movie for photobleaching
mov_remBleach = mov./repmat(reshape(pbleach, 1, 1, nframes), [nrow, ncol, 1]);

% Look at the average of the movie
img = mean(mov_remBleach, 3);
% figure(3)
% imshow(img, [], 'InitialMagnification', 'fit')

%% select ROI for analysis
[~, intens] = clicky(mov_remBleach, img, 'Select only 1 ROI, right click when done');
bkg = mean(intens(:,2));
background_reference=sort(reshape(mean(mov,3),ncol*nrow,1));
background_reference=mean(background_reference(round(0.01*length(background_reference)):round(0.05*length(background_reference))))
%intens_bkg=background_reference;
intens_bkg=bkg;
intens_fluo = intens(:,1)-intens_bkg;
% save clicky figure
saveas(gca,[pathname '\clicky analysis.fig']);
saveas(gca,[pathname '\clicky analysis.png']);

% load DAQ data
tmp = importdata([pathname DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
t_daq = [0:length(Vm)-1]*dt_daq/10^3;       % DAQ time axis in second

%% Kinetic analysis
t_window_dn = (0:frmWindow-1)*dt_mov;                    % millisecond
t_window_up = (0:frmWindow_up-1)*dt_mov;                    % millisecond
t_window=[0:0.9452:0.9452*0.5*(frmPer-1)];    
               % millisecond
% rising exponential
upStart = frmStart + frmPer*(0:frmCyc-1);       % marks the beginningfdn
indmat = repmat(upStart',1,frmWindow_up) + repmat(0:frmWindow_up-1,frmCyc,1);
upIntensMat = intens_fluo(indmat);                   % all traces
upIntens = mean(upIntensMat);                   % average over periods
fup_min = min(upIntens);           % minimal value
fup_max = max(upIntens);                        % maximum value
fup = (upIntens-fup_min)/(fup_max-fup_min);     % normalize response to 0-100%

% faling exponential
dnStart = frmStart + frmPer*(0.5:frmCyc-0.5);   % marks the beginning
indmat = repmat(dnStart',1,frmWindow) + repmat(0:frmWindow-1,frmCyc,1);
dnIntensMat = intens_fluo(indmat);                   % all traces
dnIntens = mean(dnIntensMat);                   % average over periods
fdn_max = mean(dnIntens(frmWindow-10:frmWindow));           % maximum value
fdn_min = min(dnIntens);                        % minimal value
fdn = (dnIntens-fdn_min)/(fdn_max-fdn_min);     % normalize response to 0-100%

indmat_all = repmat(upStart',1,frmPer) + repmat(0:frmPer-1,frmCyc,1);
indmat_upall=indmat_all(:,1:0.5*frmPer);
indmat_dnall=indmat_all(:,0.5*frmPer+1:800);
dnIntensMat_dnall = intens_fluo(indmat_dnall);   
dnIntensMat_upall = intens_fluo(indmat_upall); % all traces
dnIntens_dnall = mean(dnIntensMat_dnall); 
dnIntens_upall = mean(dnIntensMat_upall);% average over periods
fdnall_min = min(dnIntens_dnall);           % minimal value
fdnall_max = max(dnIntens_dnall);           % maximum value
fupall_min = min(dnIntens_upall);           % minimal value
fupall_max = max(dnIntens_upall);                        % maximum value
fdnall = (dnIntens_dnall-fdnall_min)/(fdnall_max-fdnall_min);     % normalize response to 0-100%
fupall = (dnIntens_upall-fupall_min)/(fupall_max-fupall_min);     % normalize response to 0-100%
% the built-in library function for an exponential assumens the form
% y = a*exp(b*x).  Values should be positive, so a > 0.  Photobleaching
% will always be a decaying exponential, so b < 0.



% rising expential
ftemp = 1-fup;          % invert to falling expential
Max = [1, 0, 0.5, 0];     % upper bound for [a,b,c,d]
Min = [0.5, -Inf, 0, -Inf];    % lower bound for [a,b,c,d]
Start = [max(ftemp), -1, 0, -.1];   % starting guess for [a,b,c,d]
% fit function
F = fit(t_window_up',ftemp','exp2','StartPoint',Start,'Lower',Min,'Upper',Max);
fitTraceUp = 1 - (F.a*exp(F.b*t_window_up)+F.c*exp(F.d*t_window_up));     % the fit function
% calculate coefficients
tauUp1 = -1/F.b;
tauUp2 = -1/F.d;
A1_up = F.a/(F.a+F.c);
A2_up = F.c/(F.a+F.c);

% falling expential
fdn = fdn/fdn(1); % normalized the average trace
Max = [1, 0, 0.5, 0];     % upper bound for [a,b,c,d]
Min = [0.5, -Inf, 0, -Inf];    % lower bound for [a,b,c,d]
Start = [max(fdn), -1, 0, -.1];   % starting guess for [a,b,c,d]
% fit function
F = fit(t_window_dn',fdn','exp2','StartPoint',Start,'Lower',Min,'Upper',Max);
fitTraceDn = F.a*exp(F.b*t_window_dn)+F.c*exp(F.d*t_window_dn);     % the fit function
% calculate tau and coefficients
tauDn1 = -1/F.b;
tauDn2 = -1/F.d;
A1_dn = F.a/(F.a+F.c);
A2_dn = F.c/(F.a+F.c);
%% Sensitivity analysis (?F/F per 100mV)
% Set up time axis under a step voltage change (-70mV to +30mV)
t_step = (0:frmPer-1)*dt_mov;                       % millisecond
StepdnStart = upStart -1;                           % marks the step beginning
inStepmat = repmat(StepdnStart',1,frmPer) + repmat(0:frmPer-1,frmCyc,1);
StepIntensMat = intens_fluo(inStepmat);             % all step traces 
StepIntens = mean(StepIntensMat);                   % average over periods
SteplowIntens = mean(StepIntens(770:790));          % average the steady state low fluorescence
StephiIntens = mean(StepIntens(370:390));           % average the steady state high fluorescence  
SquareSens = (StephiIntens-SteplowIntens)/SteplowIntens; % Square sensitivity per 100mV

%% Plot kinetics results
fig = figure(4);
set(fig,'units','normalized','outerposition',[0 0 1 1]);    % max window
% plot voltage trace
subplot(3,4,1:4);
plot(t_daq, Vm);
xlabel('time (s)'); ylabel('V_m (mV)');
axis tight;

% plot camera trace
subplot(3,4,5:8);
plot(t_mov, intens_fluo);
xlabel('time (s)'); ylabel('Intensity');
axis tight;

% plot falling exponential fit results
% plot traces
subplot(3,4,9);
plot(t_window_up, dnIntensMat');
axis tight;
% plot average and fit
subplot(3,4,10);
plot(t_window,fdnall,'.','Color','k'); hold on;    % average trace
plot(t_window_up,fitTraceDn,'Color','b');          % fit curve
title({'normalized F/F_0';
    'F(t) = A*exp(-t/\tau_1) + B*exp(-t/\tau_2)';
    ['\tau_1 = ',num2str(tauDn1,4),' ms,  \tau_2 = ',num2str(tauDn2,4),' ms'];
    ['A1 = ', num2str(A1_dn,3), ',  A2 = ',num2str(A2_dn,3)]});
xlabel('time (ms)'); ylabel('Normalized');
legend('data','fit')
axis tight;

% plot rising exponential fit results
% plot traces
subplot(3,4,11);
plot(t_window_dn, upIntensMat');
axis tight;
% plot average and fit
subplot(3,4,12);
plot(t_window,fupall,'.','Color','k'); hold on;        % average trace
plot(t_window_up,fitTraceUp,'Color','b');              % fit curve
title({'normalized F/F_0';
    'F(t) = 1 - A*exp(-t/\tau_1) + B*exp(-t/\tau_2)';
    ['\tau_1 = ',num2str(tauUp1,4),' ms,  \tau_2 = ',num2str(tauUp2,4),' ms'];
    ['A1 = ', num2str(A1_up,3), ',  A2 = ',num2str(A2_up,3)]});
xlabel('time (ms)'); ylabel('Normalized');
legend('data','fit','Location','southeast');
axis tight;

% save F-V figure
saveas(gca,[pathname '\kinetic analysis.fig']);
saveas(gca,[pathname '\kinetic analysis.png']);

%% plot sensitivity results
fig = figure(5);
set(fig,'units','normalized','outerposition',[0 0 1 1]);    % max window

% plot voltage trace
subplot(3,4,1:4);
plot(t_daq, Vm);
xlabel('time (s)'); ylabel('V_m (mV)');
axis tight;

% plot camera trace
subplot(3,4,5:8);
plot(t_mov, intens_fluo);
xlabel('time (s)'); ylabel('Intensity');
axis tight;

% plot all step traces
figure(5);
subplot(3,4,9:10);
plot(t_step, StepIntensMat');
xlabel('time (ms)'); ylabel('Intensity');
axis tight;

% plot average and calculate the sensitivity per 100mV
subplot(3,4,11:12);
plot(t_step, StepIntens');
title({['? F/F  = ',num2str(SquareSens*100,4),' % per 100mV'];});
xlabel('time (ms)'); ylabel('Intensity');
axis tight;

% save step sensitivity figure
saveas(gca,[pathname '\step sensitivity analysis.fig']);
saveas(gca,[pathname '\step sensitivity analysis.png']);
%% Write a xlsx file
xlswrite([pathname '\0analysis.xlsx'],{'Sensitivity', 'tau_1', 'A1','tau_2','A2','tau_3','B1','tau_4','B2','Brightness'},'Average','A1:J1');
xlswrite([pathname '\0analysis.xlsx'],SquareSens*-100,'Average','A2');
xlswrite([pathname '\0analysis.xlsx'],tauDn1,'Average','B2');
xlswrite([pathname '\0analysis.xlsx'],A1_dn,'Average','C2');
xlswrite([pathname '\0analysis.xlsx'],tauDn2,'Average','D2');
xlswrite([pathname '\0analysis.xlsx'],A2_dn,'Average','E2');
xlswrite([pathname '\0analysis.xlsx'],tauUp1,'Average','F2');
xlswrite([pathname '\0analysis.xlsx'],A1_up,'Average','G2');
xlswrite([pathname '\0analysis.xlsx'],tauUp2,'Average','H2');
xlswrite([pathname '\0analysis.xlsx'],A2_up,'Average','I2');
xlswrite([pathname '\0analysis.xlsx'],StephiIntens','Average','J2');
