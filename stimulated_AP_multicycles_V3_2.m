% pay attention to the ncol,nrow,dt_mov,dnsamp, and cycle!
% this script is for df/F probe.
% Luxin Peng, last edited in 20181106
% Luxin Peng, last edited in 20190415
% Luxin Peng, last edited in 20190415
% Luxin Peng, last edited in 20200602
%             change the methods for fingding peaks and normalization
%             
% Zoulab

%% load and select
clear all; clc;
% Movie loading path
basepath = 'H:\ZouLab\YJunqi\Sensitivity\Cepheid\20220728 Response of sensors to action potentials\Dish1\neuron1\';
subfolder = '184250_300pA-10ms-80 times_(Baseline 0 pA)\'
pathname = [basepath subfolder '\']
a = importdata([pathname '\movie_DAQ.txt']);
if exist([basepath '\patch param.txt'])
    b = importdata([basepath '\patch param.txt']);
    datab = b.data;
    Cm = datab(1,5);%pF
    Rm = datab(1,4);
    Ra = datab(1,3);
else
    Rm = 'N.A.'
    Cm = 'N.A.'
    Ra = 'N.A.'
end

dire = -1;% dire = 1 means positive GEVI; dire = -1 means negative GEVI

load ([pathname 'Current injection_300pA-10ms-80 times_(Baseline 0 pA)\matlab variables.mat']);
% constants
% cycles = 80;
mode = 484;
genre = 2;%1 = optopatch;2 = patch
if genre == 1
if mode == 484
    movname = '\movie.bin';
    ncol = 312;         % x
    nrow = 208;         % y
    camera_bias = 400;  % background due to camera bias (100 for bin 1x1)
    dt_mov = 2.0658;    % exposure time in millisecond (484 Hz)
    Fs = samprate;
else 
    movname = '\movie.bin';
    ncol = 176;         % x
    nrow = 96;         % y
    camera_bias = 400;  % background due to camera bias (100 for bin 1x1)
    dt_mov = 0.9549;    % exposure time in millisecond (??Hz)
    Fs = samprate;
end
else
    if mode == 484
    movname = '\movie.bin';
    ncol = 312;         % x
    nrow = 208;         % y
    camera_bias = 400;  % background due to camera bias (100 for bin 1x1)
    dt_mov = 2.0658;    % exposure time in millisecond (484 Hz)
    Fs = samprate;
else 
    movname = '\movie.bin';
    ncol = 648;         % x
    nrow = 648;         % y
    camera_bias = 400;  % background due to camera bias (100 for bin 1x1)
    dt_mov = 6.4020;    % exposure time in millisecond (??Hz)
    Fs = samprate;
end
end
DAQname = '\movie_DAQ.txt';
dnsamp = Fs/(1000/dt_mov);        % downsampling rate = DAQ rate/camera rate
dnsamp = round(dnsamp);

% load DAQ data
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
t_daq = [0:length(Vm)-1]*dt_daq/10^3;       % DAQ time axis in second
a=importdata([pathname '\movie_DAQ.txt']);
data=a.data;
AI_scaled=data(:,1);
AI_10Vm=data(:,2)*100;
time=(1:length(AI_scaled)')'./Fs;
figure;
set(gcf,'outerposition',get(0,'screensize'));
plot(time,AI_scaled,time,AI_10Vm);
legend('AI\_scaled','Vm (mV)','Location','Northeast');
hold on
xlim=[0,max(time)];
ylim=[-70,-60];
xL=xlim;yL=ylim;
set(gca,'xtick',[0:5:max(time)])
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
axis([xL yL])
axis tight
saveas(gca,[pathname '0 waveform of AI.fig']);
saveas(gca,[pathname '0 waveform of AI.png']);
%% loading the video movie

% load movie
fname = [pathname movname];
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = double(mov);img = mean(mov, 3);

% select ROI for analysi
% s
[~, intens_raw] = clicky(mov, img, 'elect only 1 ROI, right click when done');
background = mean(intens_raw(:,size(intens_raw,2)));
title('Please left click twice for the boundary of X axis of the signal')
[boundary,~] = ginput(2);
boundary(boundary<=1) = 1 ;
boundary(boundary>=size(intens_raw,1)) = size(intens_raw,1) ;
boundary = round(boundary);
% save clicky figure
saveas(gca,[pathname '\1 clicky analysis_stimulated_AP.fig']);
saveas(gca,[pathname '\1 clicky analysis_stimulated_AP.png']);
len = size(intens_raw,1);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second
intens = intens_raw(boundary(1):boundary(2),:);
t_mov_selected = t_mov(boundary(1):boundary(2));
figure();
plot(t_mov_selected,intens(:,1)-intens(:,size(intens_raw,2)));
box off
axis tight
title(['Selected region, from frame ' num2str(boundary(1)) ' to ' num2str(boundary(2)) ' from the raw movie'])
xlabel('Time (s)')
ylabel('Intensity, without background');
saveas(gca,[pathname '\2 selected region, from frame ' num2str(boundary(1)) ' to ' num2str(boundary(2)) '.fig']);
saveas(gca,[pathname '\2 selected region, from frame ' num2str(boundary(1)) ' to ' num2str(boundary(2)) '.png']);
period_param = 0.9; % this is for peak searching, now we have a selection step from the height, so we could increase the peak-finding frequency
%% Check the FWHM of AP with elelctrophysiological signal
nback = 40*round(dnsamp);  % number of frames to extend kernel behind AP peak
nfront = 40*round(dnsamp); % number of frames to extend kernel in front of AP peak
% The spikefind_corr algorithm finds spikes in two iterations.  First it
% uses a simple threshold-and-max procedure (spikefind.m).  It generates a
% kernel which it convolves with the original data.  The maxima in the
% convolved data are taken as the revised spike times.  This procedure
% identifies spikes by their global similarity to a template, rather than
% by the location of a single peak, which could be subject to noise.
if boundary(2)*dnsamp > length(AI_10Vm)
    AI_10Vm_interest = AI_10Vm((boundary(1)-1)*dnsamp+1:end);
    t_daq_interest = t_daq((boundary(1)-1)*dnsamp+1:end);
else
    AI_10Vm_interest = AI_10Vm((boundary(1)-1)*dnsamp+1:boundary(2)*dnsamp);
    t_daq_interest = t_daq((boundary(1)-1)*dnsamp+1:boundary(2)*dnsamp);
end

[~, spikeT1_elec] = findpeaks(AI_10Vm_interest, 'MinPeakDistance',(period_param*length(AI_10Vm))./cycles);
figure();
findpeaks(AI_10Vm_interest, 'MinPeakDistance',(period_param*length(AI_10Vm))./cycles);hold on
title('Left click one time for the threshold');
[~,threshold_for_peak_elec] = ginput(1);
plot([1 length(AI_10Vm_interest)],[threshold_for_peak_elec threshold_for_peak_elec],'r','LineWidth',2);hold on
legend('Threshold','Peak from period')
[~, spikeT2_elec] = findpeaks(AI_10Vm_interest, 'MinPeakHeight',threshold_for_peak_elec);% last 500 points are excluded now for better analysis

spikeT_real_elec = intersect(spikeT1_elec,spikeT2_elec); % find the intersection of the peaks location
%kernel = kernel1;spikeT = spikeT1;
nspike_elec = length(spikeT_real_elec);
Lk_elec = nfront+nback+1;
kernel_stack_elec = zeros(Lk_elec,size(spikeT_real_elec,1));
for n = 1:size(spikeT_real_elec,1);
    kernel_stack_elec(:,n) = AI_10Vm_interest(spikeT_real_elec(n)-nfront:spikeT_real_elec(n)+nback);
    base_elec(n) = mean(kernel_stack_elec(1:30*round(dnsamp),n));
    peak_elec(n) = max(kernel_stack_elec(:,n));
    HM_elec(n) = 0.5*(base_elec(n)+peak_elec(n));
    temp = kernel_stack_elec(:,n);
    locs_in = find(temp>= HM_elec(n));
    FWHM_elec(n) = (locs_in(end)-locs_in(1)+1)/Fs*1000;
    Amplitude_elec(:,n) = abs(base_elec(:,n)-peak_elec(:,n));
%   plot([locs_in(1) locs_in(end)+1],[temp(locs_in(1)) temp(locs_in(end)+1)])
    clear locs_in
    clear temp
end
figure()
plot(FWHM_elec)
box off
axis tight
xlabel('Spike num')
ylabel('FWHM (ms)')
title(['FWHM calculation from patch clamp data, acquisition rate = ' ,num2str(Fs) ' Hz'])
saveas(gca,[pathname '\3 FWHM calculation_from patch.fig'])
saveas(gca,[pathname '\3 FWHM calculation_from patch.png'])
%% Drawing mean AP with every AP trace
figure()
plot([0:1:size(kernel_stack_elec,1)-1]'*dt_daq,kernel_stack_elec,'Color',[0.8,0.8,0.8],'LineWidth',0.8);hold on;
plot([0:1:size(kernel_stack_elec,1)-1]'*dt_daq,mean(kernel_stack_elec,2),'Color',[0,0,0],'LineWidth',3);
box off
axis([0 dt_daq*size(kernel_stack_elec,1)-1 -70 60])
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')
title(['Average AP elelctrophysiological signal (Vm) from ' num2str(size(kernel_stack_elec,2)) ' APs'])
saveas(gca,[pathname '\3_1 average AP patch signal.fig'])
saveas(gca,[pathname '\3_1 average AP patch signal.png'])
%% AP chracterization (from electro. data)
Rp = mean(AI_10Vm(1:400));
spikelocs = [];
Curr = unique(patchAO*2000);
%Curr(2:end) = Curr(2:end)+60;
if Cm ~=0;
    Cdensity = unique(patchAO*2000./Cm);
else
    Cdensity = NaN;
end
spikelocs = [];
clear Vm
    for n = 1:size(kernel_stack_elec,2)
        Vm(:,n) = kernel_stack_elec(:,n);
        Vm_smooth(:,n) = smooth(smooth(smooth(Vm(:,n),5),5),5);
        diff3(:,n) = diff(Vm_smooth(:,n),3);
        diff1(:,n) = diff(Vm_smooth(:,n),1);
        [pks,locs]=findpeaks(diff3(766:810,n),'minpeakheight',0.02);
        if length(pks)==0
            threshold(:,n) = NaN;
        else
            threshold(:,n) = Vm(min(locs)+765+3,n);
            threshold_loc(:,n) = (min(locs)+3)+765;
        end
        rising(:,n)= Vm(min(find(diff1(:,n)==max(diff1(:,n)))),n);
        diff1_raw(:,n) =  diff(Vm(:,n),1);
        tmp2 = diff1_raw(:,n);
        rising_speed (:,n) = tmp2((min(find(diff1(:,n)==max(diff1(:,n))))))*Fs/1000;
        numofspike(:,n) = length(spikefind(Vm(:,n),0));
        if ~isnan(threshold(:,n))
            APD_threshold(:,n) = (min(find(Vm(min(locs)+3:end,n)<threshold(:,n)))+2)/Fs*1000;
            downstream_loc(:,n) = threshold_loc(:,n)*Fs+min(find(Vm(min(locs)+1:end,n)<threshold(:,n)));
        end
        if n~=1
            base(:,n) = mean(Vm(end-300:end-100,n-1));
        else
            base(:,n) = Rp;
        end
        
        if ~isnan(threshold(:,n))
            max_V_loc(:,n)= min(spikefind(Vm(:,n),0));
            max_V(:,n) = Vm(max_V_loc(:,n),n);
%             Amplitude(:,n) = abs(base_elec(:,n)-max_V(:,n));
            half_V(:,n) = (max_V(:,n)+base(:,n))*0.5;
            downstream_loc_FWHM(:,n) = min(find(Vm(max_V_loc(:,n)+1:end,n)<half_V(:,n)));
            upstream_loc_FWHM(:,n) = size(find(Vm(1:(max_V_loc(:,n)),n)>half_V(:,n)),1)-1;
            FWHM_reference(:,n) = (downstream_loc_FWHM(:,n)+upstream_loc_FWHM(:,n)+1)/Fs*1000;
        end
        spikelocs1 = spikefind(Vm(:,n),0);
        spikelocs1 = spikelocs1+length(headPts)+(hiPts+lowPts)*(n-1);
        
        spikelocs = [spikelocs spikelocs1];
        clear spikelocs1;
        clear pks;
        clear locs;
        clear tmp2;
        
    end

segment = min(find(~isnan(threshold)>0));
Rheob = Curr(segment+1);
if ~isnan(Cdensity)
    Rheob_den = Cdensity(segment+1);
else
    Rheob_den = NaN;
end

% spikesum = sum(numofspike);
% spikevalue = [(spikelocs)' (spikelocs./Fs)' AI_10Vm(spikelocs)];
% max_AP = spikevalue(1,3);
% Tp_1st = threshold((segment));
% R_speed_1st = rising_speed((segment));
% Amp_1st = Amplitude((segment));
% FWHM_1st = FWHM_reference((segment));
% APD_thres_1st = APD_threshold((segment));
% base_1st = base((segment));
% Tp_1st
% segment
xlswrite([pathname 'electro_analysis.xlsx'],{'Resting potential (mV)','Spike number','Amplitude (mV)', 'Base (mV)', 'Peak voltage (mV)','Max. rising speed (mV/ms)','Vm (Max. rising speed) ','Threshold (mV)','FWHM (ms)','Ra (MOhm)','Rm (MOhm)','Cm (pF)'},'Single AP','A1');
xlswrite([pathname 'electro_analysis.xlsx'],Rp','Single AP','A2');
xlswrite([pathname 'electro_analysis.xlsx'],[1:1:n]','Single AP','B2');
xlswrite([pathname 'electro_analysis.xlsx'],Amplitude_elec','Single AP','C2');
xlswrite([pathname 'electro_analysis.xlsx'],base_elec','Single AP','D2');
xlswrite([pathname 'electro_analysis.xlsx'],peak_elec','Single AP','E2');
xlswrite([pathname 'electro_analysis.xlsx'],rising_speed','Single AP','F2');
xlswrite([pathname 'electro_analysis.xlsx'],rising','Single AP','G2');
xlswrite([pathname 'electro_analysis.xlsx'],threshold','Single AP','H2');
xlswrite([pathname 'electro_analysis.xlsx'],FWHM_elec','Single AP','I2');
xlswrite([pathname 'electro_analysis.xlsx'],Ra,'Single AP','J2');
xlswrite([pathname 'electro_analysis.xlsx'],Rm','Single AP','K2');
xlswrite([pathname 'electro_analysis.xlsx'],Cm','Single AP','L2');

% xlswrite([pathname 'analysis.xlsx'],intens_mean','Average','I2');

xlswrite([pathname 'electro_analysis.xlsx'],{'Parameter','Amplitude (mV)', 'Base (mV)', 'Peak voltage (mV)','Max. rising speed (mV/ms)','Vm (Max. rising speed) ','Threshold (mV)','FWHM (ms)','path'}','Average','A1');
xlswrite([pathname 'electro_analysis.xlsx'],{'Average'},'Average','B1');
xlswrite([pathname 'electro_analysis.xlsx'],mean(Amplitude_elec),'Average','B2');
xlswrite([pathname 'electro_analysis.xlsx'],mean(base_elec),'Average','B3');
% xlswrite([pathname 'analysis.xlsx'],baseline_vec','single_trace','B2');
xlswrite([pathname 'electro_analysis.xlsx'],mean(peak_elec),'Average','B4');
xlswrite([pathname 'electro_analysis.xlsx'],mean(rising_speed),'Average','B5');
xlswrite([pathname 'electro_analysis.xlsx'],mean(rising),'Average','B6');
xlswrite([pathname 'electro_analysis.xlsx'],mean(threshold),'Average','B7');
xlswrite([pathname 'electro_analysis.xlsx'],mean(FWHM_elec),'Average','B8');

xlswrite([pathname 'electro_analysis.xlsx'],{'S.E.M.'},'Average','C1');
xlswrite([pathname 'electro_analysis.xlsx'],std(Amplitude_elec)./sqrt(n),'Average','C2');
xlswrite([pathname 'electro_analysis.xlsx'],std(base_elec)./sqrt(n),'Average','C3');
% xlswrite([pathname 'analysis.xlsx'],baseline_vec','single_trace','B2');
xlswrite([pathname 'electro_analysis.xlsx'],std(peak_elec)./sqrt(n),'Average','C4');
xlswrite([pathname 'electro_analysis.xlsx'],std(rising_speed)./sqrt(n),'Average','C5');
xlswrite([pathname 'electro_analysis.xlsx'],std(rising)./sqrt(n),'Average','C6');
xlswrite([pathname 'electro_analysis.xlsx'],std(threshold)./sqrt(n),'Average','C7');
xlswrite([pathname 'electro_analysis.xlsx'],std(FWHM_elec)./sqrt(n),'Average','C8');

xlswrite([pathname 'electro_analysis.xlsx'],{pathname},'Average','B9');

%% generate kernel
[intensN, pbleach] = rem_pbleach(dire*((intens(:,1)-intens(:,2))), round(4*size(intens_raw,1)/cycles)+1);% Remove the photobleaching, '-1' means find the max
pbleach = dire*pbleach;
%intensN = intens(:,1)-intens(:,2); % if we don't want to remove the
% photobleaching and do the anaalysis, run this line before line 97
intensN = intensN./mean(intensN(1:25));
figure(2)
plot(t_mov_selected, intensN)
xlabel('Time')'
ylabel('Normalized intensity')
title('Normalized trace after removing photobleach')
box off
axis tight
saveas(gca, [pathname '4 after rem_pbleach.fig']);
saveas(gca, [pathname '4 after rem_pbleach.png']);
close(gcf)

intensN_norem = (intens(:,1)-intens(:,2));
intensN_norem = intensN_norem./mean(intensN_norem(1:25));

figure(2)
plot(t_mov_selected, intensN_norem)
xlabel('Time')
ylabel('Normalized intensity')
title('Normalized trace without removing photobleach')
box off
axis tight
saveas(gca, [pathname '5 Normalize, but no rem_pbleach.fig']);
saveas(gca, [pathname '5 Normalize, but no rem_pbleach.png']);
close(gcf)

% Find positions of spikes
nback = 25;  % number of frames to extend kernel behind AP peak
nfront = 25; % number of frames to extend kernel in front of AP peak
% The spikefind_corr algorithm finds spikes in two iterations.  First it
% uses a simple threshold-and-max procedure (spikefind.m).  It generates a
% kernel which it convolves with the original data.  The maxima in the
% convolved data are taken as the revised spike times.  This procedure
% identifies spikes by their global similarity to a template, rather than
% by the location of a single peak, which could be subject to noise.
cutoff_num = 0;% this is the number of frame in the ending of the train that would be excluded
[~, spikeT1] = findpeaks(dire*(intensN(1:end-cutoff_num)),'MinPeakDistance',((period_param*length(intensN))./cycles));
figure();
findpeaks(dire*(intensN(1:end-cutoff_num)),'MinPeakDistance',((period_param*length(intensN))./cycles));hold on
title('Left click one time for the threshold');
[~,threshold] = ginput(1);
plot([1 length(intensN)],[threshold threshold],'r','LineWidth',2);hold on
legend('Threshold','Peak from period')
[~, spikeT2] = findpeaks(dire*(intensN(1:end-cutoff_num)), 'MinPeakHeight',threshold);% last 500 points are excluded now for better analysis
saveas(gca,[pathname '6 find_the_peaks.fig']);
saveas(gca,[pathname '6 find_the_peaks.png']);

spikeT_real = intersect(spikeT1,spikeT2); % find the intersection of the peaks location
%kernel = kernel1;spikeT = spikeT1;
%spikeT_real = spikeT_real(2:end);
%%
nspike = length(spikeT_real);
Lk = nfront+nback+1;

kernel_stack = zeros(Lk,size(spikeT_real,1));
for n = 1:size(spikeT_real,1);
    kernel_stack(:,n) = intensN(spikeT_real(n)-nfront:spikeT_real(n)+nback);
end
kernel_mean = mean(kernel_stack,2);
% Lk = length(kernel1); % = nback + nfront + 1

% for normalize to the base
baseframe = 15;
kernel_mean = kernel_mean-mean(kernel_mean(1:baseframe))+1;% moving the baseline to "1"
% figure();
% plot(kernel_mean);
% box off
% axis tight
% xlabel('frame')
% ylabel('delta F')
% title ('average of trails')
% saveas(gca,[pathname '8 mean_AP.fig']);
% saveas(gca,[pathname '8 mean_AP.png']);

figure();
plot(t_mov_selected,pbleach);
box off
axis tight
xlabel('Time (s)')
ylabel('Intensity')
title ('Line for removing photobleach')
saveas(gca,[pathname '8 line_of_pbleach.fig']);
saveas(gca,[pathname '8 line_of_pbleach.png']);
% Sensitivity = (max(kernel1)-mean(kernel1(1:20)))
% Brightness = (mean(intens(1:100,1)-mean(intens(:,2))))
% Average_SNR = Sensitivity./std(kernel1(1:20))
% fname = [pathname '\Sensitivity.txt'];       % file name
% fid = fopen(fname, 'w');                    % open file
% fprintf(fid,'Brightness (average) =%.5f\r\n',Brightness);
% fprintf(fid,'Sensitivity (average) =%.5f\r\n',Sensitivity);
% fprintf(fid,'SNR (average) =%.5f\r\n',Average_SNR);% output result
% fclose(fid);
%% peak-by-peak analysis
Brightness = mean(pbleach(1:100,1));% without background
% intens_matrix = zeros(Lk,nspike);
t_mov2 = [0:dt_mov:dt_mov*(size(kernel_stack,1)-1)];
sensitivity_vec = zeros(1,nspike);
baseline_vec = zeros(1,nspike);
dF_vec = zeros(1,nspike);
SNR_vec = zeros(1,nspike);
RMSE_vec = zeros(1,nspike);
for n = 1:nspike
    baseline_vec(n) = mean(kernel_stack(1:baseframe,n));
    if dire == 1
        dF_vec(n) = max(kernel_stack(:,n))- baseline_vec(n);%for positive GEVIs
    else
        dF_vec(n) = min(kernel_stack(:,n))- baseline_vec(n);%for negative GEVIs
    end
    sensitivity_vec(n) = dF_vec(n)/baseline_vec(n);
    % for the detrending
    [xData, yData] = prepareCurveData( [], kernel_stack(1:baseframe,n));
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    RMSE(n)=gof.rmse;
    SNR_vec(n) = abs(dF_vec(n))./RMSE(n);
end



SNR_single = mean(SNR_vec);
sensitivity_single = mean(sensitivity_vec);
baseline_averaged = mean(kernel_mean(1:baseframe));
if dire == 1
    dF_averaged = max(kernel_mean)- baseline_averaged;% for positive GEVIs
else
    dF_averaged = min(kernel_mean)- baseline_averaged;% for negative GEIVs
end
sensitivity_averaged = (dF_averaged)/baseline_averaged;
[xData, yData] = prepareCurveData( [], kernel_mean(1:baseframe));
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
RMSE_mean=gof.rmse;
SNR_averaged = abs(dF_averaged)/RMSE_mean;

figure();
plot(kernel_stack);
box off
axis tight
xlabel('frame')
ylabel('delta F')
title (['average of trails, \DeltaF/F = ' ,num2str(sensitivity_single*100,4), ' +- ' ,num2str(std(sensitivity_vec)./sqrt(n)*100,3) '%'])
saveas(gca,[pathname '7 trails of single APs.fig']);
saveas(gca,[pathname '7 trails of single APs.png']);

figure()
plot(t_mov2,kernel_mean);
box off
axis tight
title({['\DeltaF/F  = ',num2str(sensitivity_averaged*100,4),' % per AP (100 mV), averaged ' num2str(nspike) ' trails'];});
xlabel('Time (ms)')
ylabel('Fluoresent intenstity (w/o background)')
saveas(gca,[pathname '\9 mean AP_with_sensitivity.fig']);
saveas(gca,[pathname '\9 mean AP_with_sensitivity.png']);
%% Drawing mean fluoresence response with every fluorescence change
kernel_stack_normalized = kernel_stack./repmat(kernel_stack(1,:),size(kernel_stack,1),1);% normalization to 1st frame in every response
figure()
plot([0:1:size(kernel_stack_normalized,1)-1]'*dt_mov,kernel_stack_normalized,'Color',[0.8,0.8,0.8],'LineWidth',0.8);hold on;
plot([0:1:size(kernel_stack_normalized,1)-1]'*dt_mov,mean(kernel_stack_normalized,2),'Color',[0.99,0,0],'LineWidth',3);
box off
axis tight
xlabel('Time (ms)')
ylabel('Normalized fluorescence')
title(['Average AP elelctrophysiological signal (Vm) from ' num2str(size(kernel_stack_normalized,2)) ' APs'])
saveas(gca,[pathname '\9_1 average fluorscence response.fig'])
saveas(gca,[pathname '\9_1 average fluorscence response.png'])
%% plot Vm-t and Intensity-t
figure
subplot(3,1,1)
plot(mat2gray(img));% gray img plot
imshow(mat2gray(img));
title('selected cell');
subplot(3,1,2)
plot(t_daq_interest,AI_10Vm_interest);
title('membrane voltage signal');
xlabel('t(s)');
ylabel('Vm (mV)');
axis tight;
hold on
% xL=[0,max(t_mov)];yL=[min(Vm)-10,max(Vm)+10];% range of figure
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
subplot(3,1,3)
plot(t_mov_selected,intensN);
hold on
title('Optical signal');
xlabel('t(s)');
ylabel('Intensity (Normalized)');
axis tight;
% xL=[0,max(t_mov)];yL=[min(intens(:,1)-background)-2,max(intens(:,1)-background)+2];% range of figure
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
set(gcf,'outerposition',get(0,'screensize'));

saveas(gca,[pathname '\10 integrated analysis.fig'])
saveas(gca,[pathname '\10 integrated analysis.png'])
%% Calculation of the FWHM with linear interpolation
upsamplingrate = 20;
intensN_interp = interp1(t_mov_selected,intensN,[min(t_mov_selected):dt_mov/1000/upsamplingrate:max(t_mov_selected)]);

nback = 40*upsamplingrate;  % number of frames to extend kernel behind AP peak
nfront = 40*upsamplingrate; % number of frames to extend kernel in front of AP peak
% The spikefind_corr algorithm finds spikes in two iterations.  First it
% uses a simple threshold-and-max procedure (spikefind.m).  It generates a
% kernel which it convolves with the original data.  The maxima in the
% convolved data are taken as the revised spike times.  This procedure
% identifies spikes by their global similarity to a template, rather than
% by the location of a single peak, which could be subject to noise.
cutoff_num = 0;% this is the number of frame in the ending of the train that would be excluded
[~, spikeT1_interp] = findpeaks(dire*(intensN_interp(1:end-cutoff_num)), 'MinPeakDistance',(period_param*length(intensN_interp))./cycles);
figure();
findpeaks(dire*(intensN_interp(1:end-cutoff_num)), 'MinPeakDistance',(period_param*length(intensN_interp))./cycles);hold on
title('Left click one time for the threshold');
[~,threshold_interP] = ginput(1);
plot([1 length(intensN_interp)],[threshold_interP threshold_interP],'r','LineWidth',2);hold on
legend('Threshold','Peak from period')
[~, spikeT2_interp] = findpeaks(dire*(intensN_interp(1:end-cutoff_num)), 'MinPeakHeight',threshold_interP);% last 500 points are excluded now for better analysis

spikeT_real_interp = intersect(spikeT1_interp,spikeT2_interp); % find the intersection of the peaks location
%kernel = kernel1;spikeT = spikeT1;
nspike_interp = length(spikeT_real_interp);
Lk_interp = nfront+nback+1;
kernel_stack_interp = zeros(Lk_interp,size(spikeT_real_interp,2));
for n = 1:size(spikeT_real_interp,2);
    kernel_stack_interp(:,n) = intensN_interp(spikeT_real_interp(n)-nfront:spikeT_real_interp(n)+nback);
    base_interp(n) = mean(kernel_stack_interp(1:30*upsamplingrate,n));
    if dire == 1
        peak_interp(n) = max(kernel_stack_interp(:,n));
        HM_interp(n) = 0.5*(base_interp(n)+peak_interp(n));
        temp = kernel_stack_interp(:,n);
          temp3 = find(temp < HM_interp(n));
        locs_in = temp3(max(find(temp3<=nback+1)));
        temp2 = find(temp< HM_interp(n));
        locs_out = temp2(min(find(temp2>=nback+1)));
        FWHM_interp(n) = (locs_in(end)-locs_in(1)+1)*dt_mov/upsamplingrate;
    else
        peak_interp(n) = min(kernel_stack_interp(:,n));
        HM_interp(n) = 0.5*(base_interp(n)+peak_interp(n));
        temp = kernel_stack_interp(:,n);
        temp3 = find(temp > HM_interp(n));
        locs_in = temp3(max(find(temp3<=nback+1)));
        temp2 = find(temp> HM_interp(n));
        locs_out = temp2(min(find(temp2>=nback+1)));
        FWHM_interp(n) = (locs_out-locs_in(1)+1)*dt_mov/upsamplingrate;
    end
    
%     plot([locs_in(1) locs_in(end)+1],[temp(locs_in(1)) temp(locs_in(end)+1)])
    clear locs_in locs_out
    clear temp temp2
end
figure()
plot(FWHM_interp)
box off
axis tight
xlabel('Spike num')
ylabel('FWHM (ms)')
title(['FWHM calculation from GEVI, upsampling rate = ' ,num2str(upsamplingrate)])
saveas(gca,[pathname '\11 FWHM calculation_from GEVI.fig'])
saveas(gca,[pathname '\11 FWHM calculation_from GEVI.png'])
% kernel_stack_interp_m
ean = mean(kernel_stack_interp,2);
% plot(kernel_stack_interp)% for the check
figure()
subplot(2,1,1)
plot(FWHM_elec)
title(['FWHM calculation from patch clamp data, acquisition rate = ' ,num2str(Fs) ' Hz'])
box off
axis tight
xlabel('Spike num')
ylabel('FWHM (ms)')
subplot(2,1,2)
plot(FWHM_interp)
title(['FWHM calculation from GEVI, upsampling rate = ' ,num2str(upsamplingrate)])
box off
axis tight
xlabel('Spike num')
ylabel('FWHM (ms)')
saveas(gca,[pathname '\12 FWHM comparison.fig'])
saveas(gca,[pathname '\12 FWHM comparison.png'])

figure()
plot(FWHM_elec,FWHM_interp,'o')
box off
% axis tight
ylabel('FWHM from GEVI (ms)')
title(['FWHM display, acquisition rate = ', num2str(Fs,4), ' Hz'])
saveas(gca,[pathname '\13 FWHM 2D-display.fig'])
saveas(gca,[pathname '\13 FWHM 2D-display.png'])
%
% plot(kernel_stack_interp(:,n));hold on
% plot([locs_in(1) locs_in(end)+1],[temp(locs_in(1)) temp(locs_in(end)+1)]);hold on
% plot([1 1601],[HM_interp(1) HM_interp(1)])
% box off
% axis tight
%% Write a xlsx file
xlswrite([pathname 'analysis.xlsx'],{'Sensitivity (single trace)','S.E.M of sensitivity', 'baseline (brightness)', 'SNR from single trail','S.E.M of signle-trail SNR','spike number','Ra','Rm','Cm','path'}','Average','A1');
xlswrite([pathname 'analysis.xlsx'],sensitivity_single','Average','B1');
xlswrite([pathname 'analysis.xlsx'],std(sensitivity_vec)./sqrt(n),'Average','B2');
xlswrite([pathname 'analysis.xlsx'],Brightness','Average','B3');
xlswrite([pathname 'analysis.xlsx'],SNR_single','Average','B4');
xlswrite([pathname 'analysis.xlsx'],std(SNR_vec)./sqrt(n),'Average','B5');
xlswrite([pathname 'analysis.xlsx'],nspike','Average','B6');
xlswrite([pathname 'analysis.xlsx'],Ra,'Average','B7');
xlswrite([pathname 'analysis.xlsx'],Rm','Average','B8');
xlswrite([pathname 'analysis.xlsx'],Cm','Average','B9');
xlswrite([pathname 'analysis.xlsx'],{pathname},'Average','B10');
% xlswrite([pathname 'analysis.xlsx'],intens_mean','Average','I2');

xlswrite([pathname 'analysis.xlsx'],{'Spike num','Sensitivity ', 'SNR','Brightness','Noise (RMSE)','Delta F','FWHM from patch','FWHM from GEVI'},'single_trace','A1');
xlswrite([pathname 'analysis.xlsx'],[1:1:n]','single_trace','A2');
xlswrite([pathname 'analysis.xlsx'],sensitivity_vec','single_trace','B2');
% xlswrite([pathname 'analysis.xlsx'],baseline_vec','single_trace','B2');
xlswrite([pathname 'analysis.xlsx'],SNR_vec','single_trace','C2');
xlswrite([pathname 'analysis.xlsx'],baseline_vec','single_trace','D2');
xlswrite([pathname 'analysis.xlsx'],RMSE','single_trace','E2');
xlswrite([pathname 'analysis.xlsx'],dF_vec','single_trace','F2');
xlswrite([pathname 'analysis.xlsx'],FWHM_elec','single_trace','G2');
xlswrite([pathname 'analysis.xlsx'],FWHM_interp','single_trace','H2');
% xlswrite([pathname 'analysis.xlsx'],dF_vec','single_trace','D2');
%% generate a movie from the kernle we selected
mkdir([pathname '\movie\tiff'])
name1 = 'Sensor'
mov_APmovie = zeros(size(img,1),size(img,2),Lk,length(spikeT_real));
for n = 1:length(spikeT_real)
    mov_APmovie(:,:,:,n) = mov(:,:,spikeT_real(n)-(Lk-1)*0.5:spikeT_real(n)+(Lk-1)*0.5);
end
mov_APmovie = mean(mov_APmovie,4);
mov_APmovie = mov_APmovie(:,:,21:51);
t_APmovie = t_mov2(1:size(mov_APmovie,3));
clear M j
figure();
for j = 1:size(mov_APmovie,3);
    imshow(mov_APmovie(:,:,j),[410,800], 'Border', 'tight')
    colormap(pseudocolor(585));
    text(230,20,[num2str(t_APmovie(j),3) ' ms'], 'FontSize', 16, 'color', [0.99 0.99 0.99])
%     text(62,20,[' min'], 'FontSize', 20, 'color', [0.99 0.99 0.99])
    text(10,330,[name1], 'FontSize', 18, 'color', getrgb(585))
%     text(10,480,['180 uM GM1'], 'FontSize', 16, 'color', [0.99 0.99 0.99])
    M(j) = getframe(gca);
    %saveas(gca, [pathname '\Blink_' num2str((j)) '.tif']);  % Use this line to save each frame as a separate TIF file.
end;
% save the movie
speedup = 20;
myObj = VideoWriter([pathname '\movie'  '\' name1 ' movie_' num2str(speedup) 'x.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = speedup*(1/dt_mov)
open(myObj);
writeVideo(myObj,M)
close(myObj);
myObj = VideoWriter([pathname '\movie' '\' name1 ' movie_compressed_H.264 encoding_' num2str(speedup) 'x.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = speedup*(1/dt_mov);
open(myObj);
writeVideo(myObj,M)
close(myObj); 

save([pathname name1 ' .mat']);% save the variant
%% save a .tiff file
mkdir([pathname  '\movie\tiff\' num2str(name1) '_lapse\'])
Tiff([pathname  '\movie\tiff\' num2str(name1) '.tif'],'w')
for num_tiff=1:size(mov_APmovie,3)
    imwrite(uint16(mov_APmovie(:,:,num_tiff)),[pathname  '\movie\tiff\' num2str(name1) '.tif'],'WriteMode', 'append',  'Compression','none')
    imwrite(uint16(mov_APmovie(:,:,num_tiff)),[pathname  '\movie\tiff\' num2str(name1) '_lapse\' num2str(name1)  '_' num2str(num_tiff,2) '.tif'],'WriteMode', 'append',  'Compression','none')
end
%plot(squeeze(mean(mean(mov_APmovie))));
%% Generate the elec. mean trace
mean_elec_kernel = mean(kernel_stack_elec,2);
temp = mean_elec_kernel((21-1)*upsamplingrate+1:(62-1)*upsamplingrate);
mkdir([pathname '\movie\for GIF\'])
clear K
% set the background of the figure for printing screen
%set(0,'defaultfigurecolor','k') 
%colordef black
%set(0,'defaultfigurecolor','w') 
%colordef white
for nn = 1:820;
figure('color','k')
set(gca, 'color', [0 0 0])

plot([0:2.0658/20:819*dt_mov/upsamplingrate],temp,'w','LineWidth',1.5);hold on
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')
box off
axis tight
plot((nn-1)*dt_mov/upsamplingrate,temp(nn),'r+','LineWidth',10);
% saveas(gcf,[pathname '\movie\for GIF\' num2str(nn) '.fig'])
% saveas(gcf,[pathname '\movie\for GIF\' num2str(nn) '.bmp'])
% saveas(gcf,[pathname '\movie\for GIF\' num2str(nn) '.png'])
% saveas(gcf,[pathname '\movie\for GIF\' num2str(nn) '.jpg'])
K(nn) = getframe(gcf);

close(gcf)
end
speedup = 400;
myObj = VideoWriter([pathname '\movie'  '\' name1 ' test_' num2str(speedup) 'x.avi'],'Uncompressed AVI');%create an AVI file
myObj.FrameRate = speedup*(1/dt_mov)
open(myObj);
writeVideo(myObj,K)
close(myObj);
myObj = VideoWriter([pathname '\movie' '\' name1 ' test_compressed_H.264 encoding_' num2str(speedup) 'x.mp4'],'MPEG-4');%create an AVI file
myObj.FrameRate = speedup*(1/dt_mov);
open(myObj);
writeVideo(myObj,K)
close(myObj); 

% plot([1:20:801],temp([1:20:801]),'.')
% reset the background after printing screen
%set(0,'defaultfigurecolor','w') 
%colordef white

