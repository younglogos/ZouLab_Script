%% load and select
clear all; clc;
% Movie loading path
basepath = 'X:\91 Data and analysis\YJunqi\Cytotoxcity test\20210918 Cytotoxcity&Sensitvity of AceS81-HaloTag(HaloTag-ST2) in neurons\Dish2\cell3\';
subfolder = '220642_I = 0';
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

dire = 1;% dire = 1 means positive GEVI; dire = -1 means negative GEVI
%% I = 0 mode
load ([pathname 'Current injection_200pA-10ms-80 times_(Baseline 0 pA)\matlab variables.mat']);
Fs = samprate;

DAQname = '\movie_DAQ.txt';
% dnsamp = Fs/(1000/dt_mov);        % downsampling rate = DAQ rate/camera rate
% dnsamp = round(dnsamp);

% load DAQ data
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
% dt_daq = dt_mov/dnsamp;             % DAQ dt in millisecond
% t_daq = [0:length(Vm)-1]*dt_daq/10^3;       % DAQ time axis in second
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
ylim=[-80,-30];
xL=xlim;yL=ylim;
set(gca,'xtick',[0:5:max(time)]);
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
axis([xL yL]);
axis tight
saveas(gca,[pathname '0 waveform of AI.fig']);
saveas(gca,[pathname '0 waveform of AI.png']);

Vm = mean(AI_10Vm);
Std = std(AI_10Vm);
xlswrite([pathname 'electro_analysis.xlsx'],{'Mean resting potential (mV)','Std resting potential', 'Ra (MOhm)','Rm (MOhm)','Cm (pF)'},'I=0','A1');
xlswrite([pathname 'electro_analysis.xlsx'],Vm,'I=0','A2');
xlswrite([pathname 'electro_analysis.xlsx'],Std,'I=0','B2');
xlswrite([pathname 'electro_analysis.xlsx'],Ra,'I=0','C2');
xlswrite([pathname 'electro_analysis.xlsx'],Rm','I=0','D2');
xlswrite([pathname 'electro_analysis.xlsx'],Cm','I=0','E2');

%% Current clamp mode
load ([pathname 'Current injection_200pA-10ms-80 times_(Baseline 0 pA)\matlab variables.mat']);
Fs = samprate;

DAQname = '\movie_DAQ.txt';


% load DAQ data
% load DAQ data
tmp = importdata([pathname DAQname]);   % import data
data = tmp.data;                    % get array
Vm = data(:,2)*100;                 % Vm in millivolt, column vector
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
ylim=[-80,-30];
xL=xlim;yL=ylim;
set(gca,'xtick',[0:5:max(time)]);
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
box off
axis([xL yL]);
axis tight
saveas(gca,[pathname '0 waveform of AI.fig']);
saveas(gca,[pathname '0 waveform of AI.png']);

[boundary,~] = ginput(2);
boundary = boundary .* Fs;
boundary(boundary<=1) = 1 ;
boundary(boundary>=size(AI_10Vm, 1)) = size(AI_10Vm,1);
boundary = round(boundary);
period_param = 0.9;
%% Check the FWHM of AP with electrophysiological signal
% dnsamp = Fs/(1000/2.0658); 
nback = 800;  % number of points to extend kernel behind AP peak
nfront = 800; % number of points to extend kernel in front of AP peak

% The spikefind_corr algorithm finds spikes in two iterations.  First it
% uses a simple threshold-and-max procedure (spikefind.m).  It generates a
% kernel which it convolves with the original data.  The maxima in the
% convolved data are taken as the revised spike times.  This procedure
% identifies spikes by their global similarity to a template, rather than
% by the location of a single peak, which could be subject to noise.
if boundary(2) > length(AI_10Vm)
    AI_10Vm_interest = AI_10Vm(boundary(1):end);
    time_interest = time(boundary(1):end);
else
    AI_10Vm_interest = AI_10Vm(boundary(1):boundary(2));
    time_interest = time(boundary(1):boundary(2));
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
% kernel = kernel1;spikeT = spikeT1;
nspike_elec = length(spikeT_real_elec);
Lk_elec = nfront+nback+1;
kernel_stack_elec = zeros(Lk_elec,size(spikeT_real_elec,1));
for n = 1:size(spikeT_real_elec,1)
    kernel_stack_elec(:,n) = AI_10Vm_interest(spikeT_real_elec(n)-nfront:spikeT_real_elec(n)+nback);
    base_elec(n) = mean(kernel_stack_elec(1:600,n));
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

%% Drawing mean sitmulated potential with every stimulus
figure()
dt_daq = 1/Fs * 1000;
plot([0:1:size(kernel_stack_elec,1)-1]'*dt_daq,kernel_stack_elec,'Color',[0.8,0.8,0.8],'LineWidth',0.8);hold on;
plot([0:1:size(kernel_stack_elec,1)-1]'*dt_daq,mean(kernel_stack_elec,2),'Color',[0,0,0],'LineWidth',3);
box off
axis([0 dt_daq*size(kernel_stack_elec,1) -70 60])
axis tight
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')
title(['Average AP electrophysiological signal (Vm) from ' num2str(size(kernel_stack_elec,2)) ' APs'])
saveas(gca,[pathname '\3_1 average AP patch signal.fig'])
saveas(gca,[pathname '\3_1 average AP patch signal.png'])
%% Stimulated potenital chracterization (from electro. data)
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
        diff1_raw(:,n) = diff(Vm(:,n),1);
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
xlswrite([pathname 'electro_analysis.xlsx'],{'Resting potential (mV)','Spike number','Amplitude (mV)', 'Base (mV)', 'Peak voltage (mV)','Max. rising speed (mV/ms)','Vm (Max. rising speed) ','Threshold (mV)','FWHM (ms)','Ra (MOhm)','Rm (MOhm)','Cm (pF)'},'Single stimulus','A1');
xlswrite([pathname 'electro_analysis.xlsx'],Rp','Single Stimulus','A2');
xlswrite([pathname 'electro_analysis.xlsx'],[1:1:n]','Single Stimulus','B2');
xlswrite([pathname 'electro_analysis.xlsx'],Amplitude_elec','Single Stimulus','C2');
xlswrite([pathname 'electro_analysis.xlsx'],base_elec','Single Stimulus','D2');
xlswrite([pathname 'electro_analysis.xlsx'],peak_elec','Single Stimulus','E2');
xlswrite([pathname 'electro_analysis.xlsx'],rising_speed','Single Stimulus','F2');
xlswrite([pathname 'electro_analysis.xlsx'],rising','Single Stimulus','G2');
xlswrite([pathname 'electro_analysis.xlsx'],threshold','Single Stimulus','H2');
xlswrite([pathname 'electro_analysis.xlsx'],FWHM_elec','Single Stimulus','I2');
xlswrite([pathname 'electro_analysis.xlsx'],Ra,'Single Stimulus','J2');
xlswrite([pathname 'electro_analysis.xlsx'],Rm','Single Stimulus','K2');
xlswrite([pathname 'electro_analysis.xlsx'],Cm','Single Stimulus','L2');

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