% analysis of photocurrent in voltage clamp
% for HVI, Ace2N-2AA-mNeon and other GEVIs,
% eg, fig s2d,s2e
clear all
%%
clear all;
basepath = 'H:\ZouLab\YJunqi\Photophysical property\Cepheid\20220826 Photo-current test\Dish2\cell1\181924_1V_561nmImaging+5V_488nmPulse'
subfolder = ''
WL = 561;%
pathname = [basepath subfolder '\']
% name1 = 'Ace D81S-mOrange2';
% Quad = 'O';
a = importdata([pathname '\movie_DAQ.txt']);
% b = importdata([basepath '\patch param.txt']);
data=a.data;
% datab = b.data;
AI_scaled=data(:,1);
AI_10Vm=data(:,2)*100;
% Cm = datab(1,5);%pF
Fs = 21159.48;
time=(1:length(AI_scaled)')'./Fs;
% figure;
% set(gcf,'outerposition',get(0,'screensize'));
% plot(time,AI_scaled,time,AI_10Vm);
% legend('AI\_scaled','Vm(mV)','Location','Northeast');
% hold on
% xlim=[0,max(time)]
% ylim=[-80,70];
% xL=xlim;yL=ylim;
% set(gca,'xtick',[0:5:max(time)])
% plot(xL,[yL(2),yL(2)],'w',[xL(2),xL(2)],[yL(1),yL(2)],'w')
% box off
% axis([xL yL]) 
% axis tight
% saveas(gca,[pathname 'waveform.fig']);
% saveas(gca,[pathname 'waveform.png']);


%%
load ([pathname '\20KHz_300ms_on+700ms_off_full_power\matlab variables.mat']);
timeunit = hiPts+lowPts;
AO_PC = AO(1:size(time,1));
DO_PC = DO(1:size(time,1));
figure();
subplot(3,1,1)
plot(time,AO_PC);
xlabel('Time (s)')
ylabel('Laser intensity (AO) (V)')
box off
axis tight
subplot(3,1,2)
plot(time,DO_PC);
xlabel('Time (s)')
ylabel('Laser on/off (DO) (V)')
box off
axis tight
subplot(3,1,3)
plot(time,AI_scaled*1000);
xlabel('Time (s)')
ylabel('Photo current(pA)')
box off
axis tight
saveas(gca,[pathname 'global_photocurrent.fig']);
saveas(gca,[pathname 'global_photocurrent.png']);
%% plot
% figure();
% subplot(3,1,1)
% plot(time(0.5*timeunit+1:5.5*timeunit),AO_PC(0.5*timeunit+1:5.5*timeunit));
% xlabel('Time (s)')
% ylabel('Laser intensity (AO) (V)')
% box off
% axis tight
% subplot(3,1,2)
% plot(time(0.5*timeunit+1:5.5*timeunit),DO_PC(0.5*timeunit+1:5.5*timeunit));
% xlabel('Time (s)')
% ylabel('Laser on/off (DO) (V)')
% box off
% axis tight
% subplot(3,1,3)
% plot(time(0.5*timeunit+1:5.5*timeunit),AI_scaled(0.5*timeunit+1:5.5*timeunit)*1000);
% xlabel('Time (s)')
% ylabel('Photo current(pA)')
% box off
% axis tight
% saveas(gca,[pathname 'global_photocurrent_0.4V.fig']);
% saveas(gca,[pathname 'global_photocurrent_0.4V.png']);
kernel_stack = reshape(AI_scaled(0.8*timeunit+1:(cycles+0.8)*timeunit)*1000,timeunit,cycles);
% kernel = mean(kernel_stack,2);
% figure()
% plot([0:1/Fs:(timeunit-1)/Fs],kernel);
% box off
% axis tight
% xlabel('Time (s)')
% ylabel('Photo current(pA)')
% saveas(gca,[pathname 'Mean.fig']);
% saveas(gca,[pathname 'Mean.png']);

figure()
plot([0:1/Fs:(timeunit-1)/Fs],kernel_stack);
box off
axis tight
xlabel('Time (s)')
ylabel('Photo current(pA)')
title('stack: from the period')
saveas(gca,[pathname 'v2Mean_stack from period.fig']);
saveas(gca,[pathname 'v2Mean_stack from period.png']);

kernel = mean(kernel_stack,2);
plot(kernel)

figure()
plot([0:1/Fs:(timeunit-1)/Fs],kernel)
box off
axis tight
xlabel('Time (s)')
ylabel('Photo current(pA)')
title(['photocurrent with (' num2str(WL) ' nm,' num2str(hiPts/Fs,3) ' s illumination)'])
saveas(gca,[pathname 'v2Mean.fig']);
saveas(gca,[pathname 'v2Mean.png']);

figure()
plot([0:1/Fs:(timeunit-1)/Fs],kernel-mean(kernel(1:900)))
box off
axis tight
xlabel('Time (s)')
ylabel('Photo current(pA)')
title(['photocurrent with (' num2str(WL) ' nm,' num2str(hiPts/Fs,3) ' s illumination)'])
saveas(gca,[pathname 'v2Mean_move to 0.fig']);
saveas(gca,[pathname 'v2Mean_move to 0.png']);
%% If you want to generate the kernel by peak finding, run this knot
[~,spikeTon] = findpeaks(AI_scaled*1000,'MinPeakDistance',timeunit-100); % for outward PC
% next line for inward pc
dire = -1;% if inward (like ChR2), dire = -1; if outward (like Arch), dire = 1; 
[~,spikeTon1] = findpeaks(dire*(AI_scaled),'MinPeakDistance',timeunit-100);
findpeaks(-1*(AI_scaled),'MinPeakDistance',timeunit-100);hold on
title('Left click one time for the threshold');
[~,threshold_PC] = ginput(1);
plot([1 length(AI_scaled)],[threshold_PC threshold_PC],'r','LineWidth',2);hold on
legend('Threshold','Peak from period')
[~, spikeTon2] = findpeaks(dire*(AI_scaled), 'MinPeakHeight',threshold_PC);% last 500 points are excluded now for better analysis

spikeTon = intersect(spikeTon1,spikeTon2); % find the intersection of the peaks location





kernel_stack_align = zeros(timeunit,size(spikeTon,1));
for n = 1:size(spikeTon,1)
    kernel_stack_align(:,n) = AI_scaled(spikeTon(n)-round(timeunit*0.05)+1:spikeTon(n)+timeunit-round(timeunit*0.05));
end
kernel = mean(kernel_stack_align,2)*1000;
figure()
plot([0:1/Fs:(timeunit-1)/Fs],kernel_stack_align*1000)
box off
axis tight
xlabel('Time (s)')
ylabel('Photo current(pA)')
title('stack: align the peak on')
saveas(gca,[pathname 'v2stack from alignment.fig']);
saveas(gca,[pathname 'v2stack from alignment.png']);
figure()
plot([0:1/Fs:(timeunit-1)/Fs],kernel)
box off
axis tight
xlabel('Time (s)')
ylabel('Photo current(pA)')
title(['photocurrent with (' num2str(WL) ' nm,' num2str(hiPts/Fs,3) ' s illumination)'])
saveas(gca,[pathname 'v2Mean.fig']);
saveas(gca,[pathname 'v2Mean.png']);

% average
% meanPC = [kernel' PC_step2' PC_step3' PC_step4' PC_step5' PC_step6'];
% base = mean(mean(kernel(1:2000)));
% kernel = kernel-base;
% plot(kernel);
% figure()
% plot([0:1/Fs:(6*timeunit-1)/Fs],kernel);
% title('photocurrent of Ace(D81S)-mOrange2 (532 nm,198 ms illumination) ')
% box off
% axis tight
% xlabel('Time (s)')
% ylabel('Photo current(pA)')
% saveas(gca,[pathname 'Total.fig']);
% saveas(gca,[pathname 'Total.png']);
%% 2.8V, 4.0V and 4.8V with light
% 
% figure()
% plot([0:1/Fs:(timeunit-1)/Fs],kernel-base);hold on;
% plot([timeunit/2/Fs (timeunit/2+4200)/Fs],[max(kernel)+3-base max(kernel)+3-base],'Color',getrgb(WL),'LineWidth',5)
% box off
% axis tight
% xlabel('Time (s)')
% ylabel('Photo current(pA)')
% title(['photocurrent of Ace(D81S)-mOrange2 (' num2str(WL) ' nm,198 ms illumination)'])
% saveas(gca,[pathname 'Mean_' num2str(WL) '.fig']);
% saveas(gca,[pathname 'Mean_' num2str(WL) '.png']);
% close(gcf)

%% calculation
% base_low(1) = mean(kernel(1:round(timeunit/6)));
base_low(1) = mean(kernel(1:900));
figure()
plot(kernel)
title('Select the region of steady photocurrent')
[b0] = ginput(2);
b0 = round(b0);
stable(1) = mean(kernel(b0(1):b0(2)));%cuole
stable_pc(1) = stable(1)-base_low(1);
peak_on(1) = min(kernel)-base_low(1);
peak_off(1) = max(kernel)-base_low(1);
peakloc = find(kernel==min(kernel));
space = kernel(peakloc-1000:peakloc+1000); % choose the suitable kernel range of X axis
figure()
plot(space)
title('select two point as the boundary of current_ on')
[b1,~] = ginput(2);
peakloc = find(kernel==max(kernel));
width_on = (b1(2)-b1(1))./Fs*1000
c1 = min(find(space <= (base_low+peak_on)/2));
d1 = max(find(space <= (base_low+peak_on)/2)); 
FHWM_on = (d1-c1)/Fs*1000;
space = kernel(peakloc-2000:peakloc+10000); % choose the suitable kernel range of X axis
figure()
plot(space)
title('select two point as the boundary of current off')
[b2,~] = ginput(2);
width_off = (b2(2)-b2(1))./Fs*1000
c2 = min(find(space >= (base_low+peak_off)/2));
d2 = max(find(space >= (base_low+peak_off)/2)); 
FHWM_off = (d2-c2)/Fs*1000


xlswrite([pathname 'v2analysis.xlsx'],{'Laser (V)','base (pA)', 'stable (pA)', 'stable_pc (pA)','peak_on (pA)', 'peak off (pA)','width_on (ms)','FHWM_on (ms)','width_off (ms)','FHWM_off (ms)'},'photocurrent','A1');
xlswrite([pathname 'v2analysis.xlsx'],max(AO(:))','photocurrent','A2');
xlswrite([pathname 'v2analysis.xlsx'],base_low','photocurrent','B2');
xlswrite([pathname 'v2analysis.xlsx'],stable','photocurrent','C2');
xlswrite([pathname 'v2analysis.xlsx'],stable_pc','photocurrent','D2');
xlswrite([pathname 'v2analysis.xlsx'],peak_on','photocurrent','E2');
xlswrite([pathname 'v2analysis.xlsx'],peak_off','photocurrent','F2');

% xlswrite([pathname 'v2analysis.xlsx'],width_on','photocurrent','G2');
% xlswrite([pathname 'v2analysis.xlsx'],FHWM_on','photocurrent','H2');
% xlswrite([pathname 'v2analysis.xlsx'],width_off','photocurrent','I2');
% xlswrite([pathname 'v2analysis.xlsx'],FHWM_off','photocurrent','J2');
