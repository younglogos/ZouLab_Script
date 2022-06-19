% parameters
basepath = 'E:\ZouLab\YJunqi\Sensitivity\20220606 Kinetics, sensitivity & depolarization test\Dish2\cell1\203925_40xWideField_bin2';
subfolder = '';
path = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 312;         % x
nrow = 208;         % y
bkg = 400;          % background due to camera bias (100 for bin 1x1)
dt_mov = 0.1;       % exposure time in second
DAQname = '\movie_DAQ.txt';
dnsamp = 100;       % downsampling rate = DAQ rate
dt_daq = 0.001;     % in millisecond
ratio = 100./dnsamp;


% load movie
fname = [path movname]
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = double(mov);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov;     % time axis in second

% select ROI for analysis
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
% save clicky figure
saveas(gca,[path '\clicky analysis_gap junction.fig']);
saveas(gca,[path '\clicky analysis_gap junction.png']);
% normalize intensity wrt -100 mV (in %)
bkg=mean(intens(:,2));
background_reference=sort(reshape(mean(mov,3),ncol*nrow,1));
background_reference=mean(background_reference(10:round(0.05*length(background_reference))))
intens_norm = (intens-background_reference)/(max(intens)-background_reference)*100;
% load DAQ data
tmp = importdata([path DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
len = length(Vm);
t_daq = [0:(len-1)]*dt_daq;                 % time axis in millisecond
Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));

% plot traces
figure;
% plot F-t
subplot(2,2,1);
plot(t_mov,intens);
xlabel('time (s)');
ylabel('Intens (a.u.)');
axis tight;
% plot Vm-t
subplot(2,2,3);
plot(t_daq,Vm);
ylabel('V_m (mV)');
axis tight;
% plot F-V curve
subplot(2,2,2);
plot(Vm_dnsamp, intens);
xlabel('V_m (mV)');
ylabel('Intensity');
axis tight;
% plot normalized \DeltaF/F-V curve
subplot(2,2,4);
plot(Vm_dnsamp, intens_norm);
title('Normalized F-V curve');
xlabel('V_m (mV)');
ylabel('\DeltaF/F (%)');
axis tight;

% save F-V figure
saveas(gca,[path '\F-V analysis.fig']);
saveas(gca,[path '\F-V analysis.png']);
%%
V_trace = intens(:,1)-intens(:,2);
subplot(2,1,1)
plot(t_daq,Vm)
box off
ylim([-100 50])
xlabel('time (s)')
ylabel('Vm (mV)')
subplot(2,1,2)

plot(t_mov,V_trace)
box off
axis tight
xlabel('time (s)')
ylabel('Intensity (W/O background)')

param_pbrem = V_trace(ratio*[28 68 88 128 148 188 208])
param_pbrem = param_pbrem./max(param_pbrem)
p = polyfit(ratio*[28 68 88 128 148 188 208],param_pbrem',1)
plot(ratio*[28 68 88 128 148 188 208],param_pbrem')
correct_line = [1:p(1):1+p(1)*(size(V_trace,1)-1)]';
V_rem = V_trace./correct_line;
plot(V_rem)

subplot(2,1,1)
plot(t_mov,V_trace)
box off
axis tight
xlabel('time (s)')
ylabel('Intensity (W/O background)')
title('Before removing PB')
subplot(2,1,2)
plot(t_mov,V_rem)
box off
axis tight
xlabel('time (s)')
ylabel('Intensity (W/O background)')
title('After removing PB')
saveas(gca,[path '\F-V analysis_rem pb.fig']);
saveas(gca,[path '\F-V analysis_rem pb.png']);

figure()
% plot(Vm_dnsamp, V_rem./max(V_rem)*100-100);
plot(Vm_dnsamp, V_rem./min(V_rem)*100-100);
title('Normalized F-V curve');
xlabel('V_m (mV)');
ylabel('\DeltaF/F (%)');
axis tight;
box off
axis tight
saveas(gca,[path '\F-V analysis_rem pb_stack.fig']);
saveas(gca,[path '\F-V analysis_rem pb_stack.png']);

%%
% temp_y = V_rem./max(V_rem)*100;
temp_y = V_rem./min(V_rem)*100;
% temp_y = temp_y(46*ratio:225*ratio);
temp_y = temp_y(50*ratio:229*ratio);
average_trace = mean(reshape(temp_y,60*ratio,3),2);  
average_trace = average_trace./average_trace(1);

temp_x = Vm_dnsamp;
% temp_x = temp_x(46*ratio:225*ratio);
temp_x = temp_x(50*ratio:229*ratio);
average_trace_Vm = mean(reshape(temp_x,60*ratio,3),2);
figure()
plot(average_trace_Vm(1:30*ratio+1),100*average_trace(1:30*ratio+1)-100,'color',getrgb(450),'LineWidth',2.25);hold on
plot(average_trace_Vm(30*ratio+1:60*ratio),100*average_trace(30*ratio+1:60*ratio)-100,'color',getrgb(640),'LineWidth',2.25)
legend('Depolarization','Hyperpolarization')
legend('boxoff')         
% title(['Normalized F-V curve, -100 to 50 mV: \DeltaF/F (%) = ' num2str(100*average_trace(30*ratio+1)-100) ]);
title(['Normalized F-V curve, -100 to 100 mV: \DeltaF/F (%) = ' num2str(100*average_trace(30*ratio+1)-100) ]);
xlabel('V_m (mV)');
ylabel('\DeltaF/F (%)');
axis tight;
box off
box off
axis tight
saveas(gca,[path '\F-V analysis_rem pb_average_trace.fig']);
saveas(gca,[path '\F-V analysis_rem average_trace.png']);