% parameters
subfolder = '';
basepath = 'D:\data\test\20220401 Depolarization titration\Dish2\cell3\202733_40xWideField_bin2\';
pathname = basepath
path = [basepath subfolder];

% constants
movname = '\movie.bin';
ncol = 176;         % x
nrow = 96;         % y
bkg = 400;          % background due to camera bias (100 for bin 1x1)
dt_mov = 0.9452;       % exposure time in millisecond
DAQname = '\movie_DAQ.txt';
dnsamp = 20;       % downsampling rate = DAQ rate
dt_daq = dt_mov/dnsamp;     % in millisecond

% load movie
fname = [path movname]
[mov, nframe] = readBinMov(fname, nrow, ncol);
mov = double(mov);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov*0.001;     % time axis in second

% select ROI for analysis
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
% save clicky figure
saveas(gca,[path '\clicky analysis_gap junction.fig']);
saveas(gca,[path '\clicky analysis_gap junction.png']);
% normalize intensity wrt -100 mV (in %)
bkg=mean(intens(:,end));
background_reference=sort(reshape(mean(mov,3),ncol*nrow,1));
background_reference=mean(background_reference(10:round(0.05*length(background_reference))))
intens_norm = (intens-background_reference)/(max(intens)-background_reference)*100;
% intens_norm = intens_norm./(intens_norm(1,1));
% load DAQ data
tmp = importdata([path DAQname]);   % import data 
data = tmp.data;                    % get array
Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
len = length(Vm);
t_daq = [0:(len-1)]*dt_daq./1000;                 % time axis in second
% Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));
intens_rembkg = intens(:,1:end-1)-intens(:,end);
%intens_rembkg = intens;

%%
% plot traces
figure;
% plot F-t
subplot(2,1,1);
itens_size = size(intens)
plot(t_mov, smooth(smooth(intens_rembkg)));
% for i = 1: itens_size(2)
%     plot(t_mov, smooth(smooth(intens_rembkg(:,i))));
%     hold on;
% end
xlabel('time (s)')';
ylabel('Intens (a.u.)');
axis tight;
box off
% plot Vm-t
subplot(2,1,2);
plot(t_daq,Vm);
xlabel('time (s)');
ylabel('V_m (mV)');
axis tight;
box off
% save F-V figure
saveas(gca,[path '\F-V analysis.fig']);
saveas(gca,[path '\F-V analysis.png']);
%%
color_map = hsv(11)
color_map = flipud(color_map)
intens_stack = zeros(1000,11);
for n = 1:11
  intens_stack(:,n) = intens_rembkg((251+1000*(n-1):250+1000*n));
end;
figure();
for n = 1:11
plot([0:dt_mov:(1000-1)*dt_mov]', intens_stack(:,n),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Intensity')
saveas(gca,[path '\F-V stack.fig']);
saveas(gca,[path '\F-V stack.png']);

figure();
for n = 1:11
plot(smooth(smooth([0:dt_mov:(1000-1)*dt_mov]',1-intens_stack(:,n)./mean(intens_stack(1:ceil(100/dt_mov),n)))),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Normalized Intensity (to -70 mV)')
saveas(gca,[path '\Norm.F-V stack.fig']);
saveas(gca,[path '\Norm.F-V stack.png']);

%% for the colorbar
figure();
temp = [-100:(200/10):100]';
temp = repmat(temp',3,1);
imshow(temp',[])
colormap(hsv(11))
saveas(gca,[path '\colorbar.fig']);
saveas(gca,[path '\colorbar.png']);

%%
%pathname = 'X:\91 Data and analysis\YJunqi\Sensitivity\Screening for AceD81X-HaloTag\20190501 HEK293T cells\Dish6\cell2\191806_slow step\';
intens_plat = zeros(91,11);
intens_base = zeros(91,11);
norm_inten_plat = zeros(1,11);
for n = 1:11
  intens_satack_single = intens_stack(:,n);
  intens_plat(:,n) = intens_satack_single(307:397);
  intens_base(:,n) = intens_satack_single(100:190);
end;
norm_inten_plat = 1 - mean(intens_plat,1)./mean(intens_base,1);

figure();
for n = 1:11
plot([0:dt_mov:(90)*dt_mov]',intens_plat(:,n),'color',color_map(n,:));hold on
end
box off
axis tight
xlabel('Time (ms)')
ylabel('Intensity')
saveas(gca,[path '\F-V stack_peak.fig']);
saveas(gca,[path '\F-V stack_peak.png']);

xlswrite([pathname '0analysis.xlsx'],{'Vm', 'delta F/F0'},'A1:B1');
xlswrite([pathname '0analysis.xlsx'],{-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100}','A2:A12')
xlswrite([pathname '0analysis.xlsx'],norm_inten_plat','B2:B12')


