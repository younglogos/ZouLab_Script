name1 = '90%power532nm_250Hz_Data00085';
name1color = getrgb(550);
bin = 4;
basepath = 'H:\ZouLab\YJunqi\Application\Mouse\20220428 Aureole-ST in acute brain slice\';
subfolder = '';
pathname = [basepath subfolder];
%mkdir([pathname 'analysis']); 
FileTif=[pathname name1 '.tif'];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov=zeros(nImage,mImage,NumberImages,'uint16');
start_frame = 1;

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   mov(:,:,i)=TifLink.read();
end
TifLink.close();
%bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
dt_mov = 4;         % exposure time in millisecond
mov = double(mov);
mov = mov(:,:,start_frame:end);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov/1000;     % time axis in second
ncol = size(mov,2);         % x
nrow = size(mov,1);          % y
intens_glo = squeeze(mean(mean(mov)));
%intens_glo = intens_glo-bkg;
figure();
intens_glo_norm = intens_glo./intens_glo(1);
plot(t_mov,intens_glo_norm);
box off
axis tight
xlabel('Time (s)');
ylabel('Intensity');
saveas(gca,[pathname '\global voltage intensity.fig']);
saveas(gca,[pathname '\global voltage intensity.png']);
%  mov = mov-15000;
ColorMatrix = get(gca,'colororder');
% mov = mov(:,30:727,:);
img = mean(mov, 3);

% select ROI for analysis
[~, intens] = clicky(mov, img, 'Select only 1 ROI, right click when done');
% save clicky figure
saveas(gca,[pathname '\clicky analysis_voltage indicator.fig']);
saveas(gca,[pathname '\clicky analysis_voltage indicator.png']);
%saveas(gca,[path '\clicky analysis_calicium indicator.fig']);
%saveas(gca,[path '\clicky analysis_calicium indicator.png']);
% normalize intensity (in %)



%% Voltage imaging
intensN = zeros(size(intens,1),size(intens,2)-1);
pbleach = zeros(size(intens,1),size(intens,2)-1);
bkg=intens(:,end);
for i = 1: size(intens,2)-1
    [intensN(:,i), pbleach(:,i)] = rem_pbleach(-1*(intens(:,i)-bkg), round(size(intens(:,i),1)/5)+1);
    %intensN1(:,i) = intens(:,i) - bkg + power(pbleach1(:,i) - min(pbleach1(:,i)),0.5);
    %pbleach(:,i) = -1*power(abs(pbleach(:,i)+min(pbleach(:,i))),1);
    %[intensN(:,i), ~] = rem_pbleach(intensN1(:,i), round(size(intensN1(:,i)/50000,1)));
    %intensN(:,i) = intensN(:,i) + power(pbleach(:,i) - min(pbleach(:,i)),1);
    %intensN(:,i) = intensN(:,i)./power(abs(pbleach(:,i)+max(pbleach(:,i))),1);
end
intens_base = sort(-1*intensN, 1);
intens_base = mean(-1*intens_base(1:500,:))
intens_norm = (intensN-repmat(intens_base, size(intensN,1),1))./repmat(intens_base, size(intensN,1),1)*100;
% intens_norm = intens_norm./(intens_norm(1,1));
% load DAQ data
%tmp = importdata([path DAQname]);   % import data 
%data = tmp.data;                    % get array
%Vm = data(:,2)*100;     % Vm in millivolt, column vector
% downsample Vm
%len = length(Vm);
%t_daq = [0:(len-1)]*dt_daq./1000;                 % time axis in second
%Vm_dnsamp = mean(reshape(Vm,dnsamp,len/dnsamp));
%intens_rembkg = intens(:,1)-intens(:,2);
%intens_rembkg = intens;

% plot traces
figure;
% plot F-t      
%subplot(2,1,1);
itens_size = size(intens)
%plot(t_mov, smooth(smooth(intens_rembkg)));
for i = 1: itens_size(2)-1
    plot(t_mov, smooth(smooth(intens_norm(:,i))));
    hold on;
    
end
xlabel('time (s)');
ylabel('¦¤F/F0 (%)'); 
axis tight;
box off
saveas(gca,[pathname '\Voltage imaging.fig']);
saveas(gca,[pathname '\Voltage imaging.png']);