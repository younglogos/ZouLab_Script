% Writed by Junqi Yang, 2021-7-22
clear all
%% Load imaging data
basepath = 'X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\20210710 KCl stimulus pre-experiment\Group1\Dish1\';
subfolder = '';
filename = 'dual-view_ArcLight&HASAP1';
indicator1 = 'HASAP1-JF635';
indicator2 = 'ArcLight';
pathname = [basepath subfolder];
FileTif=[pathname name '.tif'];
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

dt_mov = 10000;         % frame period in millisecond
mov = double(mov);
mov = mov(:,:,start_frame:end);
img = mean(mov, 3);
len = size(mov,3);
t_mov = [0:(len-1)]*dt_mov/1000; 

%% Seleect ROI for analysis
mkdir([pathname 'analysis']);
[~, intens] = clicky_contrastAdjust(mov, img, 'Select only 1 ROI, right click when done', 500, 2000);
saveas(gca,[pathname 'analysis\selected regions of ' name '.fig']);
saveas(gca,[pathname 'analysis\selected regions of ' name '.png']);
right_fov = zeros(size(intens, 1), size(intens, 2)/2);
left_fov = zeros(size(intens, 1), size(intens, 2)/2);

for i = 1:size(intens, 2)/2
    right_fov(:, i) = intens(:,i*2-1);
    left_fov(:, i) = intens(:,i*2);
end
channel1 = right_fov(:, 1:end-1) - right_fov(:, end);
channel2 = left_fov(:, 1:end-1) - right_fov(:, end);
%%
figure()
for i = 1:size(channel1, 2)
    plot(t_mov, channel1(:, i), 'Color', getrgb(400+round((i-1)/size(channel1,2)*270))); hold on;
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
saveas(gca,[pathname 'analysis\' indicator1 '.fig']);
saveas(gca,[pathname 'analysis\' indicator1 '.png']);

figure()
for i = 1:size(channel2, 2)
    plot(t_mov, channel2(:, i), 'Color', getrgb(400+round((i-1)/size(channel2,2)*270))); hold on;
end
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 1','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
saveas(gca,[pathname 'analysis\' indicator2 '.fig']);
saveas(gca,[pathname 'analysis\' indicator2 '.png']);

figure()
deltaNorm_channel1 = channel1./mean(channel1(1:10, :), 1) - 1;
deltaNorm_channel2 = channel2./mean(channel2(1:10, :), 1) - 1;
plot(t_mov, deltaNorm_channel1, 'r'); hold on;
plot(t_mov, deltaNorm_channel2, 'g');

ratio_intens = zeros(size(channel1));
for i = 1:size(ratio_intens, 2)
    ratio_intens(:, i) = abs(deltaNorm_channel1(:, i))-abs(deltaNorm_channel2(:, i));
    plot(t_mov, ratio_intens(:, i), 'Color', getrgb(400+round((i-1)/size(ratio_intens,2)*270))); hold on;
end

