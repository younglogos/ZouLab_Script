
%% ª≠≥ˆŒÛ≤Ó“ı”∞°£

figure()
SizeIntens = size(StepIntensMat)
NormIntensMat = 1- StepIntensMat./repmat(StepIntensMat(:,1), [1 SizeIntens(2)])
NormIntens = mean(NormIntensMat)
SdIntens = std(NormIntensMat);
SdShadow = area(t_step',[(NormIntens - SdIntens)', (2 * SdIntens)']);hold on;
%h = area(x, 2 * e);hold on;
%h = area(x, y - e);hold on;

set(SdShadow(1),'Visible','off');
set(SdShadow(2),'EdgeColor','white','FaceColor',[0.7,0.7,1]);
plot(t_step, NormIntens, '-b', 'LineWidth', 2);

xlabel('Time (ms)');
%ylim([0,1]);
ylabel('-¶§F/F0');

%%
figure()
time = mNeondataanalysis{:,1};
mNeonNormIntens = mNeondataanalysis{:,11};
mNeonSeIntens = mNeondataanalysis{:,10};
mNeonSeShoadow = area(time, [(mNeonNormIntens - mNeonSeIntens), (2 * mNeonSeIntens)]);
hold on;
set(mNeonSeShoadow(1),'Visible','off');
set(mNeonSeShoadow(2),'EdgeColor','white','FaceColor',[0.7,0.7,1]);
plot(time, mNeonNormIntens, '-b', 'LineWidth', 2);
xlabel('Time (ms)');
%ylim([0,1]);
ylabel('Normalized Intensity');
%saveas(gca,[pathname '\.fig']);
%saveas(gca,[pathname '\.png']);
