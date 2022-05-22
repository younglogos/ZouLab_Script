h=open('X:\91 Data and analysis\Xiaotian Bi\20190518 patch\PM-ArcLight(217Q)-FL\cell 3\145757_40x CQ\F-V analysis.fig');
h=get(gca,'children')
data=get(h,{'xdata','ydata'})
close(gcf)

plot(data{1,1}(41:120),data{1,2}(41:120)./data{1,2}(120),'Color',getrgb(532),'LineWidth',2);hold on
box off
axis tight

plot(data{1,1}(73:96), data{1,2}(73:96)./data{1,2}(96)-1, 'g', 'LineWidth',2);
box off;
axis tight;
title('V-F curve of ArcLight(R217Q)');
xlabel('Membrane potential (mV)');
ylabel('¦¤F/F0');
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\V-F_ArcLight.fig']);
saveas(gca,['X:\91 Data and analysis\YJunqi\Screening method\Stimulus pre-experiments\V-F_ArcLight.png']);