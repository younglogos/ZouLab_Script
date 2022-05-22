
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
saveas(gca,[path '\F-V stack.fig']);
saveas(gca,[path '\F-V stack.png']);

xlswrite([pathname '0analysis.xlsx'],{'Vm', 'delta F/F0'},'A1:B1');
xlswrite([pathname '0analysis.xlsx'],{-100, -80, -60, -40, -20, 0, 20, 40, 60, 80}','A2:A11')
xlswrite([pathname '0analysis.xlsx'],norm_inten_plat','B2:B11')
