clear all;
% �������ǲ���ԭʼdata��������Ĥ��λ��������DF/F0, 0��һ������
h=get(gca,'children')
data=get(h,{'xdata','ydata'})
 
name1 = 'SNAP'
WVLT1 = 590;
figure()


F_matrix = zeros(5,60);
%F_matrix = zeros(4,60)
F_matrix(2:4,:) = (reshape(data{1,2}(48:227),60,3))';
%F_matrix(2:3,:) = (reshape(data{1,2}(48:167),60,2))';
F_matrix(1,14:60) = data{1,2}(1:47);
F_matrix(5,1:13) = data{1,2}(228:240);
%F_matrix(4,1:13) = data{1,2}(168:180);

% 48-227��3�������ڣ�1-47�ǡ�ǰ������ڣ�-30mV ��+50 mV�ٵ�-100 mV����228-240�ǡ��������ڣ���-100 mV ��-30 mV��
V_matrix = zeros(5,60);
%V_matrix = zeros(4,60);
V_matrix(2:4,:) = (reshape(data{1,1}(48:227),60,3))';
%V_matrix(2:3,:) = (reshape(data{1,1}(48:167),60,2))';
V_matrix(1,14:60) = data{1,1}(1:47);
V_matrix(5,1:13) = data{1,1}(168:180);
%V_matrix(4,1:13) = data{1,1}(168:180);

V_average = zeros(60,1);
F_average = zeros(60,1);
for y = 1:60;
    temp = F_matrix(:,y);
    F_average(y) = mean(temp(temp~=0));
    temp = V_matrix(:,y);
    V_average(y) = mean(temp(temp~=0));
end
plot(V_average,F_average)

%������Ϊ�˹���һ��ƽ���ģ����Ƿ����-100 mV ��-50 mV�����ǲ����� ��ÿ������4��ɨ����ƽ���������������ʹ�ã���ʱһ��mVֵ��Ӧ����Fֵ������͸���rampɨ���л�ã�
average_matrix = [V_average F_average];
average_matrix = sortrows(average_matrix,1) %����ɨ�跽�򣬽�һ�������ڵģ�V,F����Ĥ��λ�Ӹ���������

linearfit = fit(average_matrix(:,1),average_matrix(:,2)-100,'poly1')
[a,b] = corr(average_matrix(:,1),average_matrix(:,2)-100) %������ϲ���Pearson���ϵ��
% [a,b] = corr(data{1,1}',data{1,2}'-100) %������ϲ���Pearson���ϵ��
figure();
plot(data{1,1},data{1,2}-100,'.','Color',mean(pseudocolor(WVLT1)',2),'MarkerSize',12);hold on;
xlabel('Membrane potential (mV)')
ylabel('Normalized \DeltaF/F_0 (%) ')
box off
axis tight
plot([V_average],[linearfit.p1*V_average+linearfit.p2],'LineWidth',3,'Color',[0,0,0]);
title({[name1 ', Pearson correlation coefficient = ' num2str(a)];
      [ '\DeltaF/F_0 = ' num2str(linearfit.p1*100,3) '%' ' ��V_m (mV) + ' num2str(linearfit.p2,3) '%']});
  
  %% For D81Y
  linearfit = fit(average_matrix(1:40,1),average_matrix(1:40,2)-100,'poly1')
[a,b] = corr(average_matrix(1:40,1),average_matrix(1:40,2)-100) %������ϲ���Pearson���ϵ��

figure();
plot(data{1,1},data{1,2}-100,'.','Color',mean(pseudocolor(WVLT1)',2),'MarkerSize',12);hold on;
xlabel('Membrane potential (mV)')
ylabel('Normalized \DeltaF/F_0 (%) ')
box off
axis tight
plot([V_average(1:40)],[linearfit.p1*V_average(1:40)+linearfit.p2],'LineWidth',3,'Color',[0,0,0]);
title({[name1 ', Pearson correlation coefficient = ' num2str(a)];
      [ '\DeltaF/F_0 = ' num2str(linearfit.p1*100,3) '%' ' ��V_m (mV) + ' num2str(linearfit.p2,3) '%']});


