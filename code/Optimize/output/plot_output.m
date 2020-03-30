% Figure 5a
% sd = csvread('sd_output.csv', 1, 0);
% sd = sd(:,2:3);
% cg = csvread('cg_output_f_const.csv', 1, 0);
% cg = cg(:,2:3);
% nt = csvread('nt_output_f_const.csv', 1, 0);
% nt = nt(:,2:3);

% figure();
% plot(sd(1:100,1),'r', 'LineWidth', 2);
% hold on;
% plot(cg(:,1),'b', 'LineWidth', 2);
% hold on;
% plot(nt(:,1),'g', 'LineWidth', 2);
% l = legend('SDM','CG','NM');
% xlabel('Iterations','FontSize',16)
% ylabel('norm(F(x^k) - F(x^*))','FontSize',16)

% l.FontSize = 16;
% set(gca, 'YScale', 'log')
% set(gcf,'color','white');
% hLegend = findobj(gcf, 'Type', 'Legend');
% set(hLegend,'Box','off','Fontsize',16,'Location','southeast');
% haxes=findobj(gcf,'Type','axe');
% set(haxes,'Fontsize',16,'Box','on');
% hline = findobj('Type', 'line');
% savefig('compare_three_algs1.fig')
% savefig('compare_three_algs1.eps');

% Figure 5b
cg = csvread('cg_output_f_var.csv', 1, 0);
cg = cg(:,2:3);
nt = csvread('nt_output_f_var.csv', 1, 0);
nt = nt(:,2:3);

figure();
plot(cg(:,1),'b','LineWidth', 2);
hold on;
plot(nt(:,1),'g','LineWidth', 2);
l = legend('CG','NM');
xlabel('Iterations','FontSize',16)
ylabel('norm(F(x^k) - F(x^*))','FontSize',16)

l.FontSize = 16;
set(gca, 'YScale', 'log')
set(gcf,'color','white');
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Box','off','Fontsize',16,'Location','southeast');
haxes=findobj(gcf,'Type','axe');
set(haxes,'Fontsize',16,'Box','on');
hline = findobj('Type', 'line');
savefig('cg_f_var.fig')

% Other Figures
% figure();
% plot(sd(1:100,2),'r', 'LineWidth', 2);
% hold on;
% plot(cg(:,2),'b', 'LineWidth', 2);
% hold on;
% plot(nt(:,2),'g', 'LineWidth', 2);
% savefig('compare_three_algs2.fig')
% xlabel('Iterations','FontSize',12)
% ylabel('Gradient Norm', 'FontSize',12)
% savefig('compare_three_algs2.fig')
