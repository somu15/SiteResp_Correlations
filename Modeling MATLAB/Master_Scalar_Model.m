clr

% cd('C:\Users\lakshd5\Dropbox\Site Response IM Selection\Analysis 2 032619 Regression');

AF = importdata('E:\SiteResp_GMSelect\Python_SourceCode\AmpFac_0.25.txt');
SaR = importdata('E:\SiteResp_GMSelect\Python_SourceCode\SaRock_0.25.txt');

Periods = [0.01, 0.025, 0.05,0.075, 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6,0.7,0.85,1,1.5,2,2.5,3,3.75,4,5];

OPTIONS = optimset('fminsearch');
OPTIONS = optimset(OPTIONS,'MaxFunEvals',100000);
OPTIONS = optimset(OPTIONS,'MaxIter',10000);
OPTIONS = optimset(OPTIONS,'TolFun',0.0001);
OPTIONS = optimset(OPTIONS,'TolX',0.0001);

for ii = 1:length(Periods)
    
    [thetaC,~]= fminsearchbnd('Likelihood',[0.2 0.01 0.4 0.01],[0 -2 0 0],[2 0 1.5 0.1],OPTIONS,...
            SaR(ii,:)', (AF(ii,:)'));
    theta(ii,:) = thetaC;
    progressbar(ii/length(Periods))
    
end

for ii = 1:length(Periods)
ind = ii;
Sar_req = 0.0001:0.0001:1;
pred_req = exp(Pred(theta(ii,:), Sar_req));
fig = figure;
% loglog(SaR(ii,:),AF(ii,:),'o')
% hold on
plot(Sar_req,pred_req,'linewidth',3);
% ylabel(strcat('AF(',num2str(Periods(ii)),'s)'))
% set(gca,'fontsize',16)
% xlabel(strcat('SA(',num2str(Periods(ii)),'s)'))
% title(strcat('sigma=',num2str(theta(ii,3))))

% plot(t, S)
opt = [];
opt.ShowBox = 'on';
opt.LineWidth = [3];
opt.YLabel = strcat('AF(',num2str(Periods(ii)),'s)');
opt.XLabel = strcat('SA(',num2str(Periods(ii)),'s)');
opt.Title = strcat('$\sigma =~$',num2str(theta(ii,3)));opt.YLim = [1e-9 1e-2];
opt.XLim = [1e-4 1e0];
opt.YLim = [0 6];
opt.YTick = [0 1 2 3 4 5 6];
opt.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
opt.Colors = [0 0 0];
opt.BoxDim = [3.5,3.5];
opt.FontName = 'Times';
opt.XMinorTick = 'off';
opt.YMinorTick = 'off';
opt.FontSize = 14;
opt.AxisLineWidth = 0.7;
opt.LineStyle = {'-'};
opt.TickDir = 'in';
opt.XScale = 'log';
% opt.YScale = 'log';
setPlotProp(opt);
%axis tight
grid on
hold on
sz = 15;%41/256    163/256    41/256
%scatter(t(ind), Burton1(ind)*100)
scatter(SaR(ii,:),AF(ii,:),sz,'MarkerEdgeColor',[255/255    51/255    0],...
              'MarkerFaceColor',[255/255    51/255    0],...
              'LineWidth',2);
          semilogx(Sar_req,pred_req,'k','linewidth',3);

saveas(fig,strcat('E:\SiteResp_GMSelect\Data\Dummy motions\Scalar\Period_',num2str(Periods(ii)),'.eps'))
end
close all