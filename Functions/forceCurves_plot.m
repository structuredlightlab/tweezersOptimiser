function [fig] = forceCurves_plot(force,pos,dir,xlbl,ttl)

cols = get(groot,'defaultaxescolororder');

Nb = size(force,2); Nd = size(force,1);
fig = figure('position',[150 150 550 300]); hold on
k = 0;
legStr = cell(numel(force),1);
for i = 1:Nb
    for j = 1:Nd
        k = k+1;
        legStr{k} = ['[' num2str(dir(:,j)') ']'];
        if i==1
            plot(pos,force{j,i},'-','color',cols(j,:))
            legStr{k} = [legStr{k} '  Optimised'];
        else
            plot(pos,force{j,i},'--','color',cols(j,:))
            legStr{k} = [legStr{k} '  Gaussian'];
        end
    end
end
xlim([min(pos) max(pos)])
ax = axis;
line(ax(1:2),[0 0],'color',[1 1 1]*.15,'linewidth',0.5)
line([0 0],ax(3:4),'color',[1 1 1]*.15,'linewidth',0.5)
box on
xlabel(xlbl)
ylabel('Force, arb.u.')
title(ttl)
legend('string',legStr,'location','northeastoutside')



end

