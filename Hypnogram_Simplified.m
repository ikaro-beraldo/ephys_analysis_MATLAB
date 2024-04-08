function Hypnogram_Simplified(sp_plot,All_Sort,b_length, color)%% Hypnogram
%Blocks length in seconds

% Call the subplot
sp_plot;

%Define the block range (default = 1:length(All_Sort))
loop_vec = 1:length(All_Sort);


REM = find(All_Sort(loop_vec) == 1);
nREM = find(All_Sort(loop_vec) == 2);
WAKE = find(All_Sort(loop_vec) == 3);
Excluded = find(All_Sort(loop_vec) == 0);

REM = REM*b_length/60;
nREM = nREM*b_length/60;
WAKE = WAKE*b_length/60;
alt = 1*b_length/60;

for i = 1:length(REM)
    rectangle('Position',[REM(i) 0 alt 1],'FaceColor',color(1,:),'LineStyle','none')
    hold on
end

for i = 1:length(nREM)
    rectangle('Position',[nREM(i) 1 alt 1],'FaceColor',color(2,:),'LineStyle','none')
    hold on
end

for i = 1:length(WAKE)
    rectangle('Position',[WAKE(i) 2 alt 1],'FaceColor',color(3,:),'LineStyle','none')
    hold on
end

for i = 1:length(Excluded)
    rectangle('Position',[Excluded(i) 0 alt 0],'FaceColor',[.5 .5 .5],'LineStyle','none')
    hold on
end
hold off
yticks([0.5 1.5 2.5])
ylim([0 3])
yticklabels({'REM','nREM','WAKE'})
xlim([0 length(loop_vec)*b_length/60])
xlabel('Minutes')

% %Plot grey bars
% %Define the variables for the black lines dividing days
% day2 = ones(1,4) * day_timestamps(2,1) *b_length/3600;
% day3 = ones(1,4) * day_timestamps(3,1) *b_length/3600;
% 
% hold on
% plot(day2,0:3,'Color',[0 0 0],'LineWidth',2);
% hold on
% plot(day3,0:3,'Color',[0 0 0],'LineWidth',2);

%Plot patterns
set(gca,'Tickdir','out')
set(gca,'Linewidth',1.5)
set(gca,'fontname','helvetica')
set(gca,'fontsize',20)

ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
set(gcf,'color','white')
box off

end