function plotSettings(titleString,labelText,legendLabels,legendLoc)
hxlabel = xlabel(labelText{1});
hylabel = ylabel(labelText{2});
hlegend = legend(legendLabels{:},'Location',legendLoc);
htitle = title(titleString,'FontWeight','normal');
textFontSize = 32;
legendFontSize = 26;
tickFontSize = 26;
set(htitle,'FontSize',textFontSize);
set(hlegend,'FontSize',legendFontSize);
legend boxoff
set(hxlabel,'FontSize',textFontSize);
set(hylabel,'FontSize',textFontSize);
set(gca,'fontsize',tickFontSize);
end