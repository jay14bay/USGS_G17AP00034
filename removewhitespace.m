function removewhitespace(plotaxishandle,tightfrac)

ax = plotaxishandle;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + tightfrac.*ti(1);
bottom = outerpos(2) + tightfrac.*ti(2);
ax_width = outerpos(3) - tightfrac.*ti(1) - tightfrac.*ti(3);
ax_height = outerpos(4) - tightfrac.*ti(2) - tightfrac.*ti(4);
ax.Position = [left bottom ax_width ax_height];