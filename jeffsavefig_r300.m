%% save a figure to png

function [fpath] = jeffsavefig_r300(fh,path)
       
set(fh,'PaperPositionMode','auto')
print(fh,path,'-dpng','-r300')
fpath=path;