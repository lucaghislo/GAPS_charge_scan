function save_image(folder_name, extension, figure)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    image_to_save = [folder_name extension];
    if(strcmp(extension,'pdf'))
          set(figure, 'PaperUnits','centimeters');
          set(figure, 'Units','centimeters');
          pos=get(figure,'Position');
          set(figure, 'PaperSize', [pos(3) pos(4)]);
          set(figure, 'PaperPositionMode', 'manual');
          set(figure, 'PaperPosition',[0 0 pos(3) pos(4)]);

          print(image_to_save,'-dpdf')
          
    elseif(strcmp(extension,'svg'))
          print(image_to_save,'-dsvg')       
    end
end

