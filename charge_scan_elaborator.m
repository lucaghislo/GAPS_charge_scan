%% ALL CHANNELS (sei nella sezione giusta)

clear; clc;
fontsize = 12; %#ok<NASGU> 

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

clear; clc;

% Choose charge scan to compute
filename = 'IT_L4R0M5_Gigi_charge_scan_THR_205_FTH_MX';
importedData = readmatrix(['input/SSL_Berkeley/FTH/', filename, '.dat']);

f = figure;
colors = distinguishable_colors(32, 'w');
hold on
grid on

channels = strings(32, 1);
for ch = 0:31
    channels(ch+1, 1) = strcat("Ch. ", num2str(ch));
end

data = importedData(importedData(:,5)==0,1:5);
data = data(data(:,2) < 300,:);
data_out = data(:,4)/10;

data_table_out = nan(length(data), 33);
data_table_out(:, 1) = data(:,2)*0.841;

for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 300,:);
    data_out = data(:,4)/10;
    plot(data(:,2)*0.841, data_out, 'Color', [colors(ch+1, 1), colors(ch+1, 2), colors(ch+1, 3)]);
    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
    data_table_out(:, ch+2) = data_out;
end

box on
grid on
yticks([0:10:100])
ylim([0, 100])
legend(channels(1:32), 'Location', 'eastoutside', 'NumColumns', 2);
title_string = strrep(filename, '_', '\_');
title("\textbf{" + title_string + "}");

fontsize = 12;
ax = gca; 
ax.XAxis.FontSize = fontsize; 
ax.YAxis.FontSize = fontsize; 
ax.Legend.FontSize = fontsize; 

f.Position = [200 160 900  550];
exportgraphics(gcf, "output\SSL_Berkeley\FTH\" + string(filename) + ".pdf", 'ContentType','vector');
writematrix(data_table_out, "output\SSL_Berkeley\FTH\data\" + string(filename) + ".dat", "Delimiter", "\t");

disp("Exported: " + string(filename))
