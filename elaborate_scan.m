clear; clc;

fontsize = 12;

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


%% SINGLE CHANNELS (THR = 200)

importedData = readmatrix(['input/14092022/charge_scan_14092022_sens0.dat']);

ENC = nan(31, 1);
THR = nan(31, 1);

for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);
    f = figure("Visible", "off");
    box on
    grid on
    hold on
    X = data(:,2)*0.841;
    DATA =  data(:,4)/10;
    plot(X, DATA);
    [f, x, flo, fup] = ecdf(data(:,4)/10, 'Bounds','on');
    hold off
    val_min = max(X(DATA == 0))
    if max(DATA) >= 98 & min(DATA) == 0
        val_max = min(X(DATA >= 98))
        std = val_max - val_min
        THR(ch+1, 1) = val_min + std/2;
        ENC(ch+1, 1) = std/2 * 2.35;
    else
        ENC(ch+1, 1) = nan;
    end
    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
    title("\textbf{Threshold Scan - Ch. " + num2str(ch) + "}");
    exportgraphics(gcf, ['output/single_channels_THR_200/Scan di carica - ch ' num2str(ch) '.pdf'], 'ContentType', 'vector');
end

f = figure("Visible", "on");
hold on
plot([0:7], ENC(1:8))
plot([0:7], THR(1:8))
hold off

data = [[0:7]', round(THR(1:8), 3), round(ENC(1:8), 3)];
data_table = array2table(data, "VariableNames", ["Channel", "Threshold", "ENC"]);
writetable(data_table, "ENC_THR_data_charge_scan_THR_200.dat", "Delimiter", "\t")


%% SINGLE CHANNELS (THR = 214)

importedDataTemp0 = readmatrix(['input/14092022/charge_scan_14092022_sens0.dat']);
importedDataTemp1 = readmatrix(['input/14092022/charge_scan_14092022_sens1.dat']);
importedDataTemp2 = readmatrix(['input/14092022/charge_scan_14092022_sens2.dat']);
importedDataTemp3 = readmatrix(['input/14092022/charge_scan_14092022_sens3.dat']);

importedData = [importedDataTemp0; importedDataTemp1; importedDataTemp2; importedDataTemp3];

importedData = readmatrix(['input/14092022/charge_scan_14092022_sens1.dat']);

ENC = nan(31, 1);
THR = nan(31, 1);

for ch = 8:15
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);

    f = figure("Visible", "on");
    box on
    grid on

    hold on
    X = data(:,2)*0.841;
    DATA =  data(:,4)/10;
    plot(X(1:200), DATA(1:200));

    [f, x, flo, fup] = ecdf(data(:,4)/10, 'Bounds','on');
    hold off

    %val_min = max(X(DATA <= 1))

%     if max(DATA) >= 100 & min(DATA) <= 1
%         val_max = min(X(DATA >= 100))
%         std = val_max - val_min
%         THR(ch+1, 1) = val_max - std/2;
%         ENC(ch+1, 1) = std/2 * 2.35;
%     else
%         val_min = 0;
%         val_max = min(X(DATA >= 100))
%         std = val_max - val_min
%         THR(ch+1, 1) = val_max - std/2;
%         ENC(ch+1, 1) = std/2 * 2.35;
%     end

    X = X(1:200);
    DATA = DATA(1:200);

    val_min = max(X(DATA<=0));
    val_max = min(X(DATA >= 100))
    std = val_max - val_min
    THR(ch+1, 1) = val_min + std/2;
    ENC(ch+1, 1) = std/2 * 2.35;

    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
    title("\textbf{Threshold Scan - Ch. " + num2str(ch) + "}");
    exportgraphics(gcf, ['output/single_channels_THR_214/Scan di carica - ch ' num2str(ch) '.pdf'], 'ContentType', 'vector');
end

f = figure("Visible", "on");
hold on
plot([8:15], ENC(8:15))
plot([8:15], THR(8:15))
hold off

data = [[8:15]', round(THR(16:25), 3), round(ENC(8:15), 3)];
data_table = array2table(data, "VariableNames", ["Channel", "Threshold", "ENC"]);
writetable(data_table, "ENC_THR_data_charge_scan_THR_214_sens1.dat", "Delimiter", "\t")


%% ALL CHANNELS

clear; clc;

f = figure;
colors = distinguishable_colors(32, 'w');

% importedDataTemp0 = readmatrix(['input/14092022/charge_scan_14092022_sens0.dat']);
% importedDataTemp1 = readmatrix(['input/14092022/charge_scan_14092022_sens1.dat']);
% importedDataTemp2 = readmatrix(['input/14092022/charge_scan_14092022_sens2.dat']);
% importedDataTemp3 = readmatrix(['input/14092022/charge_scan_14092022_sens3.dat']);
% 
% importedData = [importedDataTemp0; importedDataTemp1; importedDataTemp2; importedDataTemp3];

importedData = readmatrix(['input/14092022/charge_scan_14092022_sens3.dat']);

hold on
grid on

channels = strings(32, 1);
for ch = 0:31
    channels(ch+1, 1) = strcat("Channel \#", num2str(ch));
end

for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 300,:);
    plot(data(:,2)*0.841,data(:,4)/10, 'Color', [colors(ch+1, 1), colors(ch+1, 2), colors(ch+1, 3)]);
    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
end

box on
grid on
xlim([0 80])
yticks([0:10:100])
%xticks([0:10:120])
legend(channels(1:32), 'Location', 'eastoutside', 'NumColumns', 2)

fontsize = 12;
ax = gca; 
ax.XAxis.FontSize = fontsize; 
ax.YAxis.FontSize = fontsize; 
ax.Legend.FontSize = fontsize; 

%title(['\textbf{Threshold Scan - Detector \#3 (Ch. 24 - 31)}']);
save_image('Threshold Scan - All detectors - TH214.','pdf',f);


%% SAVE DATA 

importedData = readmatrix(['input/14092022/charge_scan_14092022.dat']);

myFitType = fittype(@(a,b,x) 500 + 500*erf((x-a)/(sqrt(2)*b)));

results = [];
for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    myFit = fit(data(:,2)*0.841,data(:,4),myFitType,'Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[20 1]);
    results = [results; ch coeffvalues(myFit)];
end

writematrix(results,'output/data_files/a,b.dat');

f = figure()
hold on
%plot(results(:, 1), results(:, 2))
plot(results(:, 1), results(:, 3).*2.35)

results(:, 3) = results(:, 3) .* 2.35;

data_table = array2table(results, "VariableNames", ["Channel", "Threshold", "ENC"]);
writetable(data_table, "ENC_THR_data_charge_scan_THR_214_sens3.dat", "Delimiter", "\t")
