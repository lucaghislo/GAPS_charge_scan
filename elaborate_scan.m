clear; clc;

fontsize = 12;

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

folder_name = '';
importedData = readmatrix([folder_name 'scan_carica_TH_200.dat']);

for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);
    %f = figure;
    %grid on
    %plot(data(:,2)*0.841,data(:,4)/10);
    %xlabel('Incoming Energy [keV]');
    %ylabel('Hit [%]');
    
    %title(['Threshold Scan - ch ' num2str(ch)]);
    %save_image(['Scan di carica - ch ' num2str(ch) '.'],'pdf',f);
    %close
end

f = figure;
hold on
grid on
for ch = 24:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);
    plot(data(:,2)*0.841,data(:,4)/10);
    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
end

title(['\textbf{Threshold Scan - Detector \#3 (Ch. 24 - 31)}']);
save_image('Threshold Scan - Detector 3 - TH200.','pdf',f);


%% 
myFitType = fittype(@(a,b,x) 500 + 500*erf((x-a)/(sqrt(2)*b)));

results = [];
for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    myFit = fit(data(:,2)*0.841,data(:,4),myFitType,'Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[20 1]);
    results = [results; ch coeffvalues(myFit)];
end

writematrix(results,'/data_files/a,b.dat');