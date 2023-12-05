
%%Set Parameters to analyse behavioural data
dataDirectory = '\\mhsdata.anu.edu.au\mydocs$\u6764904\My Documents\PhD\Behavioiural_Studies\M1AChRAnalysis\';
% Master List for day files
% % % Mouse 1
whichDays{1}{1} = {''}; %ACSF, add session dates here
whichDays{1}{2} = {''}; %AG
whichDays{1}{3} = {''}; %AN
% %  Mouse 2
whichDays{2}{1} = {''}; %ACSF 
whichDays{2}{2} = {''}; %AG
whichDays{2}{3} = {''}; %AN
% % Mouse 3
whichDays{3}{1} = {''}; %ACSF 
whichDays{3}{2} = {''}; %AG
whichDays{3}{3} = {''}; %AN
% % % Mouse 4
whichDays{4}{1} = {''}; %ACSF 
whichDays{4}{2} = {''}; %AG
whichDays{4}{3} = {''}; %AN
% % % %  Mouse 5
whichDays{5}{1} = {''}; %ACSF 
whichDays{5}{2} = {''}; %AG
whichDays{5}{3} = {''}; %AN
% % % %  Mouse 6
whichDays{6}{1} = {''}; %ACSF 
whichDays{6}{2} = {''}; %AG
whichDays{6}{3} = {''}; %AN

%% Basic Parameters
preEventWindowTemp = 1000; 
preEventWindow = 1000;  
postEventWindow = 1000;  
psthBin = 10;
firstLickBin = 0:100:1000;
postStimWindow = 400;
colormap{1} = [0 0 0; 0.4 0.4 0.4; 0.55 0.55 0.55; 0.65 0.65 0.65; 0.8 0.8 0.8];
colormap{2} = [0 0.4 0; 0 0.55 0; 0 0.7 0; 0 0.85 0; 0 1 0];
colormap{3} = [0.4 0 0; 0.55 0 0; 0.7 0 0; 0.85 0 0; 1 0 0];
stimList = [0 15 30 60 120];

%% Load Data
for mouseID = 1:length(whichDays)
    for condition = 1:length(whichDays{mouseID})
        for dayID = 1:length(whichDays{mouseID}{condition})
            if condition == 1
                List = dir([dataDirectory 'M1_ACSF-Mouse-' num2str(mouseID) '--2021-' whichDays{mouseID}{condition}{dayID} '*.mat']); 
            elseif condition == 2
                List = dir([dataDirectory 'M1_BQCA-Mouse-' num2str(mouseID) '--2021-' whichDays{mouseID}{condition}{dayID} '*.mat']);
            else
                List = dir([dataDirectory 'M1_TD-Mouse-' num2str(mouseID) '--2021-' whichDays{mouseID}{condition}{dayID} '*.mat']); 
            end
            
            if isempty(List)
                disp('Fail to Load File')
                disp([List.name]);
                continue
            else
                
                load([dataDirectory List.name])
                pseudoBlockSize=myData.psuedoBlockSize;
                numTrials=myData.numOfTrial;
                allData=myData.allData;
                
                for b = 1:size(allData,2)/5
                    
                    block = allData(:,(b-1)*5+1:b*5); %all size rows for that particular block of 5 stim trials
                    
                    noStim = find(block(1,:)==0);
                    tmp = 1:5;
                    tmp = tmp(tmp~=noStim);
                    
                    performance(mouseID,dayID,condition,b) = sum(block(6,tmp)>0)/4;
                    hits(mouseID,dayID,condition,b) = sum(block(6,tmp)>0)/4;
                    performance2(mouseID,dayID,condition,b) = sum(block(6,tmp)>0)/4;
                    if block(6,noStim)>0
                        performance2(mouseID,dayID,condition,b)=NaN;
                    end
                    
                    fa(mouseID,dayID,condition,b)=sum(block(6,noStim)>0)/1;
                    
                    dp(mouseID,dayID,condition,b) = norminv(performance(mouseID,dayID,condition,b)-0.001)-norminv(fa(mouseID,dayID,condition,b)+0.001);

                end
                
           
                newAllData=reshape(allData(:,1:numTrials),size(allData,1),pseudoBlockSize,numTrials/pseudoBlockSize);
                %(data type) x (block-trial) x (block)
                blockIndices=[];%this is a list of trial blocks with licks for 1 stim
                for b=1:size(newAllData,3)
                 
                    idx=newAllData(1,:,b)==0; %finds the index of 0 stimulus trial within each block

                    if newAllData(6,idx,b)==0 %check if there is no lick assoc with the 0 stim
                        blockIndices=[blockIndices b];%append indices of good blocks to list
                    end
                     
                end
                newAllData=newAllData(:,:,blockIndices);%subsample only good blocks
                newAllData=reshape(newAllData,size(newAllData,1),size(newAllData,2)*size(newAllData,3));

                behaviourMatrix{condition}{mouseID}{dayID} = newAllData';
                vibrationAmplitude{condition}{mouseID}{dayID} = myData.vibrationAmp;
                lickOnsets{condition}{mouseID}{dayID} = myData.behaviourTime(myData.behaviour==1)*1000;
                
                for intensity = 1:5
                    % Generate PSTH
                    stimulusOnset{condition}{intensity}{mouseID}{dayID} = behaviourMatrix{condition}{mouseID}{dayID}((behaviourMatrix{condition}{mouseID}{dayID}(:,1)==vibrationAmplitude{condition}{mouseID}{dayID}(intensity)),5)*1000;
                    for i = 1:length(stimulusOnset{condition}{intensity}{mouseID}{dayID})
                        lickRaster{condition}{intensity}{mouseID}{dayID}{i} = lickOnsets{condition}{mouseID}{dayID}((lickOnsets{condition}{mouseID}{dayID}>(stimulusOnset{condition}{intensity}{mouseID}{dayID}(i)-preEventWindowTemp) &...
                            lickOnsets{condition}{mouseID}{dayID}<(stimulusOnset{condition}{intensity}{mouseID}{dayID}(i)+postEventWindow))) - stimulusOnset{condition}{intensity}{mouseID}{dayID}(i);
                    end
                    psth{condition}{intensity}{mouseID}{dayID} = zeros(length(stimulusOnset{condition}{intensity}{mouseID}{dayID}),(preEventWindowTemp + postEventWindow));
                    for trial = 1:length(stimulusOnset{condition}{intensity}{mouseID}{dayID})
                        psth{condition}{intensity}{mouseID}{dayID}(trial,ceil(lickRaster{condition}{intensity}{mouseID}{dayID}{trial}+preEventWindowTemp))=1;
                    end
                    psth{condition}{intensity}{mouseID}{dayID} = psth{condition}{intensity}{mouseID}{dayID}(:,[1:preEventWindowTemp-50 preEventWindowTemp:end]);
                    psth{condition}{intensity}{mouseID}{dayID} = psth{condition}{intensity}{mouseID}{dayID}(:,(size(psth{condition}{intensity}{mouseID}{dayID},2)-(preEventWindow+postEventWindow)+1):end);
                    PSTH{condition}{intensity}{mouseID}{dayID} = sum(reshape(sum(psth{condition}{intensity}{mouseID}{dayID}),psthBin,[]))/size(psth{condition}{intensity}{mouseID}{dayID},1)/psthBin.*1000;
                    
                    % Get Lick Rate + Normalized
                    postStimWindowRaster{condition}{intensity}{mouseID}{dayID} = psth{condition}{intensity}{mouseID}{dayID}(:,preEventWindow:preEventWindow+postStimWindow-1);
                    postStimWindowRate{condition}{intensity}{mouseID}(dayID) = mean(sum(reshape(sum(postStimWindowRaster{condition}{intensity}{mouseID}{dayID}),psthBin,[]))/size(postStimWindowRaster{condition}{intensity}{mouseID}{dayID},1)/psthBin.*1000);
                    preStimWindowRaster{condition}{intensity}{mouseID}{dayID} = psth{condition}{intensity}{mouseID}{dayID}(:,preEventWindow-postStimWindow:preEventWindow-1);
                    preStimWindowRate{condition}{intensity}{mouseID}(dayID) = mean(sum(reshape(sum(preStimWindowRaster{condition}{intensity}{mouseID}{dayID}),psthBin,[]))/size(preStimWindowRaster{condition}{intensity}{mouseID}{dayID},1)/psthBin.*1000);
                    lickRateNormalized{condition}{intensity}{mouseID}(dayID) = postStimWindowRate{condition}{intensity}{mouseID}(dayID)- preStimWindowRate{condition}{intensity}{mouseID}(dayID);
                    
                    % Get First Lick
                    firstLickTime{condition}{intensity}{mouseID}{dayID} = cell2mat(cellfun(@(x) x(find(x>0,1)), lickRaster{condition}{intensity}{mouseID}{dayID}, 'UniformOutput', false));
                    [firstLickHist{condition}{intensity}{mouseID}{dayID}, ~] = hist(firstLickTime{condition}{intensity}{mouseID}{dayID}, firstLickBin);
                    
                end
            end
        end
    end
end

%% Plot Average PSTH per mouse
plot_flag =input('Plot each individual mice? (0=No; 1= Yes)');
slidingWindow = 100;
for mouseID = 1:length(whichDays)
    figure
    for condition = 1:length(whichDays{mouseID})
        for intensity = 1:length(PSTH{condition})
            meanPSTH{mouseID}{condition}{intensity} = mean(cell2mat(PSTH{condition}{intensity}{mouseID}'));
            for i = 1:length(meanPSTH{mouseID}{condition}{intensity})-slidingWindow
                slidingPSTH{condition}{intensity}{mouseID}(i) = mean(meanPSTH{mouseID}{condition}{intensity}(i:i+slidingWindow));
            end
            x = linspace(-preEventWindow/1000,postEventWindow/1000,length(meanPSTH{mouseID}{condition}{intensity}));
            slidingPSTH{condition}{intensity}{mouseID} = [NaN(1,100) slidingPSTH{condition}{intensity}{mouseID}];
            
            subplot(1,3,condition)
            hold on
            plot(x, slidingPSTH{condition}{intensity}{mouseID},'color', colormap{condition}(intensity,:));
            ylim([0 4])
            xlim([-1 2])
            
            
        end
        subplot(1,3,condition)
        plot([0 0],ylim,'k--')
        set(gca,'box','off','TickDir', 'out')
        if condition == 1
            title('\color{black}aCSF')
            ylabel('Lick /s')
        elseif condition == 2
            title('\color{green}BQCA')
        else
            title('\color{red}TD')
        end
        xlabel('Time since stim onset (s)')
        axis square
    end
end
if plot_flag == 0
    close all
end

%% Performance across all session for all animals 
figure;
for m = 1:6
    subplot(2,3,m);
    hold on
    %     plot(squeeze(performance(m,:,1,:)),'k')
    plot(squeeze(mean(performance(m,:,1,1:40))),'k')
    plot(squeeze(mean(performance(m,:,2,1:40))),'g')
    plot(squeeze(mean(performance(m,:,3,1:40))),'r')
    
    ylim([0 1])
    
end


%% Plot Hit rates and FA with trial block

figure;
for m = 2
%     subplot(2,3,m);
    hold on
    grey = [.7 .7 .7];
    plot(squeeze(mean(hits(m,:,1,1:40))),'Color', grey)
    plot(squeeze(mean(hits(m,:,2,1:40))),'g')
    plot(squeeze(mean(hits(m,:,3,1:40))),'r')
    xlabel('Trial blocks')
    ylabel('Hit Rate')
    ylim([0 1])
    axis square
end

figure;
for m = 2
%     subplot(2,3,m);
    hold on
    grey = [.7 .7 .7];
    plot(squeeze(mean(fa(m,:,1,1:40))),'Color', grey)
    plot(squeeze(mean(fa(m,:,2,1:40))),'g')
    plot(squeeze(mean(fa(m,:,3,1:40))),'r')
    xlabel('Trial blocks')
    ylabel('FA')
    ylim([0 1])
    axis square
end


%% mean false alarms and hits across mice
figure;
plot(squeeze(mean(mean(hits(2:3,:,1,1:40)))),'Color', grey)
hold on;
plot(squeeze(mean(mean(hits(2:3,:,2,1:40)))),'g')
hold on;
plot(squeeze(mean(mean(hits(2:3,:,3,1:40)))),'r')
hold on;
xlabel('Trial blocks')
ylabel('Hit Rate')
ylim([0 1])

figure;
plot(squeeze(mean(mean(fa(5:6,:,1,1:40)))),'Color', grey)
hold on;
plot(squeeze(mean(mean(fa(5:6,:,2,1:40)))),'g')
hold on;
plot(squeeze(mean(mean(fa(5:6,:,3,1:40)))),'r')
hold on;
xlabel('Trial blocks')
ylabel('FA')
ylim([0 1])

figure;
plot(squeeze(mean(dp(4,:,1,1:48))),'Color', grey)
hold on;
plot(squeeze(mean(dp(4,:,2,1:48))),'g')
hold on;
plot(squeeze(mean(dp(4,:,3,1:48))),'r')
hold on;
xlabel('Trials')
ylabel('d-prime')
% ylim([0 1])
b=zeros(size(dp,1),size(dp,3));

for c=1
    for m = 4
        a=squeeze(mean(dp(m,:,c,1:40)));
        b(m,c) = mean(a,"omitnan")
        
    end
end

a=squeeze(dp(4,:,2,1:48));
c=mean(a,"omitnan")
c=mean(c,"omitnan")

%% Histogram of performance for all animals
figure;
xbins = 0:.1:1;
for m = 1:6
    subplot(2,3,m)
    hold on
    
    [y x] = hist(squeeze(mean(performance(m,:,1,1:48))),xbins);
    plot(x,y,'k')
    tmp = squeeze(mean(mean(performance(m,:,1,1:48))));
    plot([tmp tmp],ylim,'k--')
    
    
    [y x] = hist(squeeze(mean(performance(m,:,2,1:48))),xbins);
    plot(x,y,'g')
    tmp = squeeze(mean(mean(performance(m,:,2,1:48))));
    plot([tmp tmp],ylim,'g--')
    
    
    [y x] = hist(squeeze(mean(performance(m,:,3,1:48))),xbins);
    plot(x,y,'r')
    tmp = squeeze(mean(mean(performance(m,:,3,1:48))));
    plot([tmp tmp],ylim,'r--')
    
end

%% Performance for 3 conditions for all mice

endBlock = 20; %performance till how many trials
figure;
ylim ([0 1]);
hold on
errorbar(1:10:60,squeeze(mean(mean(performance(:,:,1,1:endBlock),2),4)),std(squeeze(mean(performance(:,:,1,1:endBlock),4))')/sqrt(5),'.k')
bar(1:10:60,squeeze(mean(mean(performance(:,:,1,1:endBlock),2),4)),0.2,'g')

errorbar(3:10:60,squeeze(mean(mean(performance(:,:,2,1:endBlock),2),4)),std(squeeze(mean(performance(:,:,2,1:endBlock),4))')/sqrt(5),'.k')
bar(3:10:60,squeeze(mean(mean(performance(:,:,2,1:endBlock),2),4)),0.2,'g')

errorbar(5:10:60,squeeze(mean(mean(performance(:,:,3,1:endBlock),2),4)),std(squeeze(mean(performance(:,:,3,1:endBlock),4))')/sqrt(5),'.k')
bar(5:10:60,squeeze(mean(mean(performance(:,:,3,1:endBlock),2),4)),0.2,'g')
%% FA for all 6 animals 

endBlock = 20; %performance till how many trials
figure;
ylim ([0 1]);
hold on
errorbar(1:10:60,squeeze(mean(mean(fa(:,:,1,1:endBlock),2),4)),std(squeeze(mean(fa(:,:,1,1:endBlock),4))')/sqrt(5),'.k')
bar(1:10:60,squeeze(mean(mean(fa(:,:,1,1:endBlock),2),4)),0.2,'g')

errorbar(3:10:60,squeeze(mean(mean(fa(:,:,2,1:endBlock),2),4)),std(squeeze(mean(fa(:,:,2,1:endBlock),4))')/sqrt(5),'.k')
bar(3:10:60,squeeze(mean(mean(fa(:,:,2,1:endBlock),2),4)),0.2,'g')

errorbar(5:10:60,squeeze(mean(mean(fa(:,:,3,1:endBlock),2),4)),std(squeeze(mean(fa(:,:,3,1:endBlock),4))')/sqrt(5),'.k')
bar(5:10:60,squeeze(mean(mean(fa(:,:,3,1:endBlock),2),4)),0.2,'g')

% Preallocate a matrix to store the means for each condition and mouse
numMice = 6;
numConditions = 3;
endBlock = 10; % performance till how many trials

meansMatrix = zeros(numMice, numConditions);

for mouse = 1:numMice
    for condition = 1:numConditions
        % Extract the mean false alarms and store in the matrix
        meansMatrix(mouse, condition) = mean(mean(fa(mouse,:,condition,1:endBlock), 4), 2);
    end
end

%% scatter plot based on previour figure
num_mice =6;       % Number of mice
num_days = 5;       % Number of days
num_conditions = 3; % Number of conditions
num_trials = 20;    % Number of trials
performance = fa;
% figure;
hold on;

for mouse = num_mice
    % Initialize arrays to store means and individual data points for each day
    means = zeros(num_days, num_conditions);
    err = zeros(num_days, num_conditions);
    scatter_data = cell(num_days, num_conditions);
    
    for day = 1:num_days
        for condition = 1:num_conditions
            % Extract data for the current mouse, day, and condition
            data = squeeze(performance(mouse, day, condition, 1:endBlock));
            err(day,condition) = std(squeeze(mean(performance(mouse,day,condition,1:endBlock),2))')/sqrt(5)
            %             err = std(mean(performance(mouse,day,condition,:)))/sqrt(5);
            % Calculate the mean for this day and condition
            mean_data = mean(data);
            mean_err = mean (err);
            % Store the mean
            means(day, condition) = mean_data;
            
            % Store individual data points for this day and condition
            scatter_data{day, condition} = data;
        end
    end
    
    
    % Plot grouped bars with error bars
    %     bar(x, mean(means,1), 'grouped');
    %     hold on;
    %     ylim([0 1]);
    
    % Plot individual data points as scatter points
    scatter_x = repmat([51 53 55], size(means, 1), 1); %chnage values as per the previous bar plot
    scatter_x = scatter_x(:) %+ randn(size(scatter_x(:))) * 0.05; % Add jitter for better visualization
    scatter_y = means(:) %+ randn(size(means(:))) * 0.05; % Add jitter for better visualization
    scatter(scatter_x, scatter_y, 20, 'k' , 'filled');

    hold on;
end


%% Plot scatter with errorbars  - DO NOT EDIT
figure;
hold on;

% Define x-values for the bar groups
x = 1:size(means, 2);

for mouse = 1:num_mice
    % Initialize arrays to store means and individual data points for each day
    means = zeros(num_days, num_conditions);
    scatter_data = cell(num_days, num_conditions);
    
    for day = 1:num_days
        for condition = 1:num_conditions
            % Extract data for the current mouse, day, and condition
            data = squeeze(performance(mouse, day, condition, :));
            
            % Calculate the mean for this day and condition
            mean_data = mean(data);
            
            % Store the mean
            means(day, condition) = mean_data;
            
            % Store individual data points for this day and condition
            scatter_data{day, condition} = data;
        end
    end

    
    % Plot grouped bars with error bars
    bar(x, mean(means, 1), 'grouped');
    hold on;
    
    % Plot individual data points as scatter points with error bars
    for condition = 1:num_conditions
        x_scatter = x(condition) + randn(num_days, 1) * 0.05; % Add jitter for better visualization
        y_scatter = means(:, condition);
        errorbar(x_scatter, y_scatter, err(:, condition), 'k.', 'LineWidth', 0.5);
    end
    
    hold off;
    
    % Customize the plot labels and legend
    xlabel('Condition');
    ylabel('Detection Rate');
    xticks(x);
    xticklabels(1:num_conditions);
    legend('Session');
    
    % Adjust the layout
    grid on;
    box on;
end

%% Plot Average PSTH across all Mice
figure
for condition = 1:3
    subplot(3,1,condition)
    hold on
    for intensity = 1:length(slidingPSTH{condition})
        sliding_PSTH_mean = nanmean(cell2mat(slidingPSTH{condition}{intensity}'));
        sliding_PSTH_SEM = nanstd(cell2mat(slidingPSTH{condition}{intensity}'))/sqrt(length(length(whichDays)));       
        shadedErrorBar(x,nanmean((cell2mat(slidingPSTH{condition}{intensity}')),1),sliding_PSTH_SEM,{'color', colormap{condition}(intensity,:)},1);
    end
    ylim([0 6])
    xlim([-1 4])
    plot([0 0],ylim,'k--')
    set(gca,'box','off','TickDir', 'out')
    if condition == 1
        title('\color{black}ACSF')
        ylabel('Lick /s')
    elseif condition == 2
        title('\color{green}BQCA')
    else
        title('\color{red}TD')
    end
    xlabel('Time since stim onset (s)')
%     axis square
end

%% First Lick Reaction Time across all Rats

for condition = 1:3
    for intensity = 2:5
        for mouseID = 1:length(firstLickTime{condition}{intensity})
            firstLickCombinedPerRat{condition}{intensity}{mouseID} = cell2mat(firstLickTime{condition}{intensity}{mouseID});
        end
        firstLickCombined{condition}{intensity} = cell2mat(firstLickCombinedPerRat{condition}{intensity});
        %         firstLick {mouseID} = cell2mat(firstLickCombinedPerRat{condition}{intensity})
        
    end
end
figure
hold on
errorbar(stimList,cellfun(@(x) nanmean(x), firstLickCombined{1}),cellfun(@(x) nanstd(x)/sqrt(length(x)), firstLickCombined{1}),'k.-','MarkerSize',12)
errorbar(stimList,cellfun(@(x) nanmean(x), firstLickCombined{2}),cellfun(@(x) nanstd(x)/sqrt(length(x)), firstLickCombined{2}),'g.-','MarkerSize',12)
errorbar(stimList,cellfun(@(x) nanmean(x), firstLickCombined{3}),cellfun(@(x) nanstd(x)/sqrt(length(x)), firstLickCombined{3}),'r.-','MarkerSize',12)
set(gca,'box','off','TickDir', 'out')
ylabel('First lick time (ms)')
xlabel('Stimulus intensity (um)')
% xlim([-5 85])

%%
% Initialize a cell array to store individual first lick times for each mouse
firstLickCombined = cell(1, length(stimList));

for condition = 3
    for intensity = 5
        for mouseID = 1:length(firstLickTime{condition}{intensity})
            firstLickCombinedPerRat{condition}{intensity}{mouseID} = cell2mat(firstLickTime{condition}{intensity}{mouseID});
            % Store individual data points for each mouse in the combined array
            firstLickCombined{mouseID} = firstLickCombinedPerRat{condition}{intensity}{mouseID};
        end
    end
end

% Plot individual data points for all 6 mice for each condition
figure
hold on
for mouseID = 1:length(firstLickCombined)
    scatter(length(firstLickCombined), firstLickCombined{mouseID}, 'filled')
end

% Plot mean and error bars for each condition
errorbar(stimList,cellfun(@(x) nanmean(x), firstLickCombined),cellfun(@(x) nanstd(x)/sqrt(length(x)), firstLickCombined),'k.-','MarkerSize',12)
set(gca,'box','off','TickDir', 'out')
ylabel('First lick time (ms)')
xlabel('Stimulus intensity (um)')



%% Quanitfy increases in Lick rate

for condition = 1:3
    for mouseID = 1:length(whichDays)
        faLickRateTemp{condition}{mouseID} = postStimWindowRate{condition}{5}{mouseID};
    end
    faLickRate{condition} = cell2mat(faLickRateTemp{condition});
end
figure
hold on
errorbar(1:3,cellfun(@(x) nanmean(x), faLickRate),cellfun(@(x) nanstd(x)/sqrt(length(x)), faLickRate),'k.')
bar(1, nanmean(faLickRate{1}),'k')
bar(2, nanmean(faLickRate{2}),'g')
bar(3, nanmean(faLickRate{3}),'r')
set(gca,'box','off','TickDir', 'out')
ylabel('Average Lick /s')
set(gca,'XTickLabel',{' ' 'ACSF' 'BQCA' 'TD' ' '})
set(gcf,'Position',[744 625 318 425])
% ylim ([0 2.5])

%% Plot Psychometric Function
x = stimList;
figure
for condition = 1:3
    for intensity = 1:5
        lickRateNormalizedAllRats{condition}{intensity} = cell2mat(lickRateNormalized{condition}{intensity});
    end
    hold on
    y = cellfun(@(x) nanmean(x), lickRateNormalizedAllRats{condition});
    if condition == 1
        errorbar(x,y,cellfun(@(x) std(x)/sqrt(length(x)), lickRateNormalizedAllRats{condition}),'ko-','MarkerSize', 5,'MarkerFaceColor', 'k')
    elseif condition == 2
        errorbar(x,y,cellfun(@(x) std(x)/sqrt(length(x)), lickRateNormalizedAllRats{condition}),'go-','MarkerSize', 5,'MarkerFaceColor', 'g')
    elseif condition == 3
        errorbar(x,y,cellfun(@(x) std(x)/sqrt(length(x)), lickRateNormalizedAllRats{condition}),'ro-','MarkerSize', 5,'MarkerFaceColor', 'r')
        
    end
end
set(gca,'box','off','TickDir', 'out')
axis square
xlabel('Stimulus intensity (um)')
ylabel('Mean lick /s')
xlim([-5 130])
