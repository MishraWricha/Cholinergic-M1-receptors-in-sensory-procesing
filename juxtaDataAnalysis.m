clearvars
clc
close all

%% parameters
allNeurons=[1]; %Add neuron number here 
amplitudes=[0 20 50 100 200];
conditions = {'ACSF','BQCA', 'TD'};
preStimDuration = .1;
postStimDuration = .1;
%% loading And Filtering Data
cd '\\researchdata.anu.edu.au\\...' %specify the directory to load files 
window=0.05;
for neuron=21
    for condition = 1:length(conditions)
        tmp1=dir(['SpikeTimesNeuron-' int2str(neuron) '-' conditions{condition} '*' ])
        
        for i=1:length(tmp1)
            load(tmp1(i).name)
            
            figure
           
            
            
            for amp=1:size(spikeTimes,2)%finds the length of the second dimension of spikeTimes
                for trial=1:size(spikeTimes,1)%finds the length of the first dimension of spikeTimes-30 trials
                    subplot(2,size(spikeTimes,2),amp)%creates subpplot of 1X30 and puts axes in position specified by amp
                    hold on
                    plot(spikeTimes{trial,amp},trial*ones(size(spikeTimes{trial,amp})),'.k')%
                    title(['neuron ' int2str(neuron) conditions{condition} ' session ' int2str(i)])
                    ylabel('No of trials'); xlabel('spike time');
                    xlim([-0.1 0.1])
                end
            end
            hold on
            %% PSTH
            bin = 0.01;
            edges=-0.1:bin:0.1; %-0.4 to 0.6 in steps of 0.05 - 20 bins
            psth=zeros(length(edges)-1,length(amplitudes));%create array of zeros 19*5
            for amp=1:size(spikeTimes,2)
                for trial = 1:size(spikeTimes,1)
                    a = spikeTimes{trial,amp};
                    for i=1:length(edges)-1
                        b = length(find(a>edges(i) & a<edges(i+1)));
                        psth(i,amp) = psth(i,amp)+b;
                        amp1=amp+5;
                        subplot(2,size(spikeTimes,2),amp1)
                    end
                    
                end
                bar(psth(:,amp))
                ylabel('No of neurons'); xlabel('Bin number');
                ylim([0 80])
            end
                   
        end
    end
    
    
    
end
%% neurometric plots 
for neuron =25
    figure
    title(['Neuron ' int2str(neuron)])
    for condition = 1:length(conditions)
        tmp1=dir(['SpikeTimesNeuron-' int2str(neuron) '-' conditions{condition} '*' ]);
        
        hold on
        for i=1:length(tmp1)
            load(tmp1(i).name)
            
            window=0.05;
            hold on
            for amp=1:size(spikeTimes,2)%finds the length of the second dimension of spikeTimes
                for trial=1:size(spikeTimes,1)%finds the length of the first dimension of spikeTimes-30 trials
                    spikeCounts(trial,amp)=length(find(spikeTimes{trial,amp}<window & spikeTimes{trial,amp}>0));
                end
            end
                            
            if condition==1
                plot(amplitudes,mean(spikeCounts),'ko-')
            elseif condition==2
                plot(amplitudes,mean(spikeCounts),'go-')
            elseif condition==3
                plot(amplitudes,mean(spikeCounts),'ro-')
  
            end
            %ylim([0 3]);
        end
    end
end  
