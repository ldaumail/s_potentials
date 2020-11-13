             
DataDir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\160609_I\';

ns6_filename = '160609_I_cinterocdrft013.ns6';            
                    
             
clear ext NS_header banks neural 
     % Read in NS Header
     ext          = 'ns6'; 
     NS_Header    = openNSx(strcat(DataDir,ns6_filename),'noread');
    el  = 'eD';
    % get basic info about recorded data
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
    N.neural     = sum(neural); % number of neural channels 
    NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
    Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
    nyq          = Fs/2; 
    r            = Fs/1000; 
    r2           = Fs/15000;

    % counters
    clear nct
    nct = 0;

    tic 
    % process data electrode by electrode
    for e = 1:length(neural)

        if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want

            nct = nct+1;

            % open data for this channel. 
            clear NS DAT
            electrode = sprintf('c:%u',e);
            NS        = openNSx(strcat(DataDir,ns6_filename),electrode,'read','uV');
            DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!


            % preallocate data matrices 
            if nct == 1
                N.samples = length(DAT); 
                sLFP       = zeros(ceil(N.samples/r2),N.neural); % preallocating for downsampled data
            end

            % high pass filter
            
              clear hpc hWn bwb bwa hpMUA
            hpc       = 250; %high pass cutoff
            hWn       = hpc/nyq;
            [bwb,bwa] = butter(4,hWn,'high');
            hpsLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
            
            % decimate analog MUA (aMUA) to get 1kHz samp freq
           % MUA(:,nct) = decimate(lpMUA,r);

            % decimate sLFP to get 20kHz samp freq
            sLFP(:,nct) = decimate(hpsLFP,r2); 

            clear DAT 

        end

    end
    toc
       % Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
    %As the .ns6 data was retrieved using openNSx() and not getLFP(), we need
    %to sort the channels ourselves. With getLFP(), sorting is done
    %automatically
    % sort data from top of electrode to bottom.

       % get indices
    idx = zeros(1,length(NeuralLabels));
    for i = 1:length(NeuralLabels)

       Str  = cell2mat(NeuralLabels(i));
       Key   = 'eD';
       Str(strfind(Str, '%02d')) = [];

       Index = strfind(Str, Key);
       idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);

    end
    Srt_sLFP = nan(length(sLFP(:,1)), length(NeuralLabels));
    %Srt_aMUA = nan(length(MUA(:,1)), length(NeuralLabels));


    Srt_sLFP(:,idx) = sLFP(:,:); 
    %Srt_aMUA(:,idx) = MUA(:,:);
    sortedLabels = NeuralLabels(idx);
    
    %% save filtered data to be able to reuse it in different settings in the future analyses
    

        filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';

        % save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
        RelChan = Srt_sLFP(:,11); %store data of the relevant channel
        save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_sLFP.mat'), 'RelChan','-v7.3');
        
        filtChan = load( strcat('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data\' ,ns6_filename, '_250hzhighpass_15khzdownsamp_sLFP.mat'));
   
    %% Manual Spike Sorting
    
    % we need : thresholds (the same as used on BOSS)
    th1 = 39.25; %in ADC = 157
    th2 = -41.25; %in ADC =-165
    
    waveformLen = 48/r2 ; % samples / 1600 microseconds
    preThPeriod =10/r2 ; % samples/ 333 microseconds
    refracPeriod = 30/r2; %samples/ 1000 microseconds
    
    spikes = nan(waveformLen, length(Srt_sLFP(:,11)));
    spkIdxs = nan(length(Srt_sLFP(:,11)),1);
    for s = 1:length(Srt_sLFP(:,11))
       if (Srt_sLFP(s,11) >=th1 && Srt_sLFP(s-1,11)<th1) 
           %for ref = s:s+100/r2
          %         if Srt_sLFP(ref,11) <=th1
                      spk =Srt_sLFP(s-preThPeriod:s+8/r2+refracPeriod,11);
                      spikes(:,s) = spk(1:waveformLen); 
                      spkIdxs(s) =s;
         % break
          %         end
               
         %  end
           
       else 
           if (Srt_sLFP(s,11) <= th2 && Srt_sLFP(s-1,11)>th2)
             %  for ref = s:s+100/r2
             %          if Srt_sLFP(ref,11) >=th2
                          spk =Srt_sLFP(s-preThPeriod:s+8/r2+refracPeriod,11);
                          spikes(:,s) = spk(1:waveformLen);
                          spkIdxs(s) =s;
             %  break
              %         end
                  
              % end
           end
       end
    end
    
   save('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes.mat', 'spikes','-v7.3') %for data larger than 2GB, use version 7.3
   saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes.mat');
   save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_timestamps.mat'), 'spkIdxs','-v7.3'); %save the spikes timestamps
       
   
   figure();
    plot(Srt_sLFP(1:400,11))
    
 
    spikes = spikes(:,~all(isnan(spikes)));
    
    x = linspace(0,1.6,24);
    figure();
    plot(x, 4*spikes)
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    
  %%  Perform a PCA on the data
       %PCA
cenSpikes = spikes - mean(spikes,1); %center the data on the mean

[coeff, score, latent] = pca(cenSpikes'); 
    

%% plot PCA results
x = score(:,1);
y = score(:,2);
z = score(:,3);
figure();
plot3(x,y,z,'.')
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')
grid on


%% plot 10 first components
 
x = linspace(0,1.6,24);
figure();
for i =1:10
  subplot(2,5,i)
  plot(x,coeff(:,i))
 
    xlabel('Time (ms)')
  title(sprintf('PCA %d',i))
end 

%% K means clustering algorithm

idx = kmeans(score,7);
save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx.mat'), 'idx','-v7.3'); %save the spikes timestamps
  
x = linspace(0,1.6,24);
figure();
for i =1:length(unique(idx))
    subplot(1,length(unique(idx)),i)
    plot(x, 4*spikes(:, idx == i)) 
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    title(sprintf('Cluster # %d, n = %d',[i, numel(find(idx ==i))]))
end


%% plot PCA results with K-means
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink



figure();
for i =1:length(unique(idx))
x = score( idx == i,1);
y = score( idx == i,2);
z = score( idx == i,3);

plot3(x,y,z,'.','Color',col(i,:))
hold on
end
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')
grid on


%% Plot spikes with K means clusters colors
 
x = linspace(0,1.6,24);
figure();
for i =1:length(unique(idx))
   h= subplot(1,length(unique(idx)),i)
    plot(x, 4*spikes(:, idx == i),'Color',col(i,:)) 
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    title(sprintf('Cluster # %d, n = %d',[i, numel(find(idx ==i))]))
    set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
end

%% Look at where spikes from the clustering analysis fall according to stimulation onsets and offsets



trialIndexDir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [trialIndexDir, 'selected_trials_idx']);


STIMFileName = '160609_I_p01_uclust1_cinterocdrft_stab_fft_sig';
STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',STIMFileName]);

 %1)get the stim onset times for all trials of the given penetration
    STIM_onsets = STIM_file.STIM.photo_on;
 %only keep the selected trials onsets
    selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(2).idx));
    selected_STIM_onsets = selected_STIM_onsets(1:4:end);
  
 %need to scale up the triggering onset and offset times
 %for data sampled af FS = 15kHz, ==> 15 X more than FS
 %at 1000 Hz
 pre2  = -7500;
 post2 = 22500;

 
 trigSpikeTimes  = trigData(spkIdxs,floor(selected_STIM_onsets./2),-pre2,post2); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units
 trigClustIdxs  = squeeze(trigData(idx,floor(selected_STIM_onsets./2),-pre2,post2)); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units

 trigsLFP  = trigData(Srt_sLFP(:,11),floor(selected_STIM_onsets./2),-pre2,post2); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units

 
 
 %% make raster plot
 filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
 ns6_filename = '160609_I_cinterocdrft013.ns6';  
 clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx.mat')); %load the spikes clusters timestamps

 %trigger the spike cluster idxs to stimuli onsets and offsets
 trialIndexDir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
 selected_trials_idx = load( [trialIndexDir, 'selected_trials_idx']);


STIMFileName = '160609_I_p01_uclust89_cinterocdrft_stab_fft_sig'; %this session has also isolated units 160609_I_p01_uclust1_cinterocdrft_stab_fft_sig, 160609_I_p01_uclust64_cinterocdrft_stab_fft_sig, and 160609_I_p01_uclust89_cinterocdrft_stab_fft_sig
STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',STIMFileName]);

 %1)get the stim onset times for all trials of the given penetration
    STIM_onsets = STIM_file.STIM.photo_on;
 %only keep the selected trials onsets
    selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(5).idx));
    selected_STIM_onsets = selected_STIM_onsets(1:4:end);
  
 %need to scale up the triggering onset and offset times
 %for data sampled af FS = 15kHz, ==> 15 X more than FS
 %at 1000 Hz
 pre2  = -7500;
 post2 = 22500;
 
 trigClustIdxs  = squeeze(trigData(clustIdx.idx,floor(selected_STIM_onsets./2),-pre2,post2)); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units

 %raster plot of each cluster
 figure();
 for cl =1:length(unique(clustIdx.idx(~isnan(clustIdx.idx))))
     h = subplot(2,4,cl);
     spkcnt =0;
     for tr =1:length(trigClustIdxs(1,:))
         spikeTimes =find(trigClustIdxs(:,tr)==cl)';
         x = repmat(spikeTimes,3,1);
         y = nan(size(x));

         if ~isempty(y)
             y(1,:) = tr-1;
             y(2,:) = tr;
         end
     plot(x/15-500,y,'Color','k')
     hold on
     spkcnt = spkcnt+length(spikeTimes);
     end
     title(sprintf('Cluster # %d, n = %d',[cl, spkcnt]))
     set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      xlabel('Time (ms)')
      ylabel('Trial number')
      
 end
 
 %plot corresponding clusters spikes to make sure they are the same as
 %precedently
 
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes.mat');
x = linspace(0,1.6,24);
figure();
for i =1:length(unique(clustIdx.idx(~isnan(clustIdx.idx))))
   h= subplot(1,length(unique(clustIdx.idx(~isnan(clustIdx.idx)))),i);
    plot(x, 4*saved_spikes.spikes(:, clustIdx.idx == i),'Color',col(i,:)) 
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    title(sprintf('Cluster # %d, n = %d',[i, numel(find(clustIdx.idx ==i))]))
    set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
end

% recover PCA results to compare with previous clusters
cenSpikes = saved_spikes.spikes - mean(saved_spikes.spikes,1); %center the data on the mean

[coeff, score, latent] = pca(cenSpikes'); 
figure();
for i =1:length(unique(clustIdx.idx(~isnan(clustIdx.idx))))
x = score( clustIdx.idx == i,1);
y = score( clustIdx.idx == i,2);
z = score( clustIdx.idx == i,3);

plot3(x,y,z,'.','Color',col(i,:))
hold on
end
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')
grid on
 
%raster plot of 1 cluster against every other cluster
%compare cluster 2 with all other clusters
 figure();
 for cl =1:length(unique(clustIdx.idx(~isnan(clustIdx.idx))))
     
     if ~(cl == 2)
         h = subplot(2,4,cl);
         spkcnt =0;
         for tr =1:length(trigClustIdxs(1,:))

             spikeTimes =find(trigClustIdxs(:,tr)==cl)';
             x = repmat(spikeTimes,3,1);
             y = nan(size(x));

             if ~isempty(y)
                 y(1,:) = tr-1;
                 y(2,:) = tr;
             end
         plot(x/15-500,y,'Color','k')
         hold on
         clust2spikeTimes =find(trigClustIdxs(:,tr)==2)';
         clust2x = repmat(clust2spikeTimes,3,1);
         clust2y = nan(size(clust2x));

             if ~isempty(clust2y)
                 clust2y(1,:) = tr-1;
                 clust2y(2,:) = tr;
             end
         plot(clust2x/15 -500, clust2y, 'Color','r')
         hold on
         spkcnt = spkcnt+length(spikeTimes);
         end
     end
     title(sprintf('Cluster # %d vs cluster # 2, n = %d',[cl, spkcnt]))
     set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      xlabel('Time (ms)')
      ylabel('Trial number')
      
 end

    % calculate CSD before triggering to trials OR on the trial data BUT not on
    % the mean LFP. 
%{
    CSD = mod_iCSD(Srt_sLFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                            % feed in units of microV and get back units of
                            % nA/mm^3
    % pad array if you want to keep the matrix the same size on the channel
    % dimension as the other matrices

    CSD = padarray(CSD,[0 1],NaN,'replicate');
     % trigger the neural data to the event codes of interest
    pre   = -500;
    post  = 1500; 

    %need to scale up the triggering onset and offset times
    %for data sampled af FS = 15kHz, ==> 15 X more than FS
    %at 1000 Hz
    pre2  = -7500;
    post2 = 22500;

  %  STIM_file.STIM.aMUA  = trigData(Srt_aMUA,floor(selected_STIM_onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units
    STIM_file.STIM.sCSD  = trigData(CSD,floor(selected_STIM_onsets./2),-pre2,post2); 
    %}
