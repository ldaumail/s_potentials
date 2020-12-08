             
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
    
    %% save filtered data to be able to reuse it in different settings in future analyses
    

        filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';

        % save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
        RelChan = Srt_sLFP(:,11); %store data of the relevant channel
        save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_sLFP.mat'), 'RelChan','-v7.3');
        
    %% Manual Spike Sorting
    

    ns6_filename = '160609_I_cinterocdrft013.ns6';  
    filtChan = load( strcat('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data\' ,ns6_filename, '_250hzhighpass_15khzdownsamp_sLFP.mat'));
    Srt_sLFP = filtChan.RelChan;
    % we need : thresholds (the same as used on BOSS)
    th1 = 39.25; %in ADC = 157
    th2 = -41.25; %in ADC =-165
    r2 =2;
    waveformLen = 100/r2 ; % samples /? 
    preThPeriod =40/r2 ; % samples/ ? microseconds
    refracPeriod = 30/r2; %samples/ ? microseconds
    
    spikes = nan(waveformLen, length(Srt_sLFP));
    spkIdxs = nan(length(Srt_sLFP),1);
    spkLength = nan(length(Srt_sLFP),1);
    for s = preThPeriod+1:length(Srt_sLFP(:,1))
        spk =0;
       if (Srt_sLFP(s,1) >=th1 && Srt_sLFP(s-1,1)<th1) 
           %for ref = s:s+100/r2
          %         if Srt_sLFP(ref,11) <=th1
                      spk =Srt_sLFP(s-preThPeriod:s+8/r2+refracPeriod,1);
                      if all(isnan(spikes(:,s-8:s-1)))
                      spikes(1:length(spk),s) = spk(:); 
                      spkIdxs(s) =s;
                      end
         % break
          %         end
               
         %  end
           
       else 
           if (Srt_sLFP(s,1) <= th2 && Srt_sLFP(s-1,1)>th2)
             %  for ref = s:s+100/r2
             %          if Srt_sLFP(ref,11) >=th2
                          spk =Srt_sLFP(s-preThPeriod:s+8/r2+refracPeriod,1);
                          if all(isnan(spikes(:,s-8:s-1)))
                          spikes(1:length(spk),s) = spk(:);
                          spkIdxs(s) =s;
                          end
             %  break
              %         end
                  
              % end
           end
           spkLength(s) =length(spk);
       end
       
    end

    minLength = min(spkLength);
    spikes = spikes(1:minLength,:);
    
   save('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes_11272020.mat', 'spikes','-v7.3') %for data larger than 2GB, use version 7.3
   saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes_11272020.mat');
   
   filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
   save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_timestamps_11272020.mat'), 'spkIdxs','-v7.3'); %save the spikes timestamps
       
   
   figure();
    plot(Srt_sLFP(1:400,11))
    
  spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));
    %spikes = spikes(:,~all(isnan(spikes)));
    
    x = linspace(0,2.6,40);
    figure();
    plot(x, 4*spikes(:,1:100))
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
save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx_11272020.mat'), 'idx','-v7.3'); %save the spikes timestamps
  
x = linspace(0,2.6,40);
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
y = score( idx == i,3);
z = score( idx == i,2);

plot3(x,y,z,'.','Color',col(i,:))
hold on
end
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA3')
zlabel('PCA2')
grid on


%% Plot spikes with K means clusters colors
 
x = linspace(0,2.6,40);
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

ns6_filename = '160609_I_cinterocdrft013.ns6';            
    
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes_11272020.mat');
 spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));

 filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data'; 
clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx_11272020.mat')); %load the spike cluster IDs
idx = clustIdx.idx;
%idx = idx(~isnan(idx)); 
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink


x = linspace(0,2.6,40);
figure();
for i =1:length(unique(idx))
    meanSpk = mean(4*spikes(:, idx == i),2);
    ciLow = meanSpk - 1.96*std(4*spikes(:, idx == i),0,2)/sqrt(length(spikes(1, idx == i)));
    ciHigh = meanSpk + 1.96*std(4*spikes(:, idx == i),0,2)/sqrt(length(spikes(1, idx == i)));
   h= subplot(1,length(unique(idx)),i);
    plot(x,meanSpk,'Color',col(i,:)) 
   % plot(x,4*spikes(:, idx == i),'Color',col(i,:))
    hold on
    h1 = ciplot(ciLow, ciHigh,x, col(i,:),0.5);
    set(h1, 'edgecolor','none')
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    title(sprintf('Cluster # %d, n = %d',[i, numel(find(idx ==i))]))
    set(h,'position',get(h,'position').*[1 1 1.15 1])
      
      set(gca,'box','off')
      
      set(gca, 'linewidth',2)
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
 
 %% plot corresponding clusters spikes to make sure they are the same as
 %precedently
 filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
 ns6_filename = '160609_I_cinterocdrft013.ns6';  
 clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx.mat')); %load the spikes clusters timestamps
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes.mat');
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink


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
 
 %% convert spike times into spike rate with Poisson-like sdf convolution
 
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

 %sdftm     = [-0.3*Fs/r: 0.15*Fs/r + max(diff(TP,[],2)/r)];
 %[spk,spktm,~] = intersect(sdftm,x,'stable') ;
 Fs = 30000;
 r = 2;
 k         = jnm_kernel( 'psp', (20/1000) * (Fs/r) );

 trigSdf = struct();
  for cl =1:length(unique(clustIdx.idx(~isnan(clustIdx.idx))))
         
         sua = zeros(length(clustIdx.idx),1);
         sua(find(clustIdx.idx==cl))=1;
         sdf = conv(sua,k,'same') * Fs/r;
         suaN = sprintf('sua%d', cl);
         trigSdf.(suaN) = squeeze(trigData(sdf, floor(selected_STIM_onsets./2),-pre2,post2));
  end
  
figure()
plot(mean(trigSdf.sua7,2))

%% Compare clusters to spike waveforms of single units isolated in this cluster computing the squared error between the single units spike waveforms and the clusters' spike waveforms

%load the spikes and the cluster index of each spike
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s15khz_250hzbutter_dual_spikes_11272020.mat');
 spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));

 
 %Load stim file of a single unit containing the spike waveform.
 
%this session contains isolated units 160609_I_p01_uclust1_cinterocdrft_stab_fft_sig, 160609_I_p01_uclust64_cinterocdrft_stab_fft_sig, and 160609_I_p01_uclust89_cinterocdrft_stab_fft_sig
STIM_file1 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust1_cinterocdrft_stab_fft_sig');
STIM_file2 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust64_cinterocdrft_stab_fft_sig');
STIM_file3 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust89_cinterocdrft_stab_fft_sig');


waveform(:,1) = decimate(mean(squeeze(STIM_file1.STIM.wf.waveForms(:,:,11,:)),1),2);
waveform(:,2) = decimate(mean(squeeze(STIM_file2.STIM.wf.waveForms(:,:,11,:)),1),2);
waveform(:,3) = decimate(mean(squeeze(STIM_file3.STIM.wf.waveForms(:,:,11,:)),1),2);

%{
figure()
plot(waveform(15:30,:), 'Color', 'b')
hold on
plot(4*saved_spikes.spikes(1:15,1:10000), 'Color', 'r')

  th1 = 39.25; %in ADC = 157
  th2 = -41.25; %in ADC =-165
  
 
[~,max_location] = max(abs(waveform),[],1);
sqdError = nan(11, length(spikes(1,:)), length(waveform(1,:)));
for n = 1:length(waveform(1,:))
   % for i =1:length(spikes(1,:))

   for i =1:10
       
        for g =2:24 %g for growth
            
            

            if spikes(g,i) >=th1 && spikes(g-1,i) <th1
               % plot(4*saved_spikes.spikes(max_location(i)-5:max_location(i)+5,i), 'Color', 'r')
               % hold on
               for d =1:24 %decay 
                   if (g+d<=24 && spikes(g+d,i)< spikes(g+d-1,i)) && (g+d-1+5 <= 24 && g+d-1-5 >=1)
                       %normalize spike waveforms
                       clustwf = (spikes(g+d-5-1:g+d+5-1,i)-min(spikes(g+d-5-1:g+d+5-1,i)))./(max(spikes(g+d-5-1:g+d+5-1,i))-min(spikes(g+d-5-1:g+d+5-1,i)));
                       suwf = (waveform(max_location(n)-5:max_location(n)+5,n) - min(waveform(max_location(n)-5:max_location(n)+5,n)))./(max(waveform(max_location(n)-5:max_location(n)+5,n))-min(waveform(max_location(n)-5:max_location(n)+5,n)));
                       sqdError(1:11,i,n) = (clustwf - suwf).^2;
                   
                   else
                        if (g+d<=24 && spikes(g+d,i)< spikes(g+d-1,i)) && g+d-1+5 > 24
                       clustwf = (spikes(g+d-5-1:end,i)-min(spikes(g+d-5-1:end,i)))./(max(spikes(g+d-5-1:end,i))-min(spikes(g+d-5-1:end,i)));
                       suwf = (waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n) - min(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n)))./(max(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n))-min(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n)));
     
                       sqdError(1:length(spikes(g+d-5-1:end,i)),i,n) = (clustwf - suwf).^2;
                       
                        
                           else
                               if (g+d<=24 && spikes(g+d,i)< spikes(g+d-1,i)) &&  g+d-1-5 <1
                                clustwf = (spikes(1:g+d+5-1,i)-min(spikes(1:g+d+5-1,i)))./(max(spikes(1:g+d+5-1,i))-min(spikes(1:g+d+5-1,i)));
                                suwf = (waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n) - min(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n)))./(max(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n))-min(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n)));
                          
                                sqdError(1:length(spikes(1:g+d+5-1,i)),i,n) = (clustwf - suwf).^2;

                               end
                        end
                   end
                   break
               end
               
            else
                if spikes(g,i) <=th2 && spikes(g-1,i) >th2 
                    clear d
                    for d =1:24 %decay
                       if (g+d<=24 && spikes(g+d,i)> spikes(g+d-1,i))  && (g+d-1+5 <= 24 && g+d-1-5 >=1)
                           clustwf = (spikes(g+d-5-1:g+d+5-1,i)-min(spikes(g+d-5-1:g+d+5-1,i)))./(max(spikes(g+d-5-1:g+d+5-1,i))-min(spikes(g+d-5-1:g+d+5-1,i)));
                           suwf = (waveform(max_location(n)-5:max_location(n)+5,n) - min(waveform(max_location(n)-5:max_location(n)+5,n)))./(max(waveform(max_location(n)-5:max_location(n)+5,n))-min(waveform(max_location(n)-5:max_location(n)+5,n)));
                          sqdError(1:11,i,n) = (clustwf - suwf).^2;
                            
                            else
                        if (g+d<=24 && spikes(g+d,i)> spikes(g+d-1,i)) && g+d-1+5 > 24 
                       clustwf = (spikes(g+d-5-1:end,i)-min(spikes(g+d-5-1:end,i)))./(max(spikes(g+d-5-1:end,i))-min(spikes(g+d-5-1:end,i)));
                       suwf = (waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n) - min(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n)))./(max(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n))-min(waveform(max_location(n)-5:max_location(n)+length(spikes(g+d-1:end,i))-1,n)));
     
                       sqdError(1:length(spikes(g+d-5-1:end,i)),i,n) = (clustwf - suwf).^2;
                        
                           else
                               if (g+d<=24 && spikes(g+d,i)> spikes(g+d-1,i)) &&  g+d-1-5 <1
                                   clustwf = (spikes(1:g+d+5-1,i)-min(spikes(1:g+d+5-1,i)))./(max(spikes(1:g+d+5-1,i))-min(spikes(1:g+d+5-1,i)));
                                   suwf = (waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n) - min(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n)))./(max(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n))-min(waveform(max_location(n)-(length(spikes(1:g+d-1,i))-1):max_location(n)+5,n)));
                                   sqdError(1:length(spikes(1:g+d+5-1,i)),i,n) = (clustwf - suwf).^2;

                              
                               end
                        end
                       end
                       break
                    end
                end
                
                

            end
            
        end
        
    end

end
%}

 th1 = 39.25; %in ADC = 157
  th2 = -41.25; %in ADC =-165
  
 
[~,max_location] = max(abs(waveform),[],1);
[~,max_ClustLoc] = max(abs(spikes(10:34,:)),[],1);
sqdError = nan(11, length(spikes(1,:)), length(waveform(1,:)));
for n = 1:length(waveform(1,:))
    for i =1:length(spikes(1,:))

                       %clustwf = (spikes(max_ClustLoc(i)-5+10:max_ClustLoc(i)+5+10,i)-min(spikes(max_ClustLoc(i)-5+10:max_ClustLoc(i)+5+10,i)))./(max(spikes(max_ClustLoc(i)-5+10:max_ClustLoc(i)+5+10,i))-min(spikes(max_ClustLoc(i)-5+10:max_ClustLoc(i)+5+10,i)));
                       %suwf = (waveform(max_location(n)-5:max_location(n)+5,n) - min(waveform(max_location(n)-5:max_location(n)+5,n)))./(max(waveform(max_location(n)-5:max_location(n)+5,n))-min(waveform(max_location(n)-5:max_location(n)+5,n)));
                      % sqdError(1:11,i,n) = (clustwf - suwf).^2;
                        sqdError(1:11,i,n) = (4*spikes(max_ClustLoc(i)-5+10:max_ClustLoc(i)+5+10,i) - waveform(max_location(n)-5:max_location(n)+5,n)).^2;
    end
end

%{
figure();
plot(suwf)
hold on
plot(clustwf)
hold on
%}
%% compute the sum of the squared error, keep the index of the smallest summed squared error (SSE)


SSE = nan(length(spikes(1,:)), length(waveform(1,:)));
for i =1:length(spikes(1,:))
    for n =1:length(waveform(1,:))
   SSE(i,n) = sum(sqdError(~isnan(sqdError(:,i,n)),i,n))./length(find(~isnan(sqdError(:,i,n)))); %normalize by number of non nan values, as spikes could be compared over a -5:+5 sample nb for a 15 kHz dataset with a -10:+10 in a 30kHz dataset
    end
end

%%plot for each cluster the median SSE for each single unit spike
%%waveform
filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
ns6_filename = '160609_I_cinterocdrft013.ns6';  
clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_15khzdownsamp_spikes_clustidx_11272020.mat')); %load the spikes clusters timestamps

clust = clustIdx.idx(~isnan(clustIdx.idx));

col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink

figure();
x =1:3;
lowCI = nan(length(waveform(1,:)),1);
highCI = nan(length(waveform(1,:)),1);
medianSSE = nan(length(waveform(1,:)),length(unique(clust)));
   
for cl = 1:length(unique(clust))
    allClustSSE = SSE(clust ==cl,:);
    low_idx = floor(0.5*length(allClustSSE(~isnan(allClustSSE(:,1)),1)) - 1.96*sqrt(0.5*length(allClustSSE(~isnan(allClustSSE(:,1)),1))*(1-0.5)));
    high_idx = ceil(0.5*length(allClustSSE(~isnan(allClustSSE(:,1)),1)) + 1.96*sqrt(0.5*length(allClustSSE(~isnan(allClustSSE(:,1)),1))*(1-0.5)));
 
     for n =1:length(waveform(1,:))
          clustSSE = SSE(clust ==cl,n);
          NotNanSSE =clustSSE(~isnan(clustSSE));
         [~,sortIdx]  = sort(NotNanSSE);

        lowCI(n,1) = NotNanSSE(sortIdx(low_idx));
        highCI(n,1) = NotNanSSE(sortIdx(high_idx));
        medianSSE(n,cl) = median(NotNanSSE);
    end
    subplot(1,7,cl)
    plot(x,medianSSE(:,cl),'-*','Color', col(cl,:))
    hold on
    h1=ciplot(lowCI, highCI,x, col(cl,:),0.5);
    set(h1, 'edgecolor','none')
   % ylim([0.02 0.334])
   ylim([0 100000])
   
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    if cl > 1
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off'; 
    end
  ylabel('Median SSE (microVolts^2)')
   if cl ==1
  xtickformat('%i')
  xlabel('Single Unit #')
   end
   
  title(sprintf('Cluster %d',cl)) 
  
end
%legend('median','95%CI', 'Location', 'bestoutside')

%% Plot the spike waveforms of KiloSort-ed Single Units
STIM_file1 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust1_cinterocdrft_stab_fft_sig');
STIM_file2 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust64_cinterocdrft_stab_fft_sig');
STIM_file3 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust89_cinterocdrft_stab_fft_sig');


waveform(:,1) = mean(squeeze(STIM_file1.STIM.wf.waveForms(:,:,11,:)),1);
waveform(:,2) = mean(squeeze(STIM_file2.STIM.wf.waveForms(:,:,11,:)),1);
waveform(:,3) = mean(squeeze(STIM_file3.STIM.wf.waveForms(:,:,11,:)),1);

semWaveform(:,1) = std(squeeze(STIM_file1.STIM.wf.waveForms(:,:,11,:)),0,1)./sqrt(length(squeeze(STIM_file1.STIM.wf.waveForms(:,:,11,1))));
semWaveform(:,2) = std(squeeze(STIM_file2.STIM.wf.waveForms(:,:,11,:)),0,1)./sqrt(length(squeeze(STIM_file2.STIM.wf.waveForms(:,:,11,1))));
semWaveform(:,3) = std(squeeze(STIM_file3.STIM.wf.waveForms(:,:,11,:)),0,1)./sqrt(length(squeeze(STIM_file3.STIM.wf.waveForms(:,:,11,1))));
col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue

figure();
x =linspace(0,(81/30000)*1000,81);
for i =1:3
    h =subplot(1,3,i)
    plot(x,waveform(:,i), 'Color', col(i,:))
    hold on
    h1 = ciplot(waveform(:,i)-1.96*semWaveform(:,i),waveform(:,i)+1.96*semWaveform(:,i), x, col(i,:), 0.5);
    set(h1, 'edgecolor','none') 
    set(h,'position',get(h,'position').*[1 1 1.15 1])
    set(gca, 'linewidth',2)
    set(gca,'box','off')
    ylabel('Voltage (microVolt)')
    xlabel('time (ms)')
    title(sprintf('Single Unit %d',i))
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
