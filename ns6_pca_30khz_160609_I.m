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
    r2           = Fs/30000;

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
    
    
    
    filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';

    % save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
    RelChan = Srt_sLFP(:,11); %store data of the relevant channel
    save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_sLFP.mat'), 'RelChan','-v7.3');

    %%
    ns6_filename = '160609_I_cinterocdrft013.ns6';            
                    
    filtChan = load( strcat('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data\' ,ns6_filename, '_250hzhighpass_30khz_sLFP.mat'));
    
    %% Manual Spike Sorting
    Srt_sLFP = filtChan.RelChan;
    % we need : thresholds (the same as used on BOSS)
    th1 = 39.25; %in ADC = 157
    th2 = -41.25; %in ADC =-165
    
    r2 =1;
    waveformLen = 100/r2 ; % samples / 1600 microseconds
    preThPeriod =40/r2 ; % samples/ 333 microseconds
    refracPeriod = 30/r2; %samples/ 1000 microseconds
    
    spikes = nan(waveformLen, length(Srt_sLFP));
    spkIdxs = nan(length(Srt_sLFP),1);
    spkLength = nan(length(Srt_sLFP),1);
    for s = 1:length(Srt_sLFP(:,1))
       if (Srt_sLFP(s,1) >=th1 && Srt_sLFP(s-1,1)<th1) 
           %for ref = s:s+100/r2
          %         if Srt_sLFP(ref,11) <=th1
                      spk =Srt_sLFP(s-preThPeriod:s+8/r2+refracPeriod,1);
                      if all(isnan(spikes(:,s-15:s-1)))
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
                          if all(isnan(spikes(:,s-15:s-1)))
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
    
   save('C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat', 'spikes','-v7.3') %for data larger than 2GB, use version 7.3
   saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
   filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';

   save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khzdownsamp_spikes_timestamps.mat'), 'spkIdxs','-v7.3'); %save the spikes timestamps
   
   spikes = spikes(:,~all(isnan(spikes))); 
    x = linspace(0,2.6,79);
    figure();
    plot(x,4*spikes(:,1:1000))
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    
%% Perform PCA
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
 spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));
 
 cenSpikes = spikes - nanmean(spikes,1); %center the data on the mean

[coeff, score, latent] = pca(cenSpikes'); 
    

%% plot PCA results
x = score(:,1);
y = score(:,3);
z = score(:,2);
figure();
plot3(x,y,z,'.')
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
ylabel('PCA3')
zlabel('PCA2')
grid on


%% K means clustering algorithm

idx = kmeans(score,7);
filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
ns6_filename = '160609_I_cinterocdrft013.ns6';            
save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_spikes_7clustidx_11272020.mat'), 'idx','-v7.3'); %save the spikes timestamps
  
x = linspace(0,2.6,79);
figure();
for i =1:length(unique(idx))
    subplot(1,length(unique(idx)),i)
    plot(x, 4*spikes(:, idx == i)) 
    ylabel('Voltage (microVolts)')
    xlabel('Time (ms)')
    title(sprintf('Cluster # %d, n = %d',[i, numel(find(idx ==i))]))
end
 

%% %% plot PCA results with K-means
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
z = score( idx == i,2);
y = score( idx == i,3);

plot3(x,y,z,'.','Color',col(i,:))
hold on
end
hold on
yline(0)
hold on
xline(0)
xlabel('PCA1')
zlabel('PCA2')
ylabel('PCA3')

grid on


%% Plot spikes with K means clusters colors
  ns6_filename = '160609_I_cinterocdrft013.ns6';            
    
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
 spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));

 filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data'; 
clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_spikes_7clustidx_11272020.mat')); %load the spike cluster IDs
idx = clustIdx.idx;

col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink


x = linspace(0,2.6,79);
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

%% Compare clusters to spike waveforms of single units isolated in this cluster computing the squared error between the single units spike waveforms and the clusters' spike waveforms

%load the spikes and the cluster index of each spike
 saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
 spikes = saved_spikes.spikes(:,~all(isnan(saved_spikes.spikes)));

 
 %Load stim file of a single unit containing the spike waveform.
 
%this session contains isolated units 160609_I_p01_uclust1_cinterocdrft_stab_fft_sig, 160609_I_p01_uclust64_cinterocdrft_stab_fft_sig, and 160609_I_p01_uclust89_cinterocdrft_stab_fft_sig
STIM_file1 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust1_cinterocdrft_stab_fft_sig');
STIM_file2 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust64_cinterocdrft_stab_fft_sig');
STIM_file3 = load('C:\Users\daumail\Documents\LGN_data\single_units\160609_I_p01_uclust89_cinterocdrft_stab_fft_sig');


waveform(:,1) = mean(squeeze(STIM_file1.STIM.wf.waveForms(:,:,11,:)),1);
waveform(:,2) = mean(squeeze(STIM_file2.STIM.wf.waveForms(:,:,11,:)),1);
waveform(:,3) = mean(squeeze(STIM_file3.STIM.wf.waveForms(:,:,11,:)),1);

%{
figure()
plot(waveform(15:30,:), 'Color', 'b')
hold on
plot(4*saved_spikes.spikes(1:15,1:10000), 'Color', 'r')
%}
  th1 = 39.25; %in ADC = 157
  th2 = -41.25; %in ADC =-165
  
 
[~,max_location] = max(abs(waveform),[],1);
[~,max_ClustLoc] = max(abs(spikes(20:69,:)),[],1);
sqdError = nan(21, length(spikes(1,:)), length(waveform(1,:)));
for n = 1:length(waveform(1,:))
    for i =1:length(spikes(1,:))

                      % clustwf = (spikes(max_ClustLoc(i)-10+20:max_ClustLoc(i)+10+20,i)-min(spikes(max_ClustLoc(i)-10+20:max_ClustLoc(i)+10+20,i)))./(max(spikes(max_ClustLoc(i)-10+20:max_ClustLoc(i)+10+20,i))-min(spikes(max_ClustLoc(i)-10+10:max_ClustLoc(i)+10+10,i)));
                      % suwf = (waveform(max_location(n)-10:max_location(n)+10,n) - min(waveform(max_location(n)-10:max_location(n)+10,n)))./(max(waveform(max_location(n)-10:max_location(n)+10,n))-min(waveform(max_location(n)-10:max_location(n)+10,n)));
                      % sqdError(:,i,n) = (clustwf - suwf).^2;
                      if max_ClustLoc(i)+10+20 <=79
                        sqdError(:,i,n) = (4*spikes(max_ClustLoc(i)-10+20:max_ClustLoc(i)+10+20,i) - waveform(max_location(n)-10:max_location(n)+10,n)).^2;
                      else
                        if max_ClustLoc(i)+10+20 > 79
                          sqdError(1:length(spikes(max_ClustLoc(i)-10+20:end,i)),i,n) = (4*spikes(max_ClustLoc(i)-10+20:end,i) - waveform(max_location(n)-10:max_location(n)+floor((length(spikes(max_ClustLoc(i)-10+20:end,i)))/2-1),n)).^2;
                        end
                      end
    end
end

figure();
plot(clustwf)
hold on
plot(suwf)

figure();
plot(4*spikes(max_ClustLoc(i)-10+20:max_ClustLoc(i)+10+20,i))
hold on
plot(waveform(max_location(n)-10:max_location(n)+10,n))
            %{       
       
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
   SSE(i,n) = sum(sqdError(~isnan(sqdError(:,i,n)),i,n))./length(find(~isnan(sqdError(:,i,n))));  %normalize by number of non nan values, as spikes could be compared over a -5:+5 sample nb for a 15 kHz dataset with a -10:+10 in a 30kHz dataset
    end
end

%%plot for each cluster the median SSE for each single unit spike
%%waveform
filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';
ns6_filename = '160609_I_cinterocdrft013.ns6';  
clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_spikes_7clustidx_11272020.mat')); %load the spikes clusters timestamps

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
    %ylim([0.015 0.375])
     ylim([0 110000])
    set(gca,'box','off')
    set(gca, 'linewidth',2)
    if cl > 1
      ax1 = gca;                   
      ax1.YAxis.Visible = 'off'; 
    end
  ylabel('Median SSE (normalized)')
   if cl ==1
      
  xlabel('Single Unit #')
   end
   
  title(sprintf('Cluster %d',cl)) 
  
end
%legend('median','95%CI', 'Location', 'bestoutside')



%% Inter spike time difference histograms to merge clusters

saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
spikes = saved_spikes.spikes;

filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data'; 
ns6_filename = '160609_I_cinterocdrft013.ns6';  

clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_spikes_7clustidx_11272020.mat')); %load the spike cluster IDs
idx = clustIdx.idx;

%rebuild vector with space between spike indices to find times
spkTimes = find(~all(isnan(spikes)));

%%insert cluster indices at the right times
%new vector
clustTimes = nan(1,length(spikes(1,:)));
clustTimes(spkTimes) = idx;

     %% merge clusters
r =30;
%based on PCA and SSE analysis, we could merge clusters 2,3, and 5. 
%however, with further analysis, we can try to see the histograms of
%interspike intervals

%spike indices of each cluster
clust2Idxs =find(clustTimes ==2);
clust3Idxs =find(clustTimes ==3);
clust5Idxs =find(clustTimes ==5);
clust6Idxs =find(clustTimes ==6);


allClustsIdxs = sort([find(clustTimes ==2),find(clustTimes ==3), find(clustTimes ==5),find(clustTimes ==6)]);

%compute interspike interval for each cluster
clust2diffs = (clust2Idxs(2:end)-clust2Idxs(1:end-1))./r;
clust3diffs = (clust3Idxs(2:end)-clust3Idxs(1:end-1))./r;
clust5diffs = (clust5Idxs(2:end)-clust5Idxs(1:end-1))./r;
clust6diffs = (clust6Idxs(2:end)-clust6Idxs(1:end-1))./r;

allClustsdiffs = (allClustsIdxs(2:end) - allClustsIdxs(1:end-1))./r;

col(1,:) =[86/255 86/255 86/255] ; %--dark grey 
col(2,:) = [251/255 154/255 153/255]; % -- red
col(3,:) = [146/255 197/255 222/255]; % -- blue
col(4,:) =[194/255 165/255 207/255] ; %--purple
col(5,:) = [253/255 174/255 97/255]; % -- orange
col(6,:) = [166/255 219/255 160/255]; % -- green
col(7,:) = [238/255 58/255 104/255]; % -- pink


figure();
   
subplot( 1,5,1) 
histogram(clust2diffs,10000, 'FaceColor', col(2,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
xlabel('Time difference (ms)')
ylabel('Count')
title('Interspike Interval for Cluster 2')

subplot( 1,5,2) 
histogram(clust3diffs,100000, 'FaceColor', col(3,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for Cluster 3')

subplot( 1,5,3) 
histogram(clust5diffs,10000, 'FaceColor', col(5,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for Cluster 5')

subplot( 1,5,4) 
histogram(clust6diffs,100000, 'FaceColor', col(6,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for Cluster 6')


subplot( 1,5,5) 
histogram(allClustsdiffs,10000, 'FaceColor', col(7,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for clusters 2,3,5 and 6 combined')

%%look at clusters 1, and 7
clust1Idxs =find(clustTimes ==1);
clust7Idxs =find(clustTimes ==7);

allClustsIdxs = sort([find(clustTimes ==1),find(clustTimes ==7)]);

%compute interspike interval for each cluster
clust1diffs = (clust1Idxs(2:end)-clust1Idxs(1:end-1))./r;
clust7diffs = (clust7Idxs(2:end)-clust7Idxs(1:end-1))./r;

allClustsdiffs = (allClustsIdxs(2:end) - allClustsIdxs(1:end-1))./r;


figure();
   
subplot( 1,3,1) 
histogram(clust1diffs,10000, 'FaceColor', col(1,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
xlabel('Time difference (ms)')
ylabel('Count')
title('Interspike Interval for Cluster 1')

subplot( 1,3,2) 
histogram(clust7diffs,100000, 'FaceColor', col(7,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for Cluster 7')

subplot( 1,3,3) 
histogram(allClustsdiffs,100000, 'FaceColor', col(2,:));
xlim([0 5])
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
title('Interspike Interval for Cluster 1 and 7')
%% Spiking activity during trials

    %% raster plots
%lets suppose cluster 4 might be s-potentials
%lets suppose clusters 2,3,5 are single unit 2 (uclust64) that activates according
%to cluster 4
saved_spikes = load( 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\160609_I_cinterocdrft013_spikes_s30khz_250hzbutter_dual_spikes.mat');
spikes = saved_spikes.spikes;

filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data'; 
ns6_filename = '160609_I_cinterocdrft013.ns6';  

clustIdx = load(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_spikes_7clustidx_11272020.mat')); %load the spike cluster IDs
idx = clustIdx.idx;
clustTimes = nan(1,length(spikes(1,:)));
clustTimes(spkTimes) = idx;

UnitIdxs = sort([find(clustTimes ==2),find(clustTimes ==3), find(clustTimes ==5)]);
UnitBin = nan(length(spikes(1,:)),1);
UnitBin(UnitIdxs) = ones(length(UnitIdxs),1);

SpotIdxs = sort(find(clustTimes ==4));
SpotBin = nan(length(spikes(1,:)),1);
SpotBin(SpotIdxs) = ones(length(SpotIdxs),1);

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
 pre  = -15000;
 post = 45000;
 
 
 trigUnitIdxs  = squeeze(trigData(UnitBin,floor(selected_STIM_onsets./1),-pre,post)); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units
 trigSpotIdxs  = squeeze(trigData(SpotBin,floor(selected_STIM_onsets./1),-pre,post)); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units

 %raster plot of unit
 h =figure();
     % h = subplot(2,4,cl);
     spkcnt =0;
     for tr =1:length(trigUnitIdxs(1,:))
         spikeTimes =find(trigUnitIdxs(:,tr)==1)';
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
      spotcnt =0;
     for tr =1:length(trigSpotIdxs(1,:))
         spikeTimes =find(trigSpotIdxs(:,tr)==1)';
         x = repmat(spikeTimes,3,1);
         y = nan(size(x));

         if ~isempty(y)
             y(1,:) = tr-1;
             y(2,:) = tr;
         end
     plot(x/15-500,y,'Color','r')
     hold on
     spotcnt = spotcnt+length(spikeTimes);
     end
     title(sprintf('Unit #2, n = %d with cluster 4, n =%d', [spkcnt, spotcnt]))
     set(h,'position',get(h,'position').*[1 1 1.15 1])
       
      set(gca,'box','off')
      set(gca, 'linewidth',2)

      xlabel('Time (ms)')
      ylabel('Trial number')
      



    %% Convolve Poisson-like kernel
 Fs = 30000;
 r = 1;
 k         = jnm_kernel( 'psp', (20/1000) * (Fs/r) );
 
 trigSdf = struct();
        UnitIdxs = sort([find(clustTimes ==2),find(clustTimes ==3), find(clustTimes ==5)]);
        UnitBin = zeros(length(spikes(1,:)),1);
        UnitBin(UnitIdxs) = ones(length(UnitIdxs),1);
        SpotIdxs = sort(find(clustTimes ==4));
        SpotBin = zeros(length(spikes(1,:)),1);
        SpotBin(SpotIdxs) = ones(length(SpotIdxs),1);


        sdf = conv(UnitBin,k,'same') * Fs/r;
        SpotSdf = conv(SpotBin,k,'same') * Fs/r;
        %suaN = sprintf('sua%d', cl);
         trigSdf.sua = squeeze(trigData(sdf, floor(selected_STIM_onsets./2),-pre,post));

x =linspace(-0.5,1.833,70001);
figure()
plot(x,sdf(780000:850000),'Color','k')
hold on
plot(x,SpotSdf(780000:850000), 'Color','r')
set(gca, 'Box', 'off')
set(gca, 'linewidth',2)
legend('Single Unit', 'Cluster 4')
xlabel('Time (s)')
ylabel('Spike rate')

    
    %% plot
