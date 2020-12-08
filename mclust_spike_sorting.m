%This script was started on 10/05/2020

%by Loic Daumail

%The goal is to try isolating S potentials from the LGN data using Mclust
%(MClust-4.0, A.D. Redish et al)


%load packages used to load .ns6 file
npmkdir    = 'C:\Users\daumail\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\daumail\Documents\bootcamp-selected\nbanalysis\'; 
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))
addpath 'C:\Users\daumail\Documents\loic_code'

%ns6 directory 
ns6dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\ns6_selected_units\';

%nev directory 
nevdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\nev_selected_units\';


%ssdir
ssdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\kilosorted_files\';
%ssdir content
selected_name_list = dir(ssdir);  

% dir with saved channels data
savedir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\indiv_chan_data\';
%load file with selected trial indexes for a given penetration with
%penetration name
indexdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\analysis\';
selected_trials_idx = load( [indexdir, 'selected_trials_idx']);

%load selected penetrations files list
%selected_penetrations = 'C:\Users\daumail\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_03032020_corrected\orig_peak_values\all_units\selected_orig_units_filenames';
%penetrations_names = textscan( fopen(strcat(selected_penetrations, '.txt')), '%s', 'HeaderLines', 1);
%penetrations_names = cell2struct(penetrations_names, 'penetrations_names');

for pnt = 1:length(selected_trials_idx.logicals)
 
    if ~isempty(selected_trials_idx.logicals(pnt).idx)
       
        penetration_name = erase(selected_trials_idx.logicals(pnt).penetration, 'matDE50_NDE0');
        STIM_file = load(['C:\Users\daumail\Documents\LGN_data\single_units\',penetration_name]);
    
     %1)get the stim onset times for all trials of the given penetration
        STIM_onsets = STIM_file.STIM.photo_on;
     %only keep the selected trials onsets
        selected_STIM_onsets = cell2mat(STIM_onsets(selected_trials_idx.logicals(pnt).idx));
        selected_STIM_onsets = selected_STIM_onsets(1:4:end);
     %contact index/letter
        contact = STIM_file.STIM.chan;
        
     %2) get appropriate ns6 file of penetration 
     %isolate session date/animal for ns6 folder name
        underscore = strfind(penetration_name, '_');
        session =  penetration_name(1:underscore(2)-1);
    
        
      %load ns6 file of penetration of interest  
        %loop through ss files until we find an ss file with same
        %session and penetration as the penetration of interest
        for ss =1:length(selected_name_list)
            if contains(erase(selected_name_list(ss).name,'_ss.mat'),penetration_name(1:underscore(3)-1))
                ss_file = load(strcat(ssdir,selected_name_list(ss).name));
                ss_file_fieldnames = fieldnames(ss_file.ss);

            
                 cinterocdrft_names = ss_file_fieldnames(contains(ss_file_fieldnames,'cinterocdrft'));
                 for cint =1:length(cinterocdrft_names)
                     ns6_filename = char(cinterocdrft_names(cint)); 
                     ext          = 'ns6'; 
                     
                     file_dir = strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext);
                     %[t,wv] = load_open_ephys_data(file_dir);
                     %[t,wv] = LoadIntanSpikes(file_dir,1,4);
                      file_dir = strcat(nevdir,session,'\',ns6_filename(2:end),'.nev');
                    
                     
                     
                     clear ext NS_header banks neural 
                     % Read in NS Header
                     ext          = 'ns6'; 
                     NS_Header    = openNSx(strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext),'noread');
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
                            NS        = openNSx(strcat(ns6dir,session,'\',ns6_filename(2:end),'.',ext),electrode,'read','uV');
                            DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!


                            % preallocate data matrices 
                            if nct == 1
                                N.samples = length(DAT); 
                                MUA       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
                                sLFP       = zeros(ceil(N.samples/r2),N.neural);
                            end
                    %save each channel individually from each other to save space when loading a specific channel        
                   save(strcat(savedir,'\',ns6_filename, sprintf('_channel%d.continuous', e)), 'DAT');
                        end
                    end
                    
%{
                            % extract the aMUA.    
                            % low pass at 10000 Hz
                            clear lpc lWn bwb bwa 
                            lpc       = 10000;  % cutoff
                            lWn       = lpc/nyq;
                            [bwb,bwa] = butter(4,lWn,'low');
                            lpsLFP     = filtfilt(bwb,bwa,DAT); %low pass filter &rectify

                            % decimate analog MUA (aMUA) to get 1kHz samp freq
                           % MUA(:,nct) = decimate(lpMUA,r);
                            
                            % decimate sLFP to get 20kHz samp freq
                            sLFP(:,nct) = decimate(lpsLFP,r2); 
                            
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
                    % calculate CSD before triggering to trials OR on the trial data BUT not on
                    % the mean LFP. 

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
                    STIM_file.STIM.sLFP = trigData(Srt_sLFP,floor(selected_STIM_onsets./2),-pre2,post2); 
                    %STIM_aMUA = STIM_file.STIM.aMUA;
                    STIM_sLFP = STIM_file.STIM.sLFP;
                    STIM_sCSD = STIM_file.STIM.sCSD;


                    trigdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\trig_data';

                   % save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
                    save(strcat(trigdir,'\',ns6_filename, '_10khzlowpass_15khzdownsamp_sCSD.mat'), 'STIM_sCSD');
                    save(strcat(trigdir,'\',ns6_filename, '_10khzlowpass_15khzdownsamp_sLFP.mat'), 'STIM_sLFP');
                    contrast = '_domsup50_nondom0';
                    %LimePlot(STIM_aMUA, ns6_filename, ns6dir, contrast)
                  %  LimePlot(STIM_sCSD, ns6_filename, ns6dir, contrast)
                   % LimePlot(STIM_sLFP, ns6_filename, ns6dir, contrast)
                   % data = cat(4,STIM_sLFP,  STIM_sCSD, STIM_aMUA);
                    plotdir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\plots\';
                  %  newLimePlot(data, ns6_filename, plotdir, contrast)
                  %{
                    for n = 1:length(STIM_aMUA(1,1,:))
                   figure();
                   subplot(2,1,1)
                   x = -500:1500;
                   plot(x,STIM_aMUA(:,14,n), 'Color', 'b')
                   legend('aMUA')
                   xlim([-500 1500])
                   
                   subplot(2,1,2)
                   x= -7500:22500;
                   plot(x,STIM_sLFP(:,14,n), 'Color', 'r')
                   legend('s potential')
                   xlim([-7500 22500])
                    end
                    %}
                     %}
                 end
            end
        end
    end
end
 
chan_dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\indiv_chan_data\x160602_I_cinterocdrft012_channel14.continuous';
%%


%%
nev_dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\nev_selected_files\160602_I_cinterocdrft012.nev';
nev_dir = 'C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\nev_selected_files\140718_B_attnctrls003.nev';


%NEV = loadNEV(nev_dir, 'read');
[t,wv] = BlackrockNEVLoadingEngineU(nev_dir);

NEV2 = loadNEV(nev_dir,'read');

NEV = loadNEV(nev_dir,'read');


%%
 fcTT = {}; %Pass in a cell array of filenames
 TText = '.ntt';
 StartingDirectory = pwd;
 FDDirectory = 'FD';
 channelValidity = true(4,1);
 LoadingEngine = 'BlackrockNEVLoadingEngine';
 minClusters = 20;
 maxClusters = 60;
 maxSpikesBeforeSplit = []; % if isempty then don't split
 featuresToCompute = {'feature_Energy', 'feature_EnergyD1',...
'feature_Peak', 'feature_WavePC1', 'feature_Time'};
 featuresToUse = {'feature_Energy', 'feature_EnergyD1'};
 
 
RunClustBatch('minClusters', 10, 'maxClusters', 25);
RunClustBatch();
RunClustBatch('fcTT', {'TT1.ntt', 'TT2.ntt'});
 
 
 
 
 