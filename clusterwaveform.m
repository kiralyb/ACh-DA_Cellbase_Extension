load('D:\CellBase_DATA\Cellbase_ACh_DA')

figure
for j=1:5
    i_neur=find(abs([TheMatrix{:,25}])==j);
    %i_neur=inhibitory_inx;
    clear MEAN_WAVES 
    k = 1;
    for i = 1:length(i_neur)
        try
        cellid = CELLIDLIST{i_neur(i)}
        tsegs_spont = rel2abstimes(cellid,[-2,0],'stim','BurstOn');   % extract 2s periods before bursts
        selts_spont = extractSegSpikes(cellid,tsegs_spont);     % extract spontaneous spikes
        wave_spont = extractSpikeWaveforms(cellid,selts_spont,'chans','all');    % get waveforms for the extracted spikes
        
        % Downsample
        nm_spont = size(wave_spont,1);
        rp = randperm(min(2000,nm_spont));
        wsds = wave_spont(rp,:,:);
        
        % Average waveforms
        mean_spont = squeeze(nanmedian(wave_spont,1));
        [NORM,IND]=max(max(mean_spont')-min(mean_spont'));
        %[NORM,IND]=max(max(mean_spont'));
        MEAN_WAVES(i,:) = mean_spont(IND,:)/NORM;
        k = k+1;
        end
    end
    subplot(1,5,j)
    plot(MEAN_WAVES','Color',[0.75,0.75,0.75,0.15])
    hold on
    plot(mean(MEAN_WAVES)','Color',[0,0,1],'LineWidth',2)    
end