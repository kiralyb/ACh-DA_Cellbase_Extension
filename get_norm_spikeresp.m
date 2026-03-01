function ach_b = get_norm_spikeresp(cellid1,newcues)
sigma = 0.02;
wind = [-0.0,0.4];
[~,spsth_a,~,tags,raszt_a] =  ultimate_psth(cellid1,'trial','StimulusOn',wind,'parts','#TrialCode','sigma',sigma);
for i = 1:length(tags)
    inxs(i) = str2double(tags{i}(end));
end
inx_4 = find(inxs==1);
inx_n = find(inxs==3);
[~,locs,width,prominance] = findpeaks(spsth_a(inx_4(1),:));
[~,ind] = max(prominance);
I_a_4Hz = locs(ind);
window_a_4Hz = ceil(width(ind));
[~,locs,width,prominance] = findpeaks(spsth_a(inx_n(1),:));
[~,ind] = max(prominance);
I_a_n = locs(ind);
window_a_n = ceil(width(ind));

RESP1_4Hz = sum(raszt_a{inx_4(1)}(:,max([1,I_a_n-window_a_n]):min([length(spsth_a),I_a_n+window_a_n])),2);
[psth_a,spsth_a,~,~,raszt_a] =  ultimate_psth(cellid1,'trial','StimulusOn',wind,'parts','All','sigma',sigma);
RESP1 = sum(raszt_a(newcues,max([1,I_a_n-window_a_n]):min([length(spsth_a),I_a_n+window_a_n])),2);

ach_b = RESP1/mean(RESP1_4Hz)';
