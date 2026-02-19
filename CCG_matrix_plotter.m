function CCG_matrix_plotter(path,cluster1,cluster2)

figure
lags = -500:500;
[CCR_smooth,~,PairOfCells] = CCG_normalizer([path,'ISI1\CCG_matrices.mat']);
%CCR_smooth(500:1500);
%load([path,'ISI\CCG_matrices.mat']);
% CCR_norm2 = bsxfun(@rdivide,CCR,sum(CCR,2));
% CCR_Poi_norm = bsxfun(@rdivide,CCR_Poi,sum(CCR_Poi,2));
% CCR_norm = (CCR_norm2 - CCR_Poi_norm)*length(CCR);
% CCR_smooth = nan(size(CCR_norm));
% for ind = 1:size(CCR_norm,1)
%     if sum(CCR(ind,:),2) > 0
%     CCR_smooth(ind,:) = smoothdata(CCR_norm(ind,:),2,'gaussian',(length(CCR)/sqrt(SegmentLength(ind))));
%     else
%     'ize'    
%     end
% end

% Precompute positions of each cell in PairOfCells
pos1 = arrayfun(@(c) findcellpos(c), PairOfCells(:,1));
pos2 = arrayfun(@(c) findcellpos(c), PairOfCells(:,2));

%CCR_norm = CCR_norm(:,500:1500);

for x = 1:5
    % Precompute mask for this x
    idxX = (cluster1 == x);

    for y = 1:5
        % Precompute mask for this y
        idxY = (cluster2 == y);

        crosspairs = [];
        flipi = [];

        % Check all pairs at once (vectorized)
        a1 = idxX(pos1);
        a2 = idxX(pos2);
        d1 = idxY(pos1);
        d2 = idxY(pos2);

        mask = (a1 & d2) | (a2 & d1);
        crosspairs = find(mask);
        flipi = find(a2 & d1);

        subplot(5,5,(x-1)*5+y)

        % Use precomputed CCR_norm
        CCR_use = CCR_smooth;
        CCR_use(flipi,:) = flip(CCR_use(flipi,:),2);

        ccg_plot(CCR_use(crosspairs,:),0.8)

        % ---- Ticks & Labels ----
        setmyplot_balazs
        ax = gca;
        % Only bottom row gets x-ticks
        if x < 5
            ax.XTick = [];
        end
        % Only left column gets left y-ticks
        ax.YAxis(1).TickLabels = [];   
        % Rightmost column shows right y-ticks
        if y < 5
            ax.YAxis(2).TickLabels = [];
        end
    end
end


