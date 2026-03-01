function firstNovelBehav_DREADD(AnimalGroup)
loadcb
animals = CELLIDLIST(AnimalGroup == 1);
controlanimals = CELLIDLIST(AnimalGroup == 0)
cellbasepath = getpref('cellbase','datapath');

D  = run_psychometric_batch(animals, 3, cellbasepath);
D1 = structfun(@(x) x(1,:), D, 'UniformOutput', false);
D2 = structfun(@(x) x(2,:), D, 'UniformOutput', false);
D3 = structfun(@(x) x(3,:), D, 'UniformOutput', false);
C  = run_psychometric_batch(controlanimals,1,cellbasepath);

plot_four_conditions(C.hit,  D1.hit,  D2.hit,  D3.hit,  ...
    'First new tone (7 Hz - go) hitrate')

plot_four_conditions(C.hitF, D1.hitF, D2.hitF, D3.hitF, ...
    '4 Hz - go hitrate')

plot_four_conditions(C.faF,  D1.faF,  D2.faF,  D3.faF,  ...
    '10 Hz - no-go false alarm rate')
end

function plot_four_conditions(control, d1, d2, d3, ylabel_text)

    data = {control, d1, d2, d3};
    x = 1:4;

    figure; hold on

    for i = 1:4
        scatter(ones(1,numel(data{i}))*i, data{i})
    end

    % stats
    sigstar([1 2], ranksum(data{1}, data{2}))
    sigstar([1 3], ranksum(data{1}, data{3}))
    sigstar([1 4], ranksum(data{1}, data{4}))
    sigstar([2 4], signrank(data{2}, data{4}))
    sigstar([2 3], signrank(data{2}, data{3}))
    sigstar([3 4], signrank(data{3}, data{4}))

    % paired lines
    line([2 3 4]' * ones(1,length(data{2})), ...
         [data{2}; data{3}; data{4}], ...
         'Color',[0.5 0.5 0.5 0.5])

    xticks(x)
    xticklabels({'Control\_C21','DREADD\_C21\_1','DREADD\_C21\_2','DREADD\_h2o'})
    xlim([0.5 4.5])
    ylabel(ylabel_text)

    setmyplot_balazs
end

function data = run_psychometric_batch(animals, nSessions, dataRoot)

    nAnimals = numel(animals);

    data.hit  = nan(nSessions, nAnimals);
    data.hitF = nan(nSessions, nAnimals);
    data.faF  = nan(nSessions, nAnimals);

    for a = 1:nAnimals
        animalPath = fullfile(dataRoot, animals{a});

        d = dir(animalPath);
        isSess = [d.isdir] & ~ismember({d.name},{'.','..'});
        sessions = sort({d(isSess).name});

        nUse = min(nSessions, numel(sessions));

        for s = 1:nUse
            [data.hit(s,a), data.hitF(s,a), data.faF(s,a)] = ...
                psychometric_7hzbehav_analysis(animals{a}, sessions{s});
            if s > 1 & isnan(data.hit(s-1,a))
                data.hit(s,a) = nan;
            end
        end
    end

    % return row vectors if nSessions == 1
    if nSessions == 1
        data = structfun(@(x) x(:)', data, 'UniformOutput', false);
    end
end
