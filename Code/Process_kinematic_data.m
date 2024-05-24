%% Example code for processing kinematic data for "The role of motor cortex in motor sequence execution depends on demands for flexibility"

%% Loop through rats and calculate trajectory similarities
% also speed and smoothness

corrDist_rat = {};
corrDist_across = {};
corrDist_mean = {}; % correlation to mean movement prelesion

mean_speed_all = {}; % condition x context x rat
mean_speed_all_earlylate = {}; % first 1/3 and last 1/3? or no, too many for OT?
mean_speed_all_earlylate;
peak_speed_all = {};
peak_speed_all_earlylate = {};

LDJ_all = []; speed_metric_all = []; SARCL_all = [];
nsampsave = [];
traj_speed_all_earlylate = {}; traj_speed_all = {};


dotrajplots = 1;
toss_lever_taps = 1;
doallmove = 0;
dolate = 1; 
dolocal = 1; % local warping for correlation or global warping

% tap IPI use
curContext = {'OT','Cued','WM'};

for condition = 2:3 % loop through task context
for rat = 1:length(fpath_rat) % loop through rat
    filepath = fpath_rat{rat};
    if doallmove
         load(fullfile(filepath, ['lesion_data_allmoves_' curContext{condition} '.mat']))
    else
        load(fullfile(filepath, ['lesion_data_' curContext{condition} '.mat']))
    end
    disp(filepath);
    assert(length(seqchange)==3);
    
%% skip pre mock break, or placeholder after learning

id = seqchange(1);
if dolate % resamples seq changes so get late
    seqchange = [seqchange(2), seqchange(3), seqchange(3) + floor((length(traj) - seqchange(3))/2)];
end
% now toss all the early data?
sessID(1:id) = [];
tSess(1:id) = [];
tapon(:,1:id) = [];
traj(1:id) = [];
traj_prob(1:id) = [];

seqchange = seqchange - id;

%% 0) prep data

all_titles = {'Pre','Uni','Bi'};

% smooth trajectories + fill nans + remove bad trajs
badtrackfrac = []; trajvar = [];
for count = 1:length(traj)
    % smooth and filter...
    traj{count}(traj_prob{count}<.9 , :) = nan;
    badtrackfrac(count) = sum(traj_prob{count}<.9)/length(traj_prob{count});
    traj{count} = fillmissing(traj{count}, 'linear','EndValues','nearest');
    traj{count} = imgaussfilt(traj{count}, [.6, eps]);
    trajvar(:,count) = var(traj{count});
end
badid = find(badtrackfrac > .1); %
badid = unique([badid, find(trajvar(1,:)<10 & trajvar(2,:)<10)]);
% if the third tap is too close to the end of trial, cut!
a = cellfun(@length, traj); b = tapon(end,:)+8;
badid = unique([badid, find(b>a)]);
badid = unique([badid, find(tapon(2,:) > tapon(end,:))]);

% remove bad stuff
for j = 1:length(seqchange)
    seqchange(j) = seqchange(j) - sum(badid<seqchange(j));
end
sessID(badid) = [];
tSess(badid) = [];
tapon(:,badid) = [];
traj(badid) = [];
traj_prob(badid) = [];

%% structure data samples by context and trial time
% - remove abnormally long/short trials
% pull data
seqedge = [0 seqchange, length(traj)];
numsamp = min(diff(seqedge))-1; % pull from each seqchange

% put data into structure by context
% - excluding extremem trial times
trajtemp = {};
for dim_use = 1:2
for j = 1:(length(seqedge)-1)
    rr = (seqedge(j)+1) : (seqedge(j+1));

    tap_range = round(quantile(tapon(end,rr),[.05,.5,.75]));
    tap_range = round(quantile(tapon(end,rr),[0,.5,.9])); % want to toss slowest trials?
    
    % 2) crop from similar tap times
    [ttt,idx] = sort(tapon(end,rr));
    [~,id1] = min(abs(ttt-tap_range(1))); id1 = id1(1);
    [~,id2]=min(abs(ttt-tap_range(end))); id2 = id2(end);
    
    trajsamp = traj(rr(idx(id1:round((id2-id1)/numsamp):id2)));
    tapstemp{j} = tapon(:,rr(idx(id1:round((id2-id1)/numsamp):id2)));
    
    sessIDtemp{j} = sessID(rr(idx(id1:round((id2-id1)/numsamp):id2)));        
    
    % take only 1 dim?
    trajsamp = cellfun(@(v) v(:,dim_use),trajsamp,'un',0);
    [~,sortid] = sort(cellfun(@length,trajsamp));
    % zero pad to same length...
    maxN = max(cellfun(@length,trajsamp));
    padfun = @(v) [v; zeros( maxN - numel(v),1)] ;
    trajsamp = cellfun(padfun, trajsamp , 'un', 0);
    % plot?
    trajsamp = cell2mat(trajsamp);
    trajsamp(trajsamp==0)=nan;
   
    trajtemp{j,dim_use} = trajsamp';

end
end
if length(trajtemp)==4
    trajtemp_mock = trajtemp; trajtemp(1,:) = []; 
    tapstemp_mock = tapstemp; tapstemp(1) = []; 
end


%% pull out speed between lever taps
% This was not used in the MC paper analyses I beliave

traj_speed = {};
mean_speed = []; peak_speed = [];
for count = 1:length(traj)
    % smooth and oversample
    val = traj{count};
    valx = csaps(1:size(val,1), val(:,1), .5, linspace(1,size(val,1),4*size(val,1)));
    valy = csaps(1:size(val,1), val(:,2), .5, linspace(1,size(val,1),4*size(val,1)));
    
    dt = 1/(40*4);
    traj_speed{count} = sqrt(diff(valx).^2 + diff(valy).^2)./dt; % 40 hz * 4 upsample to get pixels/second 
    % - i think the scale factor is wrong...
    
    % - need to account for shift in tap time with upsampling...
    % - or just do for each tap
    [~,newtaps] = min( abs( tapon(:,count) - linspace(1,size(val,1),4*size(val,1)) )' );
    
    % dumbest, just look at change in average velocities pre post?
    mean_speed(count) = mean(traj_speed{count}(newtaps(1):newtaps(end)));
    % 
    peak_speed(count) = max(traj_speed{count}(newtaps(1):newtaps(end)));
    
end


mean_speed_all; % condition x context x rat
for j = 1:(length(seqedge)-1)
    rr = (seqedge(j)+1) : (seqedge(j+1));
    mean_speed_all{j, condition, rat} = mean_speed(rr);
    peak_speed_all{j, condition, rat} = peak_speed(rr);
end

mean_speed_all_earlylate;
nsamp = min([300, round(diff(seqedge)/3)]);
nsampsave(condition,rat) = nsamp;
for j = 1:(length(seqedge)-1)
    rr = (seqedge(j)+1) : (seqedge(j+1));
    jj = 2*(j-1)+1;
    mean_speed_all_earlylate{jj, condition, rat} = mean_speed(rr(1:nsamp));
    mean_speed_all_earlylate{jj+1, condition, rat} = mean_speed(rr( (end-nsamp+1) : end));
    peak_speed_all_earlylate{jj, condition, rat} = peak_speed(rr(1:nsamp));
    peak_speed_all_earlylate{jj+1, condition, rat} = peak_speed(rr( (end-nsamp+1) : end));
end

frac_stop_all_earlylate;
for j = 1:(length(seqedge)-1)
    rr = (seqedge(j)+1) : (seqedge(j+1));
    jj = 2*(j-1)+1;
    frac_stop_all_earlylate{jj, condition, rat} = frac_stop(rr(1:nsamp));
    frac_stop_all_earlylate{jj+1, condition, rat} = frac_stop(rr( (end-nsamp+1) : end));

    traj_speed_all{condition,rat} = [traj_speed{:}];
    traj_speed_all_earlylate{jj,condition,rat} = [traj_speed{rr(1:nsamp)}]; % 1,3,5,7 is early (early, mock, early uni, early bi)
    traj_speed_all_earlylate{jj+1,condition,rat} = [traj_speed{rr((end-nsamp+1):end)}]; % 2,4,6,8 is premock, preuni, prebi, late)
end
% measure of smoothness
% - look at speed metric
% - look at log dimensionless jerk
speed_metric = []; LDJ = []; SAL = [];
for count = 1:length(traj)
    % speed metric
    speed_metric(count) = peak_speed(count) / mean_speed(count);
    % LDJ
    jerk = diff(diff(traj_speed{count}*100));
    t = (1:length(traj_speed{count}))/((40*4)); t = t(2:end-1); % diff 
    LDJ(count) = -log( trapz(t,jerk.^2) * (t(end)-t(1))^3/(peak_speed(count)^2));
    % spectral arc length
    SAL(count) = SpectralArcLength(traj_speed{count}', 40*4);
end
for j = 1:(length(seqedge)-1)
    rr = (seqedge(j)+1) : (seqedge(j+1));
    jj = 2*(j-1)+1;
    speed_metric_all(jj,condition,rat) = mean(speed_metric(rr(1:nsamp)));
    speed_metric_all(jj+1,condition,rat) = mean(speed_metric(rr((end-nsamp+1):end)));
    LDJ_all(jj,condition,rat) = mean(LDJ(rr(1:nsamp)));
    LDJ_all(jj+1,condition,rat) = mean(LDJ(rr((end-nsamp+1):end)));
    SARCL_all(jj,condition,rat) = mean(SAL(rr(1:nsamp)));
    SARCL_all(jj+1,condition,rat) = mean(SAL(rr((end-nsamp+1):end)));
end
SAL_save = SAL;

%% plot traj examples

if dotrajplots

% figure out IPI to use?
% - across cued and WM?
% - need to do outside of this
    
N = 8;
c = linspace(.2,.8,N);
figure; ax1=[];
cbehuse = [.1,1,.1; 1,.1,.1; .1,.1,1];
cbehuse = [.3,.9157,.7784; 1,.1,.1; .1,.1,1];
for dim_use = 1:2
for j = 1:length(trajtemp)
    ax1(j) = subplot(1,size(trajtemp,1),j); hold on;
    taxis = ((1:size(trajtemp{j,dim_use},2))-25)*1/40;
    k = round(size(trajtemp{j,dim_use},1)*1/3);
    for i = 1:N
        k = round(size(trajtemp{j,dim_use},1)*1/2);
        val = trajtemp{j,dim_use}(k+i,:);
        if dim_use==2; val = val + 300; end
        plot(taxis(1:length(val)),val, 'Color',cbehuse(condition,:)*c(i),'linewidth',2);
    end
    xlim([taxis(1),taxis(sum(~isnan(val)))])
    xlim([-.6, 3]);
    xlabel('Time (seconds)'); ylabel('Position')
    title(all_titles{j});%ylim([-3,3])
    set(gca, 'YDir','reverse')

end
end
linkaxes(ax1,'xy');

end

%% 2) traj correlation coefficient
%  - save to do over rats?

% interpolate
trajtemp = trajtemp_mock;
tapstemp = tapstemp_mock;

% plot mean correlation
% warp 
interplen = round(median(tapon(3,:) - tapon(1,:))) + 2*round(median(tapon(1,:)));

uselen = round(median(diff(tapon)'))*2;
endlen = min(cellfun(@length, traj) - tapon(3,:));
endlen = 10;

% template tap times to warp to?
tap_interplen = [];
trajsamp = {};
for dim_use = 1:2
for t = 1:length(trajtemp)
    shuffid = randperm(size(trajtemp{t},1));
    %count = 1;
    for trial = 1:size(trajtemp{t},1)
        a = trajtemp{t,dim_use}(trial,:);
        a(isnan(a)) = [];
       
        % local linear warp
        if dolocal %contains(filepaths{f}, 'L3-Rat34') % what is this
        tt = tapstemp{t}(:,trial);
        % temp really bad fix...
        if length(a)-tt(3) < endlen
            tt(3) = length(a)-endlen-1;
        end
            
        % tt = tapstemp{t}(:,shuffid(trial)); % shuffled version?
        if (tt(3)-tt(2) < 3); tt(3) = tt(3) + 3; continue; end
        if (tt(2)-tt(1) < 3); continue; end
        b = a(1:tt(1));
        b = [b, interp1(linspace(1,uselen(1),tt(2)-tt(1)), a(tt(1):(tt(2)-1)), 1:uselen(1))];
        b = [b, interp1(linspace(1,uselen(2),tt(3)-tt(2)), a(tt(2):(tt(3)-1)), 1:uselen(2))];
        b = [b, a(tt(3): (tt(3)+endlen))];
        a = b;
        
        if toss_lever_taps
            a = trajtemp{t,dim_use}(trial,:);
            a(isnan(a)) = [];
            b = [interp1(linspace(1,uselen(1),tt(2)-tt(1)), a(tt(1):(tt(2)-1)), 1:uselen(1)),...
                interp1(linspace(1,uselen(2),tt(3)-tt(2)), a(tt(2):(tt(3)-1)), 1:uselen(2))];
            b = b(10:end-10); % 10 frames toss
            a = b;
        end
        
        else
            % global warp
            a = interp1(linspace(1,interplen,length(a)),a,1:interplen);
        
        end
        
        trajsamp{t,dim_use}(:,trial) = zeros(size(a));
        
        if any(isnan(a)); continue; end;
        if var(a)<1; continue; end
        
        trajsamp{t,dim_use}(:,trial) = a;
        %count = count+1;
    end
    trajsamp{t,dim_use} = trajsamp{t,dim_use}';
end % seqchange segment
end % dim use



% calculate correlations
% - need to do within? or across?
% lets try just within first
corrDist = {};
for j = 1:length(trajsamp)
% mean normalize
a1 = trajsamp{j,1}'; a2 = trajsamp{j,2}';
badid = find(sum(a1)<=10 | sum(a2)<=10); % ah if all zeros...
a1(:,badid) = []; a2(:,badid) = [];
a1 = (a1 - mean(a1))./std(a1);
a2 = (a2 - mean(a2))./std(a2);
C = corrcoef([a1;a2]);
corrDist{j} = squareform(C - eye(size(C)),'tovector');
end

% save 
corrDist_rat(:,rat,condition) = corrDist;

% ok time to calculate correlations across...
%corrDist_across = {};
nedge = cellfun(@(v) size(v,1), trajsamp(:,1)); nedge = [0; cumsum(nedge)];
% mean normalize
a1 = vertcat(trajsamp{:,1})'; a2 = vertcat(trajsamp{:,2})';
badid = find(sum(a1)<=10 | sum(a2)<=10); 
a1(:,badid) = []; a2(:,badid) = [];
for jjj=length(badid):-1:1; nedge(nedge>=badid(jjj)) = nedge(nedge>=badid(jjj))-1; end
a1 = (a1 - mean(a1))./std(a1);
a2 = (a2 - mean(a2))./std(a2);
C = corrcoef([a1;a2]);
% which comparison to make?
% - say pre-lesion to post-bilateral lesion
c1 = C((nedge(2)+1):(nedge(2+1)),(nedge(4)+1):(nedge(4+1))); c1 = c1(:);
corrDist_across{rat,condition} = c1;

% do late preuni, prebi, earlypost, late post
% 
if dolate
    c1 = C((nedge(1)+1) : (nedge(1+1)), (nedge(3)+1) : (nedge(3+1))); c1 = c1(:);
    corrDist_across{rat,condition} = c1; 
end


% ok correlations to means
% - currently its to the pre-lesion mean?
a_mean = mean([a1(:,(nedge(2)+1):(nedge(2+1))); a2(:,(nedge(2)+1):(nedge(2+1)))],2);
C = corrcoef([a_mean, [a1;a2]]);
for jjj = 1:length(trajsamp)
    corrDist_mean{jjj,rat,condition} = C(1,(nedge(jjj)+1):(nedge(jjj+1)));
end

%% subsample traj samp for OT comparisons
p = cellfun(@(v) size(v,1), trajsamp);
n = min([200; p(:)]);

for i = 1:size(trajsamp,1)
    for j = 1:size(trajsamp,2)
        nuse = randsample(p(i,j), n);
        trajsamp{i,j} = trajsamp{i,j}(nuse,:);
    end
end

end % over rats
end % over context

%% plot

% comparing  across lesion conditions
val = cellfun(@median, corrDist_rat); 
val2 = cellfun(@median, corrDist_across); 
val(3,:,:) = val(4,:,:); % replace unilateral with bilateral
val(4,:,:) = val2; % replace bilateral with across


figure; hold on;

bar(1, nanmean(val(2,:,2),2),'FaceColor',[1,.2,.2]);
bar(2, nanmean(val(3,:,2),2),'FaceColor',[1,.4,.4]);
bar(3, nanmean(val(4,:,2),2),'FaceColor',[1,.7,.7]);

bar(5, nanmean(val(2,:,3),2),'FaceColor',[.2,.2,1]);
bar(6, nanmean(val(3,:,3),2),'FaceColor',[.4,.4,1]);
bar(7, nanmean(val(4,:,3),2),'FaceColor',[.7,.7,1]);

bar(9, nanmean(val(2,:,1),2),'FaceColor',[.2,1,.2]);
bar(10, nanmean(val(3,:,1),2),'FaceColor',[.4,1,.4]);
bar(11, nanmean(val(4,:,1),2),'FaceColor',[.7,1,.7]);

% plot individual rats
for rat = 1:(length(fpath_rat) - length(tossrats))
    if contains(fpath_rat{rat}, {'159','151','157','152'})
        plot(1:3, val(2:4,rat,2),'Color',[1,0,0,.4]);
        plot(5:7, val(2:4,rat,3),'Color',[1,0,0,.4]);
    else
    plot(1:3, val(2:4,rat,2),'Color',[0,0,0,.4]);
    plot(5:7, val(2:4,rat,3),'Color',[0,0,0,.4]);
    end
    plot(9:11, val(2:4,rat,1),'Color',[0,0,0,.4]);
end


[~,ppreunicued] = ttest(val(2,:,2), val(3,:,2));
[~,pprebicued] =  ttest(val(2,:,2), val(4,:,2));
[~,ppreuniwm] = ttest(val(2,:,3), val(3,:,3));
[~,pprebiwm] =  ttest(val(2,:,3), val(4,:,3));
[~,ppreuniot] = ttest(val(2,:,1), val(3,:,1));
[~,pprebiot] =  ttest(val(2,:,1), val(4,:,1));

sigstar({[1,2],[1,3],[5,6],[5,7],[9,10],[9,11]},...
    [ppreunicued,pprebicued,ppreuniwm,pprebiwm,ppreuniot,pprebiot],0,1);
% only ot
% sigstar({[9,10],[9,11]},...
%     [ppreuniot,pprebiot],0,1);

% im not testing the thign i want to be testing
[~,p23] = ttest(val(3,:,2), val(4,:,2));
[~,p67] = ttest(val(3,:,3), val(4,:,3));
[~,p1011] = ttest(val(3,:,1), val(4,:,1));
% 
[p23] = signrank(val(3,:,2), val(4,:,2));
[p67] = signrank(val(3,:,3), val(4,:,3));
[p1011] = signrank(val(3,:,1), val(4,:,1));
sigstar({[2,3],[6,7],[10,11]},[p23,p67,p1011],0,1);

xticks([2,6,10]); %ylim([0,1])
xticklabels({'Cued','WM','OT'});
