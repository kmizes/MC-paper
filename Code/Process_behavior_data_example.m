%% example script to process behavior data


%% Load behavioral data
load('\MC-task-interference\Behavior Data\Full task\J4-Rat126\ratBEHstruct.mat')

%% Pull example session
s = 1403; % CUE and WM
s = 1404; % OT

%% Calculate basic session metrics

if ratBEHstruct(s).protocol==7
    % Find index of CUE and WM trials
    blocklength = max(ratBEHstruct(s).blocknumRepair);
    idxCued = find(ratBEHstruct(s).blocknumRepair < (blocklength+1)/2 );
    idxWM = find(ratBEHstruct(s).blocknumRepair >= (blocklength+1)/2  );
    
    % calculate accuracy
    accCUE = mean(ratBEHstruct(s).Hit(idxCued));
    accWM  = mean(ratBEHstruct(s).Hit(idxWM));
    
    % calculate entropy (generally ran across sessions)
    sequse = 'RCR';
    idxCuedEntropy = intersect(idxCued, find(strcmp(ratBEHstruct(s).targetNames, sequse)));
    idxWMEntropy = intersect(idxWM, find(strcmp(ratBEHstruct(s).targetNames, sequse)));
    pokeNamesCued = ratBEHstruct(s).pokeNames(idxCuedEntropy);
    pokeNamesWM  = ratBEHstruct(s).pokeNames(idxWMEntropy);
    % - filter only pokes that are present >0.1% of the time
    a = unique([pokeNamesCued, pokeNamesWM],'stable');
    b = cellfun(@(x) sum(ismember([pokeNamesCued, pokeNamesWM],x)),a,'un',0);
    b = cell2mat(b) / sum(cell2mat(b));
    goodnames = a(b>.001);
    % - convert to probabilities
    pCued = cellfun(@(x) sum(ismember(pokeNamesCued,x)), goodnames,'un',0);
    pCued = cellfun(@(v) v+1, pCued, 'un',0); % laplace smoothing
    pCued = cell2mat(pCued) / sum(cell2mat(pCued));
    pWM = cellfun(@(x) sum(ismember(pokeNamesWM,x)), goodnames,'un',0);
    pWM = cellfun(@(v) v+1, pWM, 'un',0); 
    pWM = cell2mat(pWM) / sum(cell2mat(pWM));
    
    entropy_cued =  -sum(pCued .*log(pCued)); 
    entropy_wm = -sum(pWM .*log(pWM));
    
    % calculate trial times in seconds
    % - for successful trials, and all sequences
    idxCuedTrialTimes = intersect(idxCued,  find(cellfun(@length,ratBEHstruct(s).pokeTimes)==3 & ratBEHstruct(s).Hit));
    idxWMTrialTimes = intersect(idxWM,  find(cellfun(@length,ratBEHstruct(s).pokeTimes)==3 & ratBEHstruct(s).Hit));
    
    trialTimesCUE = median( cellfun(@(v) v(3)-v(1), ratBEHstruct(s).pokeTimes(idxCuedTrialTimes)) ) / 1000;
    trialTimesWM = median( cellfun(@(v) v(3)-v(1), ratBEHstruct(s).pokeTimes(idxWMTrialTimes)) ) / 1000;
    
    % calculate horizontal speeds for a sequence
    Ldist = 2.5; % cm distance of lever separation
    sequse = 'RCR'; % sequence to get speed
    duse12 = Ldist * (1+contains(sequse(1:2),{'LR','RL'}));
    duse23 = Ldist * (1+contains(sequse(2:3),{'LR','RL'}));
    idxCuedHorzSpeed = intersect(idxCued, find(strcmp(ratBEHstruct(s).pokeNames, sequse)));
    idxWMHorzSpeed = intersect(idxWM, find(strcmp(ratBEHstruct(s).pokeNames, sequse)));
    
    horzSpeedCUE = median( duse12 ./ cellfun(@(v) (v(2) - v(1))/1000, ratBEHstruct(s).pokeTimes(idxCuedHorzSpeed)) ...
        + duse23 ./ cellfun(@(v) (v(3) - v(2))/1000, ratBEHstruct(s).pokeTimes(idxCuedHorzSpeed)) );
    horzSpeedWM = median( duse12 ./ cellfun(@(v) (v(2) - v(1))/1000, ratBEHstruct(s).pokeTimes(idxWMHorzSpeed)) ...
        + duse23 ./ cellfun(@(v) (v(3) - v(2))/1000, ratBEHstruct(s).pokeTimes(idxWMHorzSpeed)) );
    
elseif ratBEHstruct(s).protocol==8 
    % calculate accuracy!
    accOT = mean(ratBEHstruct(s).Hit);
    % NOTE: success rate fell below 20% for 30 trials, or if rats failed to
    % press the lever once in the entire session, a cue is added back in
    % - below filters these trials
    % https://www.mathworks.com/matlabcentral/fileexchange/41813-runlength
    cuedNames = ratBEHstruct(s).cuedNames;
    if ~isempty(cuedNames)
    for trial = 1:length(cuedNames)
    [b,n] = RunLength_M(cuedNames{trial});
    cuedNames{trial} = b(n~=2); 
    end
    end
    % * end filter of cued names
    idxOT = find(cellfun(@length, cuedNames) <= 1);
    accOT = mean(ratBEHstruct(s).Hit(idxOT));
    
    
    % calculate entropy
    pokeNamesOT = ratBEHstruct(s).pokeNames(idxOT);
     % - filter only pokes that are present >0.1% of the time
    a = unique(pokeNamesOT,'stable');
    b = cellfun(@(x) sum(ismember(pokeNamesOT,x)),a,'un',0);
    b = cell2mat(b) / sum(cell2mat(b));
    goodnames = a(b>.001);
    % - convert to probabilities
    pOT = cellfun(@(x) sum(ismember(pokeNamesOT,x)), goodnames,'un',0);
    pOT = cellfun(@(v) v+1, pOT, 'un',0); % laplace smoothing
    pOT = cell2mat(pOT) / sum(cell2mat(pOT));
   
    entropy_OT =  -sum(pOT .*log(pOT)); 
    
    
    % calculate trial time in seconds
    idxOTTrialTimes = intersect(idxOT,  find(cellfun(@length,ratBEHstruct(s).pokeTimes)==3 & ratBEHstruct(s).Hit));
    trialTimesOT = median( cellfun(@(v) v(3)-v(1), ratBEHstruct(s).pokeTimes(idxOTTrialTimes)) ) / 1000;
    
    % calculate horizontal speed
    Ldist = 2.5; % cm distance of lever separation
    sequse = unique(ratBEHstruct(s).targetNames); sequse = sequse{1};
    duse12 = Ldist * (1+contains(sequse(1:2),{'LR','RL'}));
    duse23 = Ldist * (1+contains(sequse(2:3),{'LR','RL'}));
    
    horzSpeedOT = median( duse12 ./ cellfun(@(v) (v(2) - v(1))/1000, ratBEHstruct(s).pokeTimes(idxOTTrialTimes)) ...
        + duse23 ./ cellfun(@(v) (v(3) - v(2))/1000, ratBEHstruct(s).pokeTimes(idxOTTrialTimes)) );

    
end


