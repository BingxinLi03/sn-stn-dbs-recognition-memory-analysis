%% firstlevel_oldnew_reviewer.m
% First-level task-fMRI GLM for valid old/new retrieval trials.
% OldValid = Hit + Miss.
% NewValid = FalseAlarm + CorrectRejection.
% Invalid retrieval, encoding, and interference events are modeled as nuisance conditions.
% Six realignment parameters are included as nuisance regressors.

clear; clc;

SPM_DIR = prompt_existing_dir('Enter the full path to the SPM12 directory: ');
RAW_ROOT = prompt_existing_dir('Enter the full path to the first-level RawData directory: ');
RESULT_ROOT = prompt_output_dir('Enter the full path for first-level output: ');
GM_MASK_PATH = prompt_existing_file('Enter the full path to the group-level mask NIfTI file: ');

addpath(SPM_DIR);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

TR = 2;
HPF = 128;
STIM_LIST = {'stim0', 'stim10', 'stim130'};
EXCLUDE_SUBJECTS = {'sub36'};
EXCLUDE_SESSIONS = cell(0, 2);
SKIP_IF_SPM_EXISTS = true;
OVERWRITE_INCOMPLETE_FOLDER = false;

QC.validRetrievalMin = 45;
QC.validOldMin = 12;
QC.validNewMin = 12;

ANALYSIS_NAME = 'Simple_Old_vs_New';
ANALYSIS_ROOT = fullfile(RESULT_ROOT, ANALYSIS_NAME);
if exist(ANALYSIS_ROOT, 'dir') ~= 7
    mkdir(ANALYSIS_ROOT);
end

LOG_FILE = fullfile(RESULT_ROOT, 'firstlevel_oldnew_QC_log.xlsx');
logRows = {};
logVars = {
    'analysis', 'stim', 'subject', 'status', 'reason', ...
    'n_scans', 'run_length_sec', 'max_included_onset_sec', ...
    'n_Hit', 'n_Miss', 'n_FalseAlarm', 'n_CorrectRejection', ...
    'n_InvalidOld', 'n_InvalidNew', ...
    'n_valid_retrieval', 'n_valid_old', 'n_valid_new', ...
    'qc_pass', 'output_dir'
};

for s = 1:numel(STIM_LIST)
    stimName = STIM_LIST{s};
    stimRawDir = fullfile(RAW_ROOT, stimName);
    stimOutRoot = fullfile(ANALYSIS_ROOT, stimName);

    if exist(stimRawDir, 'dir') ~= 7
        warning('Raw stimulation directory not found: %s', stimRawDir);
        continue;
    end
    if exist(stimOutRoot, 'dir') ~= 7
        mkdir(stimOutRoot);
    end

    subjectDirs = dir(stimRawDir);
    subjectDirs = subjectDirs([subjectDirs.isdir]);
    subjectDirs = subjectDirs(~ismember({subjectDirs.name}, {'.', '..'}));

    for i = 1:numel(subjectDirs)
        rawSubjectName = subjectDirs(i).name;
        subject = char(normalize_subject_id(rawSubjectName));

        if ismember(string(subject), normalize_subject_list(EXCLUDE_SUBJECTS))
            continue;
        end
        if is_excluded_session(subject, stimName, EXCLUDE_SESSIONS)
            continue;
        end

        funcDir = fullfile(stimRawDir, rawSubjectName);
        outDir = fullfile(stimOutRoot, subject);
        spmmat = fullfile(outDir, 'SPM.mat');

        fprintf('\nProcessing %s | %s\n', stimName, subject);

        try
            if exist(spmmat, 'file') == 2 && SKIP_IF_SPM_EXISTS
                logRows = add_log(logRows, ANALYSIS_NAME, stimName, subject, ...
                    'skipped_existing', 'SPM.mat already exists', ...
                    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, true, outDir);
                continue;
            end

            scans = cellstr(spm_select('ExtFPList', funcDir, '^swr.*\.nii$', Inf));
            scans = scans(~cellfun(@(x) isempty(strtrim(x)), scans));
            if isempty(scans)
                error('No swr*.nii functional scans found.');
            end

            nScans = numel(scans);
            runLengthSec = nScans * TR;

            rpFile = find_latest_rp_file(funcDir);
            if isempty(rpFile)
                error('No rp*.txt motion-parameter file found.');
            end

            onsetFile = find_onset_file(funcDir, subject);
            if isempty(onsetFile)
                error('No onset_%s.xlsx file found.', subject);
            end

            events = read_onset_events(onsetFile);
            counts = get_event_counts(events);
            qcResult = evaluate_qc(counts, QC);
            condDefs = make_condition_defs(events);

            if isempty(condDefs)
                qcResult.pass = false;
                qcResult.reason = append_reason(qcResult.reason, 'No model conditions were created.');
            end

            maxIncludedOnset = get_max_onset(condDefs);
            if ~isnan(maxIncludedOnset) && maxIncludedOnset > runLengthSec + 5
                qcResult.reason = append_reason(qcResult.reason, ...
                    sprintf('Maximum modeled onset %.3f s exceeds nominal run length %.3f s.', ...
                    maxIncludedOnset, runLengthSec));
            end

            if ~qcResult.pass
                logRows = add_log(logRows, ANALYSIS_NAME, stimName, subject, ...
                    'skipped_QC', qcResult.reason, ...
                    nScans, runLengthSec, maxIncludedOnset, ...
                    counts.nHit, counts.nMiss, counts.nFA, counts.nCR, ...
                    counts.nInvalidOld, counts.nInvalidNew, ...
                    counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, ...
                    false, outDir);
                continue;
            end

            if exist(outDir, 'dir') == 7
                if OVERWRITE_INCOMPLETE_FOLDER && exist(spmmat, 'file') ~= 2
                    rmdir(outDir, 's');
                    mkdir(outDir);
                end
            else
                mkdir(outDir);
            end

            finalMask = make_final_mask(funcDir, outDir, scans, subject, GM_MASK_PATH);
            run_firstlevel_spm(outDir, scans, rpFile, finalMask, condDefs, TR, HPF);
            create_contrasts(fullfile(outDir, 'SPM.mat'));

            logRows = add_log(logRows, ANALYSIS_NAME, stimName, subject, ...
                'success', qcResult.reason, ...
                nScans, runLengthSec, maxIncludedOnset, ...
                counts.nHit, counts.nMiss, counts.nFA, counts.nCR, ...
                counts.nInvalidOld, counts.nInvalidNew, ...
                counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, ...
                true, outDir);

        catch ME
            logRows = add_log(logRows, ANALYSIS_NAME, stimName, subject, ...
                'failed', ME.message, ...
                NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                NaN, NaN, NaN, false, outDir);
            fprintf('[FAILED] %s | %s | %s\n', stimName, subject, ME.message);
        end
    end
end

logTable = cell2table_or_empty(logRows, logVars);
if exist(LOG_FILE, 'file') == 2
    delete(LOG_FILE);
end
writetable(logTable, LOG_FILE, 'Sheet', 'session_log');

notes = {
    'First-level model', 'OldValid = Hit + Miss; NewValid = FalseAlarm + CorrectRejection';
    'Nuisance conditions', 'InvalidRetrieval, Encoding, Interference';
    'Motion regressors', 'Six parameters from rp*.txt';
    'High-pass filter', sprintf('%d s', HPF);
    'Contrast 1', 'OldValid > NewValid';
    'Contrast 2', 'NewValid > OldValid';
    'Contrast 3', 'AllValidRetrieval > baseline';
};
writetable(cell2table(notes, 'VariableNames', {'item', 'value'}), LOG_FILE, 'Sheet', 'notes');

fprintf('\nFirst-level analysis completed.\nLog: %s\n', LOG_FILE);

function pathOut = prompt_existing_dir(message)
    pathOut = clean_path_input(input(message, 's'));
    assert(exist(pathOut, 'dir') == 7, 'Directory not found: %s', pathOut);
end

function pathOut = prompt_output_dir(message)
    pathOut = clean_path_input(input(message, 's'));
    if exist(pathOut, 'dir') ~= 7
        mkdir(pathOut);
    end
end

function pathOut = prompt_existing_file(message)
    pathOut = clean_path_input(input(message, 's'));
    assert(exist(pathOut, 'file') == 2, 'File not found: %s', pathOut);
end

function pathOut = clean_path_input(pathIn)
    pathOut = strtrim(pathIn);
    if numel(pathOut) >= 2
        if (pathOut(1) == '"' && pathOut(end) == '"') || ...
           (pathOut(1) == '''' && pathOut(end) == '''')
            pathOut = pathOut(2:end-1);
        end
    end
end

function normList = normalize_subject_list(listIn)
    normList = strings(numel(listIn), 1);
    for i = 1:numel(listIn)
        normList(i) = normalize_subject_id(listIn{i});
    end
end

function tf = is_excluded_session(subject, stim, exSessions)
    tf = false;
    for i = 1:size(exSessions, 1)
        if string(subject) == normalize_subject_id(exSessions{i, 1}) && ...
                string(stim) == string(exSessions{i, 2})
            tf = true;
            return;
        end
    end
end

function subject = normalize_subject_id(value)
    value = strtrim(lower(char(string(value))));
    token = regexp(value, '\d+', 'match', 'once');
    if isempty(token)
        subject = string(value);
    else
        subject = string(sprintf('sub%02d', str2double(token)));
    end
end

function rpFile = find_latest_rp_file(funcDir)
    rpFile = '';
    files = dir(fullfile(funcDir, 'rp*.txt'));
    if isempty(files)
        return;
    end
    [~, idx] = max([files.datenum]);
    rpFile = fullfile(files(idx).folder, files(idx).name);
end

function onsetFile = find_onset_file(funcDir, subject)
    onsetFile = fullfile(funcDir, sprintf('onset_%s.xlsx', subject));
    if exist(onsetFile, 'file') == 2
        return;
    end
    candidates = dir(fullfile(funcDir, 'onset_sub*.xlsx'));
    if isempty(candidates)
        onsetFile = '';
    else
        onsetFile = fullfile(candidates(1).folder, candidates(1).name);
    end
end

function events = read_onset_events(onsetFile)
    events = readtable(onsetFile, 'Sheet', 'all_glm_events', 'VariableNamingRule', 'preserve');
    required = {'trial_type', 'onset'};
    for i = 1:numel(required)
        assert(ismember(required{i}, events.Properties.VariableNames), ...
            'Required column not found: %s', required{i});
    end
    if ~ismember('duration', events.Properties.VariableNames)
        events.duration = zeros(height(events), 1);
    end
    events.trial_type = string(events.trial_type);
    events.onset = double(events.onset);
    events.duration = double(events.duration);
    events = events(~isnan(events.onset), :);
end

function counts = get_event_counts(events)
    trialType = string(events.trial_type);
    counts.nHit = sum(trialType == "Hit");
    counts.nMiss = sum(trialType == "Miss");
    counts.nFA = sum(trialType == "FalseAlarm");
    counts.nCR = sum(trialType == "CorrectRejection");
    counts.nInvalidOld = sum(trialType == "InvalidOld");
    counts.nInvalidNew = sum(trialType == "InvalidNew");
    counts.nValidOld = counts.nHit + counts.nMiss;
    counts.nValidNew = counts.nFA + counts.nCR;
    counts.nValidRetrieval = counts.nValidOld + counts.nValidNew;
end

function result = evaluate_qc(counts, qc)
    reasons = {};
    if counts.nValidRetrieval < qc.validRetrievalMin
        reasons{end+1} = sprintf('Valid retrieval trials < %d.', qc.validRetrievalMin); %#ok<AGROW>
    end
    if counts.nValidOld < qc.validOldMin
        reasons{end+1} = sprintf('Valid old trials < %d.', qc.validOldMin); %#ok<AGROW>
    end
    if counts.nValidNew < qc.validNewMin
        reasons{end+1} = sprintf('Valid new trials < %d.', qc.validNewMin); %#ok<AGROW>
    end
    result.pass = isempty(reasons);
    if result.pass
        result.reason = 'pass';
    else
        result.reason = strjoin(reasons, ' ');
    end
end

function condDefs = make_condition_defs(events)
    condDefs = struct('name', {}, 'onset', {}, 'duration', {});
    condDefs = add_condition(condDefs, events, 'OldValid', ["Hit", "Miss"]);
    condDefs = add_condition(condDefs, events, 'NewValid', ["FalseAlarm", "CorrectRejection"]);
    condDefs = add_condition(condDefs, events, 'InvalidRetrieval', ["InvalidOld", "InvalidNew"]);
    condDefs = add_condition(condDefs, events, 'Encoding', "Encoding");
    condDefs = add_condition(condDefs, events, 'Interference', "Interference");
end

function condDefs = add_condition(condDefs, events, name, sourceTypes)
    idx = ismember(string(events.trial_type), string(sourceTypes));
    if ~any(idx)
        return;
    end
    subset = sortrows(events(idx, :), 'onset');
    condition.name = name;
    condition.onset = double(subset.onset(:));
    condition.duration = double(subset.duration(:));
    condDefs(end+1) = condition;
end

function value = get_max_onset(condDefs)
    value = NaN;
    onsets = [];
    for i = 1:numel(condDefs)
        onsets = [onsets; condDefs(i).onset(:)]; %#ok<AGROW>
    end
    if ~isempty(onsets)
        value = max(onsets);
    end
end

function finalMask = make_final_mask(funcDir, outDir, scans, subject, gmMaskPath)
    maskCandidates = dir(fullfile(funcDir, sprintf('glm_mask_%s_*.nii', subject)));
    if isempty(maskCandidates)
        maskCandidates = dir(fullfile(funcDir, 'glm_mask_*.nii'));
    end
    if isempty(maskCandidates)
        finalMask = gmMaskPath;
        return;
    end

    artifactMask = fullfile(maskCandidates(1).folder, maskCandidates(1).name);
    gmCopy = fullfile(outDir, 'GM_mask.nii');
    artifactCopy = fullfile(outDir, 'Artifact_mask.nii');
    copyfile(gmMaskPath, gmCopy, 'f');
    copyfile(artifactMask, artifactCopy, 'f');

    inputImages = char(scans{1}, [gmCopy ',1'], [artifactCopy ',1']);
    flags = struct('mask', 0, 'mean', 0, 'interp', 0, 'which', 1, ...
        'wrap', [0 0 0], 'prefix', 'r');
    spm_reslice(inputImages, flags);

    reslicedGM = fullfile(outDir, 'rGM_mask.nii');
    reslicedArtifact = fullfile(outDir, 'rArtifact_mask.nii');
    assert(exist(reslicedGM, 'file') == 2 && exist(reslicedArtifact, 'file') == 2, ...
        'Mask reslicing failed.');

    finalMask = fullfile(outDir, 'explicit_mask_GM_AND_ARTIFACT.nii');
    V1 = spm_vol(reslicedGM);
    V2 = spm_vol(reslicedArtifact);
    Vo = V1;
    Vo.fname = finalMask;
    Vo.dt = [spm_type('uint8') 0];
    flags = struct('dmtx', 0, 'mask', 0, 'interp', 0, 'dtype', spm_type('uint8'));
    spm_imcalc([V1; V2], Vo, '(i1>0) & (i2>0)', flags);
end

function run_firstlevel_spm(outDir, scans, rpFile, finalMask, condDefs, TR, HPF)
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {outDir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 60;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans(:);

    for i = 1:numel(condDefs)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).name = condDefs(i).name;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).onset = condDefs(i).onset;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).duration = condDefs(i).duration;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(i).orth = 1;
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {rpFile};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = HPF;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {finalMask};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(outDir, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);
end

function create_contrasts(spmmat)
    load(spmmat, 'SPM');
    SPM.xCon = [];
    definitions = {
        'OldValid_gt_NewValid', {'OldValid', 'NewValid'}, [1 -1];
        'NewValid_gt_OldValid', {'OldValid', 'NewValid'}, [-1 1];
        'AllValidRetrieval_gt_Baseline', {'OldValid', 'NewValid'}, [0.5 0.5]
    };

    for i = 1:size(definitions, 1)
        weights = build_contrast_vector(SPM, definitions{i, 2}, definitions{i, 3});
        contrast = spm_FcUtil('Set', definitions{i, 1}, 'T', 'c', weights', SPM.xX.xKXs);
        if i == 1
            SPM.xCon = contrast;
        else
            SPM.xCon(end+1) = contrast;
        end
    end

    save(spmmat, 'SPM');
    SPM = spm_contrasts(SPM, 1:numel(SPM.xCon));
    save(spmmat, 'SPM');
end

function weights = build_contrast_vector(SPM, conditionNames, conditionWeights)
    weights = zeros(1, size(SPM.xX.X, 2));
    designNames = string(SPM.xX.name);
    for i = 1:numel(conditionNames)
        idx = find(contains(designNames, string(conditionNames{i})) & contains(designNames, 'bf(1)'));
        assert(~isempty(idx), 'Condition not found in design matrix: %s', conditionNames{i});
        weights(idx) = conditionWeights(i);
    end
end

function reason = append_reason(reason, addition)
    if isempty(reason) || strcmpi(reason, 'pass')
        reason = addition;
    else
        reason = [reason ' ' addition];
    end
end

function rows = add_log(rows, analysis, stim, subject, status, reason, ...
        nScans, runLength, maxOnset, nHit, nMiss, nFA, nCR, nInvalidOld, ...
        nInvalidNew, nValidRetrieval, nValidOld, nValidNew, qcPass, outputDir)
    rows(end+1, :) = {
        analysis, stim, subject, status, reason, ...
        nScans, runLength, maxOnset, ...
        nHit, nMiss, nFA, nCR, nInvalidOld, nInvalidNew, ...
        nValidRetrieval, nValidOld, nValidNew, qcPass, outputDir
    };
end

function tableOut = cell2table_or_empty(rows, variables)
    if isempty(rows)
        tableOut = cell2table(cell(0, numel(variables)), 'VariableNames', variables);
    else
        tableOut = cell2table(rows, 'VariableNames', variables);
    end
end
