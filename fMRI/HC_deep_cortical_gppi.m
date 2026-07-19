%% Healthy-control deep-to-cortical gPPI
% HRF-deconvolved ROI-to-ROI gPPI for healthy-control old/new retrieval.
% Each subject directory must contain functional NIfTI data, motion parameters,
% and an onset_HCxx.xlsx file with an all_glm_events sheet.

clear; clc;
rng(20260614);

SPM_DIR = select_directory('Select the SPM12 directory');
RAW_ROOT = select_directory('Select the HC functional-data root directory');
SEED_DIR = select_directory('Select the seed ROI directory');
TARGET_DIR = select_directory('Select the cortical target ROI directory');
OUTPUT_PARENT = select_directory('Select the output parent directory');

addpath(SPM_DIR);
assert(exist('spm', 'file') == 2, 'SPM12 was not found in the selected directory.');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

RESULT_ROOT = fullfile(OUTPUT_PARENT, 'hc_deep_cortical_gppi');
if ~exist(RESULT_ROOT, 'dir')
    mkdir(RESULT_ROOT);
end
OUT_XLSX = fullfile(RESULT_ROOT, 'hc_deep_cortical_gppi_results.xlsx');

TR = 2;
HPF_CUTOFF_SEC = 128;

FUNC_PATTERNS = {
    '^task_sub\d+\.nii$'
    '^task_HC\d+\.nii$'
    '^swr.*\.nii$'
    '^wr.*\.nii$'
};

EXCLUDE_SUBJECTS = {};

MIN_SEED_VOXELS = 1;
MIN_TARGET_VOXELS = 5;
MIN_SEED_COVERAGE_FRAC = 0;
MIN_TARGET_COVERAGE_FRAC = 0;
APPLY_ARTIFACT_MASK_TO_SEED = false;
APPLY_ARTIFACT_MASK_TO_TARGET = false;

MIN_VALID_RETRIEVAL = 45;
MIN_VALID_OLD = 12;
MIN_VALID_NEW = 12;

SKIP_SAME_ROI = true;
DECONV_RIDGE_LAMBDA = 0.20;
MIN_CIRCUIT_EDGES = 2;

PRIMARY_MEASURES = {'ppi_new', 'ppi_old_minus_new'};
ALL_MEASURES = {'ppi_all', 'ppi_old', 'ppi_new', 'ppi_old_minus_new'};
N_PERM = 10000;
RUN_PERMUTATION = true;


SEED_SYSTEM_RULES = {
    'SN',                  {'^SN(_|$)', '^SN_L$', '^SN_R$'};
    'STN',                 {'^STN(_|$)', '^STN_L$', '^STN_R$'};
    'DBSTargetNuclei',     {'^SN(_|$)', '^STN(_|$)'};
    'GPi',                 {'^GPi(_|$)', '^Pallidum_Internum'};
    'GPe',                 {'^GPe(_|$)', '^Pallidum_Externum'};
    'Putamen',             {'^Putamen(_|$)'};
    'Caudate',             {'^Caudate(_|$)'};
    'NAcc',                {'^NAC(_|$)', '^NAcc(_|$)', '^Accumbens'};
    'BasalGanglia',        {'^Putamen(_|$)', '^Caudate(_|$)', '^GPi(_|$)', '^GPe(_|$)', '^Pallid', '^NAC(_|$)', '^NAcc(_|$)', '^Accumbens'};
    'ThalamicNuclei',      {'^Thalam', '^VIM(_|$)', '^Vim(_|$)'};
    'Hippocampus',         {'^Hippocampus(_|$)', '^Hippocampal', '^HIP(_|$)'};
    'PHG',                 {'^PHIP(_|$)', '^ParaHippocampal', '^Parahippocampal'};
    'HippocampalSystem',   {'^Hippocampus(_|$)', '^Hippocampal', '^HIP(_|$)', '^PHIP(_|$)', '^ParaHippocampal', '^Parahippocampal'};
    'Amygdala',            {'^Amygdala(_|$)'};
    'DeepNucleiAll',       {'^SN(_|$)', '^STN(_|$)', '^Putamen(_|$)', '^Caudate(_|$)', '^GPi(_|$)', '^GPe(_|$)', '^Pallid', '^NAC(_|$)', '^NAcc(_|$)', '^Accumbens', '^Thalam', '^VIM(_|$)', '^Vim(_|$)'}
};

TARGET_SYSTEM_RULES = {
    'VisualOccipital',     {'^VisualOccipital(_|$)', '^VisualOccipital_'};
    'Precuneus',           {'^Precuneus(_|$)', '^Precuneus_'};
    'PosteriorRetrieval',  {'^VisualOccipital(_|$)', '^VisualOccipital_', '^Precuneus(_|$)', '^Precuneus_'};
    'MotorSensorimotor',   {'^MotorSensorimotor(_|$)', '^MotorSensorimotor_'};
    'Parietal',            {'^Parietal(_|$)', '^Parietal_', '^Angular(_|$)', '^Angular_', '^SupraMarginal(_|$)', '^SupraMarginal_'};
    'Fusiform',            {'^Fusiform(_|$)', '^Fusiform_'};
    'Temporal',            {'^Temporal(_|$)', '^Temporal_'};
    'DLPFC',               {'^DLPFC(_|$)', '^DLPFC_'};
    'InferiorFrontal',     {'^InferiorFrontal(_|$)', '^InferiorFrontal_'};
    'MedialFrontal',       {'^MedialFrontal(_|$)', '^MedialFrontal_'};
    'FrontalControl',      {'^DLPFC(_|$)', '^DLPFC_', '^InferiorFrontal(_|$)', '^InferiorFrontal_', '^MedialFrontal(_|$)', '^MedialFrontal_'}
};

PRIMARY_CIRCUITS = {
    'SN_to_VisualOccipital',                  'SN',                 'VisualOccipital';
    'STN_to_VisualOccipital',                 'STN',                'VisualOccipital';
    'DBSTargetNuclei_to_VisualOccipital',     'DBSTargetNuclei',    'VisualOccipital';
    'SN_to_PosteriorRetrieval',               'SN',                 'PosteriorRetrieval';
    'DBSTargetNuclei_to_PosteriorRetrieval',  'DBSTargetNuclei',    'PosteriorRetrieval';
    'HippocampalSystem_to_FrontalControl',    'HippocampalSystem',  'FrontalControl'
};

SECONDARY_CIRCUITS = {
    'SN_to_Precuneus',                        'SN',                 'Precuneus';
    'STN_to_Precuneus',                       'STN',                'Precuneus';
    'SN_to_FrontalControl',                   'SN',                 'FrontalControl';
    'STN_to_FrontalControl',                  'STN',                'FrontalControl';
    'DBSTargetNuclei_to_FrontalControl',      'DBSTargetNuclei',    'FrontalControl';
    'BasalGanglia_to_FrontalControl',         'BasalGanglia',       'FrontalControl';
    'BasalGanglia_to_VisualOccipital',        'BasalGanglia',       'VisualOccipital';
    'DeepNucleiAll_to_VisualOccipital',       'DeepNucleiAll',      'VisualOccipital';
    'DeepNucleiAll_to_PosteriorRetrieval',    'DeepNucleiAll',      'PosteriorRetrieval';
    'ThalamicNuclei_to_FrontalControl',       'ThalamicNuclei',     'FrontalControl';
    'HippocampalSystem_to_VisualOccipital',   'HippocampalSystem',  'VisualOccipital';
    'HippocampalSystem_to_PosteriorRetrieval','HippocampalSystem',  'PosteriorRetrieval';
    'PHG_to_FrontalControl',                  'PHG',                'FrontalControl';
    'PHG_to_MedialFrontal',                   'PHG',                'MedialFrontal';
    'PHG_to_DLPFC',                           'PHG',                'DLPFC';
    'PHG_to_InferiorFrontal',                 'PHG',                'InferiorFrontal';
    'Hippocampus_to_VisualOccipital',         'Hippocampus',        'VisualOccipital';
    'SN_to_MotorSensorimotor',                'SN',                 'MotorSensorimotor';
    'STN_to_MotorSensorimotor',               'STN',                'MotorSensorimotor';
    'DBSTargetNuclei_to_MotorSensorimotor',   'DBSTargetNuclei',    'MotorSensorimotor'
};


fprintf('\nLoading seed ROIs:\n%s\n', SEED_DIR);
seedList = load_roi_list(SEED_DIR);
fprintf('Seeds found: %d\n', numel(seedList));

fprintf('\nLoading cortical target ROIs:\n%s\n', TARGET_DIR);
targetList = load_roi_list(TARGET_DIR);
fprintf('Targets found: %d\n', numel(targetList));

seedMembership = build_system_membership(seedList, SEED_SYSTEM_RULES, 'seed');
targetMembership = build_system_membership(targetList, TARGET_SYSTEM_RULES, 'target');

sourceSeedSystems = unique(string([PRIMARY_CIRCUITS(:,2); SECONDARY_CIRCUITS(:,2)]));
sourceTargetSystems = unique(string([PRIMARY_CIRCUITS(:,3); SECONDARY_CIRCUITS(:,3)]));
seedNamesToKeep = membership_names(seedMembership, sourceSeedSystems);
targetNamesToKeep = membership_names(targetMembership, sourceTargetSystems);

if isempty(seedNamesToKeep)
    error('No seed ROI matched SEED_SYSTEM_RULES. Check seed ROI filenames.');
end
if isempty(targetNamesToKeep)
    error('No target ROI matched TARGET_SYSTEM_RULES. Check target ROI filenames.');
end

seedList = filter_roi_list_by_names(seedList, seedNamesToKeep);
targetList = filter_roi_list_by_names(targetList, targetNamesToKeep);
seedMembership = build_system_membership(seedList, SEED_SYSTEM_RULES, 'seed');
targetMembership = build_system_membership(targetList, TARGET_SYSTEM_RULES, 'target');

fprintf('Seeds kept: %d\n', numel(seedList));
fprintf('Targets kept: %d\n', numel(targetList));

fprintf('\nBuilding HC session table...\n');
sessionTable = build_hc_session_table(RAW_ROOT, EXCLUDE_SUBJECTS);
fprintf('HC sessions found: %d\n', height(sessionTable));


edgeRows = {};
edgeVars = {
    'subject', 'stim', 'target_group', ...
    'seed', 'target_roi', ...
    'ppi_all', 'ppi_old', 'ppi_new', 'ppi_old_minus_new', ...
    'n_valid_retrieval', 'n_valid_old', 'n_valid_new', ...
    'n_seed_vox', 'n_target_vox', ...
    'seed_coverage_frac', 'target_coverage_frac', ...
    'func_pattern_used', 'deconv_lambda', 'status', 'reason'
};

roiQcRows = {};
roiQcVars = {
    'subject', 'stim', 'roi_role', 'roi_name', ...
    'n_vox', 'coverage_frac', 'status', 'reason'
};

sessionRows = {};
sessionVars = {
    'subject', 'stim', 'target_group', ...
    'status', 'reason', ...
    'n_scans', 'n_valid_retrieval', 'n_valid_old', 'n_valid_new', ...
    'func_pattern_used'
};

fprintf('\n=================================================\n');
fprintf('Running HC ROI-to-ROI HRF-deconvolved gPPI...\n');
fprintf('=================================================\n');

for r = 1:height(sessionTable)

    subj = char(sessionTable.subject(r));
    stim = char(sessionTable.stim(r));
    targetGroup = char(sessionTable.target_group(r));
    funcDir = char(sessionTable.func_dir(r));

    fprintf('\nSession: %s | %s\n', subj, funcDir);

    try
        [scans, patternUsed] = find_func_scans(funcDir, FUNC_PATTERNS);
        nScans = numel(scans);

        rpFile = find_latest_rp_file(funcDir);
        if isempty(rpFile)
            error('rp*.txt not found');
        end

        onsetFile = find_hc_onset_file(funcDir, subj);
        if isempty(onsetFile)
            error('onset_HCxx.xlsx not found');
        end

        events = read_onset_events(onsetFile);
        counts = get_retrieval_counts(events);

        if counts.nValidRetrieval < MIN_VALID_RETRIEVAL || ...
           counts.nValidOld < MIN_VALID_OLD || ...
           counts.nValidNew < MIN_VALID_NEW

            reason = sprintf('trial QC failed: valid=%d, old=%d, new=%d', ...
                counts.nValidRetrieval, counts.nValidOld, counts.nValidNew);

            sessionRows = add_session_row(sessionRows, subj, stim, targetGroup, ...
                'skipped_QC', reason, nScans, counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, patternUsed);

            fprintf('[SKIP QC] %s\n', reason);
            continue;
        end

        motion = read_motion_regressors(rpFile, nScans);
        psych = build_psychological_regressors(events, nScans, TR);
        nuisance = build_nuisance_matrix(psych, motion, nScans, TR, HPF_CUTOFF_SEC);

        cacheDir = fullfile(RESULT_ROOT, 'ROI_resliced_cache', subj);
        if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end

        seedCache = containers.Map();
        targetCache = containers.Map();

        for si = 1:numel(seedList)
            roi = seedList(si);
            [ts, nVox, coverageFrac, status, reason] = extract_roi_ts_for_session( ...
                roi.file, roi.name, scans, funcDir, cacheDir, ...
                MIN_SEED_VOXELS, MIN_SEED_COVERAGE_FRAC, APPLY_ARTIFACT_MASK_TO_SEED);

            seedCache(roi.name) = struct('ts', ts, 'nVox', nVox, ...
                'coverageFrac', coverageFrac, 'status', status, 'reason', reason);

            roiQcRows = add_roi_qc_row(roiQcRows, subj, stim, 'seed', roi.name, ...
                nVox, coverageFrac, status, reason);
        end

        for ti = 1:numel(targetList)
            roi = targetList(ti);
            [ts, nVox, coverageFrac, status, reason] = extract_roi_ts_for_session( ...
                roi.file, roi.name, scans, funcDir, cacheDir, ...
                MIN_TARGET_VOXELS, MIN_TARGET_COVERAGE_FRAC, APPLY_ARTIFACT_MASK_TO_TARGET);

            targetCache(roi.name) = struct('ts', ts, 'nVox', nVox, ...
                'coverageFrac', coverageFrac, 'status', status, 'reason', reason);

            roiQcRows = add_roi_qc_row(roiQcRows, subj, stim, 'target', roi.name, ...
                nVox, coverageFrac, status, reason);
        end

        for si = 1:numel(seedList)

            seedName = seedList(si).name;
            seedInfo = seedCache(seedName);

            if ~strcmp(seedInfo.status, 'ok')
                continue;
            end

            for ti = 1:numel(targetList)

                targetName = targetList(ti).name;

                if SKIP_SAME_ROI && strcmpi(seedName, targetName)
                    continue;
                end

                targetInfo = targetCache(targetName);

                if ~strcmp(targetInfo.status, 'ok')
                    continue;
                end

                try
                    beta = fit_gppi_edge(seedInfo.ts, targetInfo.ts, psych, nuisance, TR, DECONV_RIDGE_LAMBDA);

                    edgeRows = add_edge_row(edgeRows, subj, stim, targetGroup, ...
                        seedName, targetName, ...
                        beta.ppi_all, beta.ppi_old, beta.ppi_new, beta.ppi_old_minus_new, ...
                        counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, ...
                        seedInfo.nVox, targetInfo.nVox, ...
                        seedInfo.coverageFrac, targetInfo.coverageFrac, ...
                        patternUsed, DECONV_RIDGE_LAMBDA, 'success', '');

                catch ME
                    edgeRows = add_edge_row(edgeRows, subj, stim, targetGroup, ...
                        seedName, targetName, NaN, NaN, NaN, NaN, ...
                        counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, ...
                        seedInfo.nVox, targetInfo.nVox, ...
                        seedInfo.coverageFrac, targetInfo.coverageFrac, ...
                        patternUsed, DECONV_RIDGE_LAMBDA, 'failed', ME.message);
                end
            end
        end

        sessionRows = add_session_row(sessionRows, subj, stim, targetGroup, ...
            'success', '', nScans, counts.nValidRetrieval, counts.nValidOld, counts.nValidNew, patternUsed);

    catch ME
        sessionRows = add_session_row(sessionRows, subj, stim, targetGroup, ...
            'failed', ME.message, NaN, NaN, NaN, NaN, '');

        fprintf('[FAILED] %s | %s\n', subj, ME.message);
    end
end

edgeTable = cell2table_or_empty(edgeRows, edgeVars);
roiQcTable = cell2table_or_empty(roiQcRows, roiQcVars);
sessionLog = cell2table_or_empty(sessionRows, sessionVars);


fprintf('\n=================================================\n');
fprintf('Building HC circuit-level betas...\n');
fprintf('=================================================\n');

primaryCircuitTable = circuits_cell_to_table(PRIMARY_CIRCUITS, true);
secondaryCircuitTable = circuits_cell_to_table(SECONDARY_CIRCUITS, false);
circuitDefTable = [primaryCircuitTable; secondaryCircuitTable];

circuitTable = build_circuit_beta_table(edgeTable, circuitDefTable, seedMembership, targetMembership, ALL_MEASURES, MIN_CIRCUIT_EDGES);


fprintf('\n=================================================\n');
fprintf('Running HC one-sample circuit and edge statistics...\n');
fprintf('=================================================\n');

primaryCircuitStats = run_hc_one_sample_circuit_stats(circuitTable, PRIMARY_MEASURES, true);
primaryCircuitStats = add_fdr_one_family(primaryCircuitStats, 'p', 'p_fdr_primary');
if RUN_PERMUTATION && height(primaryCircuitStats) > 0
    primaryCircuitStats = add_signflip_maxstat_circuit(primaryCircuitStats, circuitTable, PRIMARY_MEASURES, N_PERM);
end

allCircuitStats = run_hc_one_sample_circuit_stats(circuitTable, ALL_MEASURES, false);
allCircuitStats = add_fdr_by_measure_family(allCircuitStats, 'p', 'p_fdr_by_measure');

primaryEdgeStats = run_hc_edge_followup_one_sample(edgeTable, primaryCircuitTable, seedMembership, targetMembership, PRIMARY_MEASURES);
primaryEdgeStats = add_fdr_by_circuit_measure(primaryEdgeStats, 'p', 'p_fdr_within_circuit_measure');


fprintf('\n=================================================\n');
fprintf('Saving outputs...\n');
fprintf('=================================================\n');

if exist(OUT_XLSX, 'file') == 2
    delete(OUT_XLSX);
end

writetable(edgeTable, OUT_XLSX, 'Sheet', 'edge_betas');
writetable(circuitTable, OUT_XLSX, 'Sheet', 'circuit_betas');
writetable(circuitDefTable, OUT_XLSX, 'Sheet', 'circuit_definitions');
writetable(seedMembership, OUT_XLSX, 'Sheet', 'seed_system_membership');
writetable(targetMembership, OUT_XLSX, 'Sheet', 'target_system_membership');
writetable(sessionLog, OUT_XLSX, 'Sheet', 'session_log');
writetable(roiQcTable, OUT_XLSX, 'Sheet', 'roi_qc');
writetable(primaryCircuitStats, OUT_XLSX, 'Sheet', 'HC_primary_circuit_stats');
writetable(allCircuitStats, OUT_XLSX, 'Sheet', 'HC_all_circuit_stats');
writetable(primaryEdgeStats, OUT_XLSX, 'Sheet', 'HC_primary_edge_stats');

notes = {
    'Analysis label', 'HC deep-to-cortical HRF-deconvolved ROI-to-ROI gPPI';
    'Raw root', RAW_ROOT;
    'Seed ROI directory', SEED_DIR;
    'Target ROI directory', TARGET_DIR;
    'Functional image priority', strjoin(FUNC_PATTERNS, ' then ');
    'Primary measures', strjoin(PRIMARY_MEASURES, ', ');
    'Primary correction family', 'six planned deep-to-cortical circuits x two primary measures';
    'Permutation correction', sprintf('sign-flip max-stat across primary family; N_PERM=%d; enabled=%d', N_PERM, RUN_PERMUTATION);
    'Deconvolution', sprintf('ridge-regularized canonical HRF inversion; lambda=%.4f', DECONV_RIDGE_LAMBDA);
    'Trial model', 'Hit+Miss=Old_valid; FalseAlarm+CorrectRejection=New_valid; Encoding/Interference/Invalid/ResponseOld/ResponseNew/motion/DCT as nuisance';
    'Output file', OUT_XLSX
};
notesTable = cell2table(notes, 'VariableNames', {'item', 'value'});
writetable(notesTable, OUT_XLSX, 'Sheet', 'notes');

fprintf('\n=============================================\n');
fprintf('HC deep-cortical gPPI finished.\n');
fprintf('Results saved to:\n%s\n', OUT_XLSX);
fprintf('=============================================\n');


function roiList = load_roi_list(roiDir)

    if exist(roiDir, 'dir') ~= 7
        error('ROI directory not found: %s', roiDir);
    end

    d = dir(fullfile(roiDir, '*.nii'));
    roiList = struct('name', {}, 'file', {});

    for i = 1:numel(d)
        [~, nm] = fileparts(d(i).name);
        roiList(end+1).name = matlab.lang.makeValidName(nm);
        roiList(end).file = fullfile(d(i).folder, d(i).name);
    end

    if isempty(roiList)
        error('No .nii ROI masks found in %s', roiDir);
    end
end

function roiList = filter_roi_list_by_names(roiList, namesToKeep)

    namesToKeep = string(namesToKeep(:));
    keep = false(numel(roiList), 1);

    for i = 1:numel(roiList)
        keep(i) = any(string(roiList(i).name) == namesToKeep);
    end

    roiList = roiList(keep);
end

function names = membership_names(membership, systems)

    if isempty(membership)
        names = strings(0, 1);
        return;
    end

    systems = string(systems(:));
    idx = ismember(string(membership.system_name), systems);
    names = unique(string(membership.roi_name(idx)));
end

function membership = build_system_membership(roiList, systemRules, role)

    rows = {};
    vars = {'roi_role', 'system_name', 'roi_name', 'matched'};

    for s = 1:size(systemRules, 1)
        sysName = string(systemRules{s, 1});
        patterns = systemRules{s, 2};

        for r = 1:numel(roiList)
            roiName = string(roiList(r).name);
            isHit = false;

            for p = 1:numel(patterns)
                if ~isempty(regexp(char(roiName), patterns{p}, 'once'))
                    isHit = true;
                    break;
                end
            end

            if isHit
                rows(end+1, :) = {char(role), char(sysName), char(roiName), true};
            end
        end
    end

    membership = cell2table_or_empty(rows, vars);
end

function T = circuits_cell_to_table(C, isPrimary)

    if isempty(C)
        T = cell2table(cell(0, 5), 'VariableNames', {'circuit', 'seed_system', 'target_system', 'is_primary', 'n_expected'});
        return;
    end

    n = size(C, 1);
    T = table();
    T.circuit = strings(n, 1);
    T.seed_system = strings(n, 1);
    T.target_system = strings(n, 1);
    T.is_primary = false(n, 1);
    T.n_expected = NaN(n, 1);

    for i = 1:n
        T.circuit(i) = string(C{i, 1});
        T.seed_system(i) = string(C{i, 2});
        T.target_system(i) = string(C{i, 3});
        T.is_primary(i) = isPrimary;
    end
end

function T = build_hc_session_table(rawRoot, excludeSubjects)

    if exist(rawRoot, 'dir') ~= 7
        error('RAW_ROOT not found: %s', rawRoot);
    end

    d = dir(fullfile(rawRoot, 'sub*'));
    d = d([d.isdir]);

    rows = {};
    for i = 1:numel(d)
        subj = normalize_subject_id(d(i).name);
        if strlength(subj) == 0
            continue;
        end
        rows(end+1, :) = {char(subj), 'HC_task', 'HC', fullfile(d(i).folder, d(i).name)};
    end

    T = cell2table_or_empty(rows, {'subject','stim','target_group','func_dir'});

    exSubs = strings(numel(excludeSubjects), 1);
    for i = 1:numel(excludeSubjects)
        exSubs(i) = normalize_subject_id(excludeSubjects{i});
    end
    if ~isempty(T) && ~isempty(exSubs)
        T = T(~ismember(string(T.subject), exSubs), :);
    end

    if height(T) == 0
        error('No HC subject folders found in %s', rawRoot);
    end
end

function sub = normalize_subject_id(x)

    x = char(string(x));
    x = strtrim(lower(x));
    tok = regexp(x, '\d+', 'match');

    if isempty(tok)
        sub = string(x);
    else
        sub = string(sprintf('sub%02d', str2double(tok{1})));
    end
end

function [scans, patternUsed] = find_func_scans(funcDir, patterns)

    scans = {};
    patternUsed = '';

    if exist(funcDir, 'dir') ~= 7
        error('Functional directory not found: %s', funcDir);
    end

    for p = 1:numel(patterns)
        pat = patterns{p};
        raw = spm_select('ExtFPList', funcDir, pat, Inf);

        if isempty(raw)
            continue;
        end

        tmp = cellstr(raw);
        tmp = tmp(~cellfun(@(x) isempty(strtrim(x)), tmp));
        keep = false(numel(tmp), 1);

        for i = 1:numel(tmp)
            one = tmp{i};
            commaPos = strfind(one, ',');
            if ~isempty(commaPos)
                fileOnly = one(1:commaPos(1)-1);
            else
                fileOnly = one;
            end
            keep(i) = exist(fileOnly, 'file') == 2;
        end

        tmp = tmp(keep);
        if ~isempty(tmp)
            scans = tmp;
            patternUsed = pat;
            return;
        end
    end

    error('No valid functional scans found in %s using patterns: %s', funcDir, strjoin(patterns, ', '));
end

function rpFile = find_latest_rp_file(funcDir)

    rpFile = '';
    d = dir(fullfile(funcDir, 'rp*.txt'));

    if isempty(d)
        return;
    end

    [~, idx] = max([d.datenum]);
    rpFile = fullfile(d(idx).folder, d(idx).name);
end

function onsetFile = find_hc_onset_file(funcDir, subj)

    n = regexp(char(subj), '\d+', 'match', 'once');
    onsetFile = '';

    candidates = {
        fullfile(funcDir, sprintf('onset_HC%s.xlsx', n));
        fullfile(funcDir, sprintf('onset_Hc%s.xlsx', n));
        fullfile(funcDir, sprintf('onset_hc%s.xlsx', n));
        fullfile(funcDir, sprintf('onset_%s.xlsx', subj))
    };

    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            onsetFile = candidates{i};
            return;
        end
    end

    d = dir(fullfile(funcDir, 'onset*.xlsx'));
    if isempty(d)
        return;
    end

    for i = 1:numel(d)
        if contains(d(i).name, n)
            onsetFile = fullfile(d(i).folder, d(i).name);
            return;
        end
    end

    onsetFile = fullfile(d(1).folder, d(1).name);
end

function events = read_onset_events(onsetFile)

    events = readtable(onsetFile, 'Sheet', 'all_glm_events', 'VariableNamingRule', 'preserve');
    events.trial_type = string(events.trial_type);
    events.onset = double(events.onset);

    if ~ismember('duration', events.Properties.VariableNames)
        events.duration = zeros(height(events), 1);
    else
        events.duration = double(events.duration);
    end

    events = events(~isnan(events.onset), :);
end

function counts = get_retrieval_counts(events)

    tt = string(events.trial_type);
    nHit = sum(tt == "Hit");
    nMiss = sum(tt == "Miss");
    nFA = sum(tt == "FalseAlarm");
    nCR = sum(tt == "CorrectRejection");

    counts.nValidOld = nHit + nMiss;
    counts.nValidNew = nFA + nCR;
    counts.nValidRetrieval = counts.nValidOld + counts.nValidNew;
end

function motion = read_motion_regressors(rpFile, nScans)

    rp = load(rpFile);

    if size(rp, 1) > nScans
        rp = rp(1:nScans, :);
    elseif size(rp, 1) < nScans
        error('rp rows (%d) < nScans (%d)', size(rp, 1), nScans);
    end

    motion = rp(:, 1:min(6, size(rp, 2)));
end

function psych = build_psychological_regressors(events, nScans, TR)

    tt = string(events.trial_type);

    idxAll = tt == "Hit" | tt == "Miss" | tt == "FalseAlarm" | tt == "CorrectRejection";
    idxOld = tt == "Hit" | tt == "Miss";
    idxNew = tt == "FalseAlarm" | tt == "CorrectRejection";

    psych.allValid_neural = make_neural_regressor(events, idxAll, nScans, TR);
    psych.oldValid_neural = make_neural_regressor(events, idxOld, nScans, TR);
    psych.newValid_neural = make_neural_regressor(events, idxNew, nScans, TR);

    psych.allValid_bold = convolve_regressor(psych.allValid_neural, TR);
    psych.oldValid_bold = convolve_regressor(psych.oldValid_neural, TR);
    psych.newValid_bold = convolve_regressor(psych.newValid_neural, TR);

    psych.encoding_neural = make_neural_regressor(events, tt == "Encoding", nScans, TR);
    psych.interference_neural = make_neural_regressor(events, tt == "Interference", nScans, TR);
    psych.invalid_neural = make_neural_regressor(events, tt == "InvalidOld" | tt == "InvalidNew", nScans, TR);
    psych.responseOld_neural = make_neural_regressor(events, tt == "ResponseOld", nScans, TR);
    psych.responseNew_neural = make_neural_regressor(events, tt == "ResponseNew", nScans, TR);

    psych.encoding_bold = convolve_regressor(psych.encoding_neural, TR);
    psych.interference_bold = convolve_regressor(psych.interference_neural, TR);
    psych.invalid_bold = convolve_regressor(psych.invalid_neural, TR);
    psych.responseOld_bold = convolve_regressor(psych.responseOld_neural, TR);
    psych.responseNew_bold = convolve_regressor(psych.responseNew_neural, TR);
end

function neural = make_neural_regressor(events, idx, nScans, TR)

    onsets = events.onset(idx);
    durations = events.duration(idx);
    neural = zeros(nScans, 1);

    for i = 1:numel(onsets)
        onset = onsets(i);
        dur = durations(i);

        if isnan(onset) || onset < 0
            continue;
        end

        if isnan(dur) || dur <= 0
            dur = TR;
        end

        startIdx = floor(onset / TR) + 1;
        endIdx = max(startIdx, ceil((onset + dur) / TR));
        startIdx = max(startIdx, 1);
        endIdx = min(endIdx, nScans);

        if startIdx <= nScans
            neural(startIdx:endIdx) = neural(startIdx:endIdx) + 1;
        end
    end
end

function boldReg = convolve_regressor(neural, TR)

    hrf = spm_hrf(TR);
    tmp = conv(neural(:), hrf(:));
    boldReg = tmp(1:numel(neural));
    boldReg = safe_zscore(boldReg);
end

function nuisance = build_nuisance_matrix(psych, motion, nScans, TR, hpfCutoff)

    X = [
        psych.encoding_bold(:), ...
        psych.interference_bold(:), ...
        psych.invalid_bold(:), ...
        psych.responseOld_bold(:), ...
        psych.responseNew_bold(:), ...
        motion
    ];

    dct = make_dct_basis(nScans, TR, hpfCutoff);
    X = [X, dct];

    nuisance = standardize_design(X);
end

function dct = make_dct_basis(nScans, TR, cutoffSec)

    runDur = nScans * TR;
    nBases = floor(2 * runDur / cutoffSec);

    if nBases <= 0
        dct = [];
        return;
    end

    t = (0:nScans-1)' / nScans;
    dct = zeros(nScans, nBases);

    for k = 1:nBases
        dct(:, k) = cos(pi * k * t);
    end
end

function [ts, nVox, coverageFrac, status, reason] = extract_roi_ts_for_session( ...
    roiFile, roiName, scans, funcDir, roiOutDir, minVox, minCoverage, applyArtifactMask)

    status = 'ok';
    reason = '';

    [roiMask, roiOriginalVox] = reslice_roi_to_func(roiFile, roiName, scans{1}, roiOutDir);

    if applyArtifactMask
        artMaskFile = find_artifact_mask(funcDir);
        if ~isempty(artMaskFile)
            [artMask, ~] = reslice_roi_to_func(artMaskFile, ['artifact_' roiName], scans{1}, roiOutDir);
            roiMask = roiMask & artMask;
        end
    end

    nVox = sum(roiMask(:) > 0);
    coverageFrac = nVox / max(roiOriginalVox, 1);

    if nVox < minVox
        ts = [];
        status = 'failed';
        reason = sprintf('ROI voxels %d < %d', nVox, minVox);
        return;
    end

    if minCoverage > 0 && coverageFrac < minCoverage
        fprintf('[WARNING] ROI coverage %.4f < %.4f, not excluded: %s\n', coverageFrac, minCoverage, roiName);
    end

    V = spm_vol(char(scans));
    idx = find(roiMask(:) > 0);
    ts = zeros(numel(V), 1);

    for i = 1:numel(V)
        Y = spm_read_vols(V(i));
        vals = Y(idx);
        vals = vals(isfinite(vals));

        if isempty(vals)
            ts(i) = NaN;
        else
            ts(i) = mean(vals);
        end
    end

    if any(isnan(ts))
        ts = fillmissing(ts, 'linear', 'EndValues', 'nearest');
    end

    ts = safe_zscore(ts);
end

function [maskBin, roiOriginalVox] = reslice_roi_to_func(roiFile, roiName, refScan, outDir)

    if isempty(roiFile) || exist(roiFile, 'file') ~= 2
        error('ROI file does not exist: %s', roiFile);
    end

    if isempty(refScan)
        error('Reference functional scan is empty.');
    end

    commaPos = strfind(refScan, ',');
    if ~isempty(commaPos)
        refFileOnly = refScan(1:commaPos(1)-1);
    else
        refFileOnly = refScan;
    end

    if exist(refFileOnly, 'file') ~= 2
        error('Reference functional scan does not exist: %s', refFileOnly);
    end

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    roiNameSafe = matlab.lang.makeValidName(roiName);
    roiCopy = fullfile(outDir, [roiNameSafe '.nii']);

    if exist(roiCopy, 'file') ~= 2
        copyfile(roiFile, roiCopy);
    end

    resliced = fullfile(outDir, ['r' roiNameSafe '.nii']);

    if exist(resliced, 'file') ~= 2
        P = char(refScan, [roiCopy ',1']);
        flags = struct('mask', 0, 'mean', 0, 'interp', 0, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
        spm_reslice(P, flags);
    end

    if exist(resliced, 'file') ~= 2
        error('Resliced ROI was not created: %s', resliced);
    end

    V0 = spm_vol(roiCopy);
    Y0 = spm_read_vols(V0);
    roiOriginalVox = sum(Y0(:) > 0);

    Vr = spm_vol(resliced);
    Yr = spm_read_vols(Vr);
    maskBin = Yr > 0.5;
end

function artFile = find_artifact_mask(funcDir)

    artFile = '';
    d = dir(fullfile(funcDir, 'glm_mask_*.nii'));

    if isempty(d)
        return;
    end

    artFile = fullfile(d(1).folder, d(1).name);
end

function beta = fit_gppi_edge(seedTs, targetTs, psych, nuisance, TR, lambda)

    y = safe_zscore(targetTs(:));
    seedBold = safe_zscore(seedTs(:));
    seedNeural = hrf_deconvolve_ridge(seedBold, TR, lambda);
    seedNeural = safe_zscore(seedNeural);

    allReg = safe_zscore(psych.allValid_bold(:));
    oldReg = safe_zscore(psych.oldValid_bold(:));
    newReg = safe_zscore(psych.newValid_bold(:));

    ppiAll = make_gppi_interaction(seedNeural, psych.allValid_neural(:), TR);
    ppiOld = make_gppi_interaction(seedNeural, psych.oldValid_neural(:), TR);
    ppiNew = make_gppi_interaction(seedNeural, psych.newValid_neural(:), TR);

    Xall = [ones(numel(y), 1), seedBold, allReg, ppiAll, nuisance];
    beta.ppi_all = get_beta_with_protected_columns(Xall, y, 4);

    Xoldnew = [ones(numel(y), 1), seedBold, oldReg, newReg, ppiOld, ppiNew, nuisance];
    beta.ppi_old = get_beta_with_protected_columns(Xoldnew, y, 5);
    beta.ppi_new = get_beta_with_protected_columns(Xoldnew, y, 6);
    beta.ppi_old_minus_new = beta.ppi_old - beta.ppi_new;
end

function neural = hrf_deconvolve_ridge(bold, TR, lambda)

    bold = safe_zscore(bold(:));
    n = numel(bold);
    h = spm_hrf(TR);
    h = h(:);

    H = zeros(n, n);
    for c = 1:n
        len = min(numel(h), n - c + 1);
        H(c:(c+len-1), c) = h(1:len);
    end

    if nargin < 3 || isempty(lambda)
        lambda = 0.20;
    end

    neural = (H' * H + lambda * eye(n)) \ (H' * bold);
    neural = safe_zscore(neural);
end

function ppiBold = make_gppi_interaction(seedNeural, psychNeural, TR)

    seedNeural = safe_zscore(seedNeural(:));
    psychNeural = double(psychNeural(:));

    psychCentered = psychNeural - mean(psychNeural, 'omitnan');
    interactionNeural = seedNeural .* psychCentered;

    hrf = spm_hrf(TR);
    tmp = conv(interactionNeural(:), hrf(:));
    ppiBold = tmp(1:numel(seedNeural));
    ppiBold = safe_zscore(ppiBold);
end

function betaVal = get_beta_with_protected_columns(Xraw, y, targetCol)

    keep = true(1, size(Xraw, 2));

    for c = 2:size(Xraw, 2)
        if std(Xraw(:, c), 'omitnan') < 1e-8
            keep(c) = false;
        end
    end

    if ~keep(targetCol)
        betaVal = NaN;
        return;
    end

    X = Xraw(:, keep);
    b = pinv(X) * y;

    rawIdx = find(keep);
    newIdx = find(rawIdx == targetCol);
    betaVal = b(newIdx);
end

function X = standardize_design(X)

    if isempty(X)
        return;
    end

    for c = 1:size(X, 2)
        X(:, c) = safe_zscore(X(:, c));
    end

    keep = true(1, size(X, 2));
    for c = 1:size(X, 2)
        if std(X(:, c), 'omitnan') < 1e-8
            keep(c) = false;
        end
    end
    X = X(:, keep);
end

function z = safe_zscore(x)

    x = double(x(:));
    mu = mean(x, 'omitnan');
    sd = std(x, 'omitnan');

    if sd < 1e-8 || isnan(sd)
        z = zeros(size(x));
    else
        z = (x - mu) ./ sd;
    end

    z(~isfinite(z)) = 0;
end

function T = cell2table_or_empty(rows, vars)

    if isempty(rows)
        T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
    else
        T = cell2table(rows, 'VariableNames', vars);
    end
end

function rows = add_edge_row(rows, subject, stim, targetGroup, seed, targetROI, ...
    ppiAll, ppiOld, ppiNew, ppiOldMinusNew, ...
    nValid, nOld, nNew, nSeedVox, nTargetVox, seedCoverage, targetCoverage, ...
    patternUsed, deconvLambda, status, reason)

    rows(end+1, :) = {
        subject, stim, targetGroup, ...
        seed, targetROI, ...
        ppiAll, ppiOld, ppiNew, ppiOldMinusNew, ...
        nValid, nOld, nNew, ...
        nSeedVox, nTargetVox, ...
        seedCoverage, targetCoverage, ...
        patternUsed, deconvLambda, status, reason
    };
end

function rows = add_roi_qc_row(rows, subject, stim, role, roiName, nVox, coverageFrac, status, reason)

    rows(end+1, :) = {subject, stim, role, roiName, nVox, coverageFrac, status, reason};
end

function rows = add_session_row(rows, subject, stim, targetGroup, status, reason, nScans, nValid, nOld, nNew, patternUsed)

    rows(end+1, :) = {subject, stim, targetGroup, status, reason, nScans, nValid, nOld, nNew, patternUsed};
end

function circuitTable = build_circuit_beta_table(edgeTable, circuitDefTable, seedMembership, targetMembership, measures, minEdges)

    rows = {};
    vars = [{'subject','stim','target_group','circuit','seed_system','target_system','is_primary','n_edges','status','reason'}, measures];

    if isempty(edgeTable) || isempty(circuitDefTable)
        circuitTable = cell2table_or_empty(rows, vars);
        return;
    end

    sessionList = unique(edgeTable(:, {'subject','stim','target_group'}));

    for s = 1:height(sessionList)
        subj = string(sessionList.subject(s));
        stim = string(sessionList.stim(s));
        tg = string(sessionList.target_group(s));

        for c = 1:height(circuitDefTable)
            circuit = string(circuitDefTable.circuit(c));
            seedSystem = string(circuitDefTable.seed_system(c));
            targetSystem = string(circuitDefTable.target_system(c));
            isPrimary = logical(circuitDefTable.is_primary(c));

            seedNames = system_members(seedMembership, seedSystem);
            targetNames = system_members(targetMembership, targetSystem);

            idx = string(edgeTable.subject) == subj & ...
                  string(edgeTable.stim) == stim & ...
                  string(edgeTable.target_group) == tg & ...
                  string(edgeTable.status) == "success" & ...
                  ismember(string(edgeTable.seed), seedNames) & ...
                  ismember(string(edgeTable.target_roi), targetNames);

            nEdges = sum(idx);
            vals = cell(1, numel(measures));

            if nEdges < minEdges
                status = 'skipped';
                reason = sprintf('n_edges %d < %d', nEdges, minEdges);
                for m = 1:numel(measures)
                    vals{m} = NaN;
                end
            else
                status = 'success';
                reason = '';
                for m = 1:numel(measures)
                    x = double(edgeTable.(measures{m})(idx));
                    vals{m} = mean(x, 'omitnan');
                end
            end

            rows(end+1, :) = [{char(subj), char(stim), char(tg), char(circuit), char(seedSystem), char(targetSystem), isPrimary, nEdges, status, reason}, vals];
        end
    end

    circuitTable = cell2table_or_empty(rows, vars);
end

function names = system_members(membership, systemName)

    if isempty(membership)
        names = strings(0, 1);
        return;
    end

    idx = string(membership.system_name) == string(systemName);
    names = unique(string(membership.roi_name(idx)));
end

function statsTable = run_hc_one_sample_circuit_stats(circuitTable, measures, primaryOnly)

    rows = {};
    vars = {'family','measure','circuit','seed_system','target_system','is_primary','n','mean_beta','sd_beta','sem_beta','t','df','p','p_fdr_primary','p_fdr_by_measure','p_perm_max'};

    if isempty(circuitTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    if primaryOnly
        circuits = unique(circuitTable(circuitTable.is_primary == true, {'circuit','seed_system','target_system','is_primary'}));
        familyName = 'HC_primary_circuit_one_sample';
    else
        circuits = unique(circuitTable(:, {'circuit','seed_system','target_system','is_primary'}));
        familyName = 'HC_all_circuit_one_sample';
    end

    for m = 1:numel(measures)
        measure = measures{m};
        for c = 1:height(circuits)
            circuit = string(circuits.circuit(c));
            idx = string(circuitTable.circuit) == circuit & string(circuitTable.status) == "success";
            vals = double(circuitTable.(measure)(idx));
            vals = vals(isfinite(vals));
            [tval, df, pval] = one_sample_ttest(vals);
            n = numel(vals);
            rows(end+1, :) = { ...
                familyName, measure, char(circuit), char(circuits.seed_system(c)), char(circuits.target_system(c)), logical(circuits.is_primary(c)), ...
                n, mean(vals, 'omitnan'), std(vals, 'omitnan'), std(vals, 'omitnan') / max(sqrt(n), 1), ...
                tval, df, pval, NaN, NaN, NaN};
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function edgeStats = run_hc_edge_followup_one_sample(edgeTable, circuitDefTable, seedMembership, targetMembership, measures)

    rows = {};
    vars = {'family','measure','circuit','seed','target_roi','n','mean_beta','sd_beta','sem_beta','t','df','p','p_fdr_within_circuit_measure'};

    if isempty(edgeTable) || isempty(circuitDefTable)
        edgeStats = cell2table_or_empty(rows, vars);
        return;
    end

    for c = 1:height(circuitDefTable)
        circuit = string(circuitDefTable.circuit(c));
        seedSystem = string(circuitDefTable.seed_system(c));
        targetSystem = string(circuitDefTable.target_system(c));
        seedNames = system_members(seedMembership, seedSystem);
        targetNames = system_members(targetMembership, targetSystem);

        for m = 1:numel(measures)
            measure = measures{m};
            for si = 1:numel(seedNames)
                for ti = 1:numel(targetNames)
                    seed = seedNames(si);
                    targetROI = targetNames(ti);

                    vals = edgeTable.(measure)(string(edgeTable.seed) == seed & ...
                        string(edgeTable.target_roi) == targetROI & ...
                        string(edgeTable.status) == "success");
                    vals = double(vals);
                    vals = vals(isfinite(vals));
                    [tval, df, pval] = one_sample_ttest(vals);
                    n = numel(vals);

                    rows(end+1, :) = { ...
                        'HC_primary_edge_one_sample', measure, char(circuit), char(seed), char(targetROI), ...
                        n, mean(vals, 'omitnan'), std(vals, 'omitnan'), std(vals, 'omitnan') / max(sqrt(n), 1), ...
                        tval, df, pval, NaN};
                end
            end
        end
    end

    edgeStats = cell2table_or_empty(rows, vars);
end

function [tval, df, pval] = one_sample_ttest(vals)

    vals = vals(isfinite(vals));
    if numel(vals) < 3
        tval = NaN; df = NaN; pval = NaN;
        return;
    end

    [~, pval, ~, st] = ttest(vals);
    tval = st.tstat;
    df = st.df;
end

function T = add_fdr_one_family(T, pName, qName)

    if isempty(T) || ~ismember(pName, T.Properties.VariableNames)
        return;
    end

    T.(qName) = NaN(height(T), 1);
    idx = isfinite(T.(pName));
    T.(qName)(idx) = bh_fdr(T.(pName)(idx));
end

function T = add_fdr_by_measure_family(T, pName, qName)

    if isempty(T) || ~ismember(pName, T.Properties.VariableNames)
        return;
    end

    T.(qName) = NaN(height(T), 1);
    measures = unique(string(T.measure));

    for m = 1:numel(measures)
        idx = string(T.measure) == measures(m) & isfinite(T.(pName));
        T.(qName)(idx) = bh_fdr(T.(pName)(idx));
    end
end

function T = add_fdr_by_circuit_measure(T, pName, qName)

    if isempty(T) || ~ismember(pName, T.Properties.VariableNames)
        return;
    end

    T.(qName) = NaN(height(T), 1);
    keys = unique(strcat(string(T.circuit), "___", string(T.measure)));

    for k = 1:numel(keys)
        parts = split(keys(k), "___");
        idx = string(T.circuit) == parts(1) & string(T.measure) == parts(2) & isfinite(T.(pName));
        T.(qName)(idx) = bh_fdr(T.(pName)(idx));
    end
end

function q = bh_fdr(p)

    p = p(:);
    q = NaN(size(p));
    valid = isfinite(p);
    pv = p(valid);

    if isempty(pv)
        return;
    end

    [ps, order] = sort(pv);
    m = numel(ps);
    qs = ps .* m ./ (1:m)';

    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end

    qv = NaN(size(pv));
    qv(order) = min(qs, 1);
    q(valid) = qv;
end

function T = add_signflip_maxstat_circuit(T, circuitTable, measures, nPerm)

    if isempty(T)
        return;
    end

    fprintf('Sign-flip max-stat permutation for HC primary circuit tests: %d permutations\n', nPerm);
    obs = abs(T.t);
    maxStats = NaN(nPerm, 1);

    data = cell(height(T), 1);
    for r = 1:height(T)
        measure = char(T.measure(r));
        circuit = string(T.circuit(r));
        vals = double(circuitTable.(measure)(string(circuitTable.circuit) == circuit & string(circuitTable.status) == "success"));
        vals = vals(isfinite(vals));
        data{r} = vals(:);
    end

    for b = 1:nPerm
        maxT = 0;
        for r = 1:height(T)
            vals = data{r};
            if numel(vals) < 3
                continue;
            end
            signs = (rand(numel(vals), 1) > 0.5) * 2 - 1;
            vp = vals .* signs;
            [tval, ~, ~] = one_sample_ttest(vp);
            if isfinite(tval)
                maxT = max(maxT, abs(tval));
            end
        end
        maxStats(b) = maxT;
    end

    T.p_perm_max = NaN(height(T), 1);
    for r = 1:height(T)
        if isfinite(obs(r))
            T.p_perm_max(r) = (1 + sum(maxStats >= obs(r))) / (1 + sum(isfinite(maxStats)));
        end
    end

end

function folder = select_directory(prompt)
    folder = uigetdir(pwd, prompt);
    if isequal(folder, 0)
        error('Directory selection was cancelled.');
    end
end
