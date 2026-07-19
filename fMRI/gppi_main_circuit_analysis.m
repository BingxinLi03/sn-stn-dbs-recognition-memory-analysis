% gppi_main_circuit_analysis.m
% Main ROI-to-ROI HRF-deconvolved gPPI analysis with circuit-first inference.
% Valid old trials are Hit + Miss; valid new trials are FalseAlarm + CorrectRejection.
% Behavior metrics are supplied in a separate input table.
% ROI masks are external atlas-derived inputs and their provenance should be documented separately.
% Requires MATLAB R2022b or later and SPM12.

clear; clc;
rng(20260512);

SPM_DIR = request_existing_directory('Enter the full path to the SPM12 directory: ');
addpath(SPM_DIR);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

RAW_ROOT = request_existing_directory('Enter the functional-data root directory: ');
QC_LOG_FILE = request_existing_file('Enter the first-level session QC log (.xlsx): ');
BEHAVIOR_FILE = request_existing_file('Enter the behavior table (.xlsx): ');
SEED_DIR = request_existing_directory('Enter the seed-mask directory: ');
TARGET_DIR = request_existing_directory('Enter the target-mask directory: ');
RESULT_ROOT = request_output_directory('Enter the output directory: ');
OUT_XLSX = fullfile(RESULT_ROOT, 'gppi_main_circuit_analysis_results.xlsx');

TR = 2;
HPF_CUTOFF_SEC = 128;

STIM_OFF = 'stim0';
ACTIVE_STIMS = {'stim10', 'stim130'};
TARGET_GROUPS = {'SN', 'STN'};

EXCLUDE_SUBJECTS = {'sub36'};
EXCLUDE_SESSIONS = cell(0, 2);

FUNC_PATTERNS = {
    '^swr.*\.nii$'
    '^wr.*\.nii$'
};

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
PRIMARY_BEHAVIOR = 'false_alarm_reduction';
SECONDARY_BEHAVIORS = {'d_prime_improvement', 'hit_rate_change', 'Pr_improvement', 'bias_c_change'};
N_PERM = 10000;
RUN_PERMUTATION = true;

SEED_SYSTEM_RULES = {
    'SN',                  {'^SN_L$', '^SN_R$'};
    'STN',                 {'^STN_L$', '^STN_R$'};
    'SN_STN',              {'^SN_L$', '^SN_R$', '^STN_L$', '^STN_R$'};
    'PallidoStriatal',     {'^Putamen_L$', '^Putamen_R$', '^GPi_L$', '^GPi_R$', '^NAC_L$', '^NAC_R$'};
    'PHIP',                {'^PHIP_L$', '^PHIP_R$'};
    'SubcorticalAll',      {'^SN_', '^STN_', '^Putamen_', '^GPi_', '^NAC_', '^PHIP_'}
};

TARGET_SYSTEM_RULES = {
    'VisualOccipital',     {'^VisualOccipital_L$', '^VisualOccipital_R$'};
    'MotorSensorimotor',   {'^MotorSensorimotor_L$', '^MotorSensorimotor_R$'};
    'FrontalControl',      {'^MedialFrontal_', '^DLPFC_', '^InferiorFrontal_'};
    'Precuneus',           {'^Precuneus_L$', '^Precuneus_R$'};
    'PosteriorRetrieval',  {'^VisualOccipital_', '^Precuneus_'}
};

PRIMARY_CIRCUITS = {
    'SN_to_VisualOccipital',              'SN',              'VisualOccipital';
    'SN_to_MotorSensorimotor',            'SN',              'MotorSensorimotor';
    'PallidoStriatal_to_MotorSensorimotor','PallidoStriatal','MotorSensorimotor';
    'PallidoStriatal_to_FrontalControl',  'PallidoStriatal', 'FrontalControl';
    'PHIP_to_MotorSensorimotor',          'PHIP',            'MotorSensorimotor'
};

SECONDARY_CIRCUITS = {
    'SN_STN_to_VisualOccipital',          'SN_STN',          'VisualOccipital';
    'SubcorticalAll_to_MotorSensorimotor','SubcorticalAll',  'MotorSensorimotor';
    'SubcorticalAll_to_FrontalControl',   'SubcorticalAll',  'FrontalControl';
    'SN_to_PosteriorRetrieval',           'SN',              'PosteriorRetrieval';
    'PHIP_to_VisualOccipital',            'PHIP',            'VisualOccipital'
};

fprintf('\nLoading seed ROIs:\n%s\n', SEED_DIR);
seedList = load_roi_list(SEED_DIR);
fprintf('Seeds found: %d\n', numel(seedList));

fprintf('\nLoading target ROIs:\n%s\n', TARGET_DIR);
targetList = load_roi_list(TARGET_DIR);
fprintf('Targets found: %d\n', numel(targetList));

seedMembership = build_system_membership(seedList, SEED_SYSTEM_RULES, 'seed');
targetMembership = build_system_membership(targetList, TARGET_SYSTEM_RULES, 'target');

fprintf('\nReading behavior file:\n%s\n', BEHAVIOR_FILE);
[targetMap, behaviorTable] = read_behavior_file(BEHAVIOR_FILE);
behaviorTable = add_behavior_change_columns(behaviorTable, STIM_OFF, ACTIVE_STIMS);

fprintf('\nBuilding usable session table...\n');
sessionTable = build_session_table(RAW_ROOT, QC_LOG_FILE, targetMap, EXCLUDE_SUBJECTS, EXCLUDE_SESSIONS);
fprintf('Usable sessions: %d\n', height(sessionTable));

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
fprintf('Running ROI-to-ROI HRF-deconvolved gPPI edge estimation...\n');
fprintf('=================================================\n');

for r = 1:height(sessionTable)

    subj = char(sessionTable.subject(r));
    stim = char(sessionTable.stim(r));
    targetGroup = char(sessionTable.target_group(r));
    funcDir = char(sessionTable.func_dir(r));

    fprintf('\nSession: %s | %s | %s\n', subj, stim, targetGroup);

    try
        [scans, patternUsed] = find_func_scans(funcDir, FUNC_PATTERNS);
        nScans = numel(scans);

        rpFile = find_latest_rp_file(funcDir);
        if isempty(rpFile)
            error('rp*.txt not found');
        end

        onsetFile = find_new_onset_file(funcDir, subj);
        if isempty(onsetFile)
            error('onset_subxx.xlsx not found');
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

        cacheDir = fullfile(RESULT_ROOT, 'ROI_resliced_cache', stim, subj);
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

        fprintf('[FAILED] %s | %s | %s\n', subj, stim, ME.message);
    end
end

edgeTable = cell2table_or_empty(edgeRows, edgeVars);
roiQcTable = cell2table_or_empty(roiQcRows, roiQcVars);
sessionLog = cell2table_or_empty(sessionRows, sessionVars);

fprintf('\n=================================================\n');
fprintf('Building circuit-level betas...\n');
fprintf('=================================================\n');

primaryCircuitTable = circuits_cell_to_table(PRIMARY_CIRCUITS, true);
secondaryCircuitTable = circuits_cell_to_table(SECONDARY_CIRCUITS, false);
circuitDefTable = [primaryCircuitTable; secondaryCircuitTable];

circuitTable = build_circuit_beta_table(edgeTable, circuitDefTable, seedMembership, targetMembership, ALL_MEASURES, MIN_CIRCUIT_EDGES);

fprintf('\n=================================================\n');
fprintf('Running circuit-level group statistics...\n');
fprintf('=================================================\n');

primaryDeltaStats = run_primary_circuit_delta_interaction_stats( ...
    circuitTable, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF);
primaryDeltaStats = add_fdr_one_family(primaryDeltaStats, 'p', 'p_fdr_primary');

if RUN_PERMUTATION && height(primaryDeltaStats) > 0
    primaryDeltaStats = add_permutation_maxstat_delta_interaction( ...
        primaryDeltaStats, circuitTable, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF, N_PERM);
end

simpleDeltaStats = run_circuit_simple_delta_stats( ...
    circuitTable, ALL_MEASURES, ACTIVE_STIMS, STIM_OFF, TARGET_GROUPS);
simpleDeltaStats = add_fdr_by_measure_family(simpleDeltaStats, 'p', 'p_fdr_by_measure');

freqStats = run_circuit_frequency_stats(circuitTable, ALL_MEASURES, TARGET_GROUPS);
freqStats = add_fdr_by_measure_family(freqStats, 'p', 'p_fdr_by_measure');

fprintf('\n=================================================\n');
fprintf('Running circuit-level behavior models...\n');
fprintf('=================================================\n');

primaryBehaviorStats = run_primary_behavior_interaction_stats( ...
    circuitTable, behaviorTable, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF, PRIMARY_BEHAVIOR);
primaryBehaviorStats = add_fdr_one_family(primaryBehaviorStats, 'p_interaction', 'p_fdr_primary');

if RUN_PERMUTATION && height(primaryBehaviorStats) > 0
    primaryBehaviorStats = add_permutation_maxstat_behavior_interaction( ...
        primaryBehaviorStats, circuitTable, behaviorTable, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF, PRIMARY_BEHAVIOR, N_PERM);
end

secondaryBehaviorStats = run_secondary_behavior_stats( ...
    circuitTable, behaviorTable, ALL_MEASURES, ACTIVE_STIMS, STIM_OFF, SECONDARY_BEHAVIORS);
secondaryBehaviorStats = add_fdr_by_behavior_measure_family(secondaryBehaviorStats, 'p_interaction', 'p_fdr_by_behavior_measure');

fprintf('\n=================================================\n');
fprintf('Running edge-level follow-up within primary circuits...\n');
fprintf('=================================================\n');

edgeFollowupDelta = run_edge_followup_delta(edgeTable, primaryCircuitTable, seedMembership, targetMembership, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF);
edgeFollowupDelta = add_fdr_by_circuit_measure(edgeFollowupDelta, 'p', 'p_fdr_within_circuit_measure');

edgeFollowupBehavior = run_edge_followup_behavior(edgeTable, behaviorTable, primaryCircuitTable, seedMembership, targetMembership, PRIMARY_MEASURES, ACTIVE_STIMS, STIM_OFF, PRIMARY_BEHAVIOR);
edgeFollowupBehavior = add_fdr_by_circuit_measure(edgeFollowupBehavior, 'p_interaction', 'p_fdr_within_circuit_measure');

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
writetable(behaviorTable, OUT_XLSX, 'Sheet', 'behavior_with_changes');
writetable(primaryDeltaStats, OUT_XLSX, 'Sheet', 'primary_delta_interaction');
writetable(simpleDeltaStats, OUT_XLSX, 'Sheet', 'circuit_simple_delta');
writetable(freqStats, OUT_XLSX, 'Sheet', 'circuit_freq_10_130');
writetable(primaryBehaviorStats, OUT_XLSX, 'Sheet', 'primary_behavior_FAR');
writetable(secondaryBehaviorStats, OUT_XLSX, 'Sheet', 'secondary_behavior');
writetable(edgeFollowupDelta, OUT_XLSX, 'Sheet', 'edge_followup_delta');
writetable(edgeFollowupBehavior, OUT_XLSX, 'Sheet', 'edge_followup_FAR');

notes = {
    'Analysis label', 'Circuit-first HRF-deconvolved ROI-to-ROI gPPI';
    'Functional image priority', strjoin(FUNC_PATTERNS, ' then ');
    'Primary measures', strjoin(PRIMARY_MEASURES, ', ');
    'Primary behavior', PRIMARY_BEHAVIOR;
    'Primary behavior definition', 'false_alarm_rate_OFF - false_alarm_rate_active; positive values mean false alarm reduction/improvement';
    'Primary correction family', 'primary circuits x primary measures x active frequencies';
    'Permutation correction', sprintf('max-stat; N_PERM=%d; enabled=%d', N_PERM, RUN_PERMUTATION);
    'Deconvolution', sprintf('ridge-regularized canonical HRF inversion; lambda=%.4f', DECONV_RIDGE_LAMBDA);
    'Seed artifact mask intersection', num2str(APPLY_ARTIFACT_MASK_TO_SEED);
    'Target artifact mask intersection', num2str(APPLY_ARTIFACT_MASK_TO_TARGET);
    'Output file', OUT_XLSX
};
notesTable = cell2table(notes, 'VariableNames', {'item', 'value'});
writetable(notesTable, OUT_XLSX, 'Sheet', 'notes');

fprintf('\n=============================================\n');
fprintf('Circuit-first gPPI finished.\n');
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
        roiList(end+1).name = matlab.lang.makeValidName(nm); %#ok<AGROW>
        roiList(end).file = fullfile(d(i).folder, d(i).name);
    end

    if isempty(roiList)
        error('No .nii ROI masks found in %s', roiDir);
    end
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
                rows(end+1, :) = {char(role), char(sysName), char(roiName), true}; %#ok<AGROW>
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

function [targetMap, behTable] = read_behavior_file(behaviorFile)

    if exist(behaviorFile, 'file') ~= 2
        error('behavior.xlsx not found: %s', behaviorFile);
    end

    [~, sheets] = xlsfinfo(behaviorFile);
    targetMap = containers.Map();
    behTable = table();

    for s = 1:numel(sheets)

        try
            T = readtable(behaviorFile, 'Sheet', sheets{s}, 'VariableNamingRule', 'preserve');
        catch
            continue;
        end

        if isempty(T) || height(T) == 0
            continue;
        end

        varNames = string(T.Properties.VariableNames);
        normNames = normalize_var_names(varNames);

        subjectCol = find_first_col(normNames, {'subject','sub','subid','participant','participantid'});
        targetCol  = find_first_col(normNames, {'target','group','dbstarget'});
        stimCol    = find_first_col(normNames, {'stim','stimulation','condition','frequency','freq'});

        if ~isempty(subjectCol) && ~isempty(targetCol)
            subs = column_to_string(T.(char(varNames(subjectCol))));
            tars = column_to_string(T.(char(varNames(targetCol))));

            for i = 1:numel(subs)
                sub = normalize_subject_id(subs(i));
                tar = normalize_target_label(tars(i));

                if strlength(sub) == 0 || strlength(tar) == 0
                    continue;
                end

                if ~isKey(targetMap, char(sub))
                    targetMap(char(sub)) = char(tar);
                end
            end
        end

        metricCols = struct();
        metricCols.d_prime = find_first_col(normNames, {'dprime','d_prime','dprim','dprimescore','d'});
        metricCols.hit_rate = find_first_col(normNames, {'hitrate','hit_rate','hit','hr'});
        metricCols.false_alarm_rate = find_first_col(normNames, {'falsealarmrate','false_alarm_rate','far','fa_rate','falsealarm'});
        metricCols.Pr = find_first_col(normNames, {'pr','p_r'});
        metricCols.bias_c = find_first_col(normNames, {'biasc','bias_c','criterion','bias'});

        hasMetric = ~isempty(metricCols.d_prime) || ~isempty(metricCols.hit_rate) || ...
            ~isempty(metricCols.false_alarm_rate) || ~isempty(metricCols.Pr) || ~isempty(metricCols.bias_c);

        if ~isempty(subjectCol) && ~isempty(stimCol) && hasMetric
            nRows = height(T);
            out = table();
            out.subject = strings(nRows, 1);
            out.stim = strings(nRows, 1);
            out.target_group = strings(nRows, 1);

            subjVals = column_to_string(T.(char(varNames(subjectCol))));
            stimVals = column_to_string(T.(char(varNames(stimCol))));

            if ~isempty(targetCol)
                targetVals = column_to_string(T.(char(varNames(targetCol))));
            else
                targetVals = strings(nRows, 1);
            end

            for i = 1:nRows
                out.subject(i) = normalize_subject_id(subjVals(i));
                out.stim(i) = normalize_stim_label(stimVals(i));

                if ~isempty(targetCol)
                    out.target_group(i) = normalize_target_label(targetVals(i));
                else
                    subKey = char(out.subject(i));
                    if isKey(targetMap, subKey)
                        out.target_group(i) = string(targetMap(subKey));
                    else
                        out.target_group(i) = "";
                    end
                end
            end

            metricNames = fieldnames(metricCols);
            for m = 1:numel(metricNames)
                mn = metricNames{m};
                col = metricCols.(mn);
                if isempty(col)
                    out.(mn) = NaN(nRows, 1);
                else
                    rawCol = T.(char(varNames(col)));
                    out.(mn) = safe_numeric_column(rawCol);
                end
            end

            keep = strlength(out.subject) > 0 & strlength(out.stim) > 0;
            out = out(keep, :);

            if isempty(behTable)
                behTable = out;
            else
                behTable = [behTable; out]; %#ok<AGROW>
            end
        end
    end

    if targetMap.Count == 0
        error('No subject target labels were found in behavior.xlsx');
    end

    if isempty(behTable)
        warning('No behavior metric table found. behavior analyses will be empty.');
    else
        try
            behTable = unique(behTable);
        catch
        end
    end
end

function behTable = add_behavior_change_columns(behTable, offStim, activeStims)

    if isempty(behTable)
        return;
    end

    behTable.false_alarm_reduction_note = repmat("computed as OFF-active during paired analysis", height(behTable), 1);
    behTable.d_prime_improvement_note = repmat("computed as active-OFF during paired analysis", height(behTable), 1);
    behTable.hit_rate_change_note = repmat("computed as active-OFF during paired analysis", height(behTable), 1);
    behTable.Pr_improvement_note = repmat("computed as active-OFF during paired analysis", height(behTable), 1);
    behTable.bias_c_change_note = repmat("computed as active-OFF during paired analysis", height(behTable), 1);

    offStim = offStim;
    activeStims = activeStims;
end

function normNames = normalize_var_names(varNames)
    normNames = lower(string(varNames));
    normNames = regexprep(normNames, '[^a-zA-Z0-9]', '');
end

function idx = find_first_col(normNames, aliases)

    idx = [];

    for a = 1:numel(aliases)
        hit = find(normNames == lower(aliases{a}), 1);
        if ~isempty(hit)
            idx = hit;
            return;
        end
    end

    for a = 1:numel(aliases)
        hit = find(contains(normNames, lower(aliases{a})), 1);
        if ~isempty(hit)
            idx = hit;
            return;
        end
    end
end

function s = column_to_string(col)

    if isstring(col)
        s = col(:);
    elseif ischar(col)
        s = string(cellstr(col));
    elseif isnumeric(col)
        s = string(col(:));
    elseif iscategorical(col)
        s = string(col(:));
    elseif iscell(col)
        s = strings(numel(col), 1);
        for i = 1:numel(col)
            v = col{i};
            if isempty(v)
                s(i) = "";
            elseif isnumeric(v)
                s(i) = string(v(1));
            else
                s(i) = string(v);
            end
        end
    else
        s = string(col(:));
    end
end

function x = safe_numeric_column(col)

    if isnumeric(col)
        x = double(col(:));
        return;
    end

    if iscategorical(col)
        col = string(col(:));
    end

    if ischar(col)
        col = string(cellstr(col));
    end

    if isstring(col)
        x = NaN(numel(col), 1);
        for i = 1:numel(col)
            x(i) = parse_one_numeric_value(col(i));
        end
        return;
    end

    if iscell(col)
        x = NaN(numel(col), 1);
        for i = 1:numel(col)
            v = col{i};
            if isempty(v)
                x(i) = NaN;
            elseif isnumeric(v)
                x(i) = double(v(1));
            else
                x(i) = parse_one_numeric_value(string(v));
            end
        end
        return;
    end

    s = string(col(:));
    x = NaN(numel(s), 1);
    for i = 1:numel(s)
        x(i) = parse_one_numeric_value(s(i));
    end
end

function val = parse_one_numeric_value(v)

    if ismissing(v)
        val = NaN;
        return;
    end

    s = strtrim(string(v));

    if strlength(s) == 0
        val = NaN;
        return;
    end

    lowerS = lower(s);
    if lowerS == "nan" || lowerS == "na" || lowerS == "n/a" || lowerS == "missing" || lowerS == "<missing>"
        val = NaN;
        return;
    end

    hasPercent = contains(s, "%");
    s = erase(s, "%");
    s = erase(s, ",");

    token = regexp(char(s), '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');

    if isempty(token)
        val = NaN;
        return;
    end

    val = str2double(token);
    if hasPercent
        val = val / 100;
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

function tar = normalize_target_label(x)

    x = upper(strtrim(char(string(x))));

    if contains(x, 'STN')
        tar = "STN";
    elseif contains(x, 'SN')
        tar = "SN";
    else
        tar = "";
    end
end

function stim = normalize_stim_label(x)

    x = lower(strtrim(char(string(x))));

    if contains(x, '130')
        stim = "stim130";
    elseif contains(x, '10')
        stim = "stim10";
    elseif contains(x, 'off') || strcmp(x, '0') || strcmp(x, 'stim0')
        stim = "stim0";
    else
        stim = string(x);
    end
end

function T = build_session_table(rawRoot, qcLogFile, targetMap, excludeSubjects, excludeSessions)

    if exist(qcLogFile, 'file') ~= 2
        error('QC log not found: %s', qcLogFile);
    end

    q = readtable(qcLogFile, 'Sheet', 'session_log', 'VariableNamingRule', 'preserve');

    q.analysis = string(q.analysis);
    q.status = string(q.status);
    q.subject = string(q.subject);
    q.stim = string(q.stim);

    q = q(q.analysis == "Simple_Old_vs_New" & q.status == "success", :);

    subjects = strings(height(q), 1);
    stims = strings(height(q), 1);
    targetGroups = strings(height(q), 1);
    funcDirs = strings(height(q), 1);

    for i = 1:height(q)
        subjects(i) = normalize_subject_id(q.subject(i));
        stims(i) = q.stim(i);

        subKey = char(subjects(i));
        if isKey(targetMap, subKey)
            targetGroups(i) = string(targetMap(subKey));
        else
            targetGroups(i) = "";
        end

        funcDirs(i) = string(fullfile(rawRoot, char(stims(i)), char(subjects(i))));
    end

    T = table(subjects, stims, targetGroups, funcDirs, ...
        'VariableNames', {'subject','stim','target_group','func_dir'});

    exSubs = strings(numel(excludeSubjects), 1);
    for i = 1:numel(excludeSubjects)
        exSubs(i) = normalize_subject_id(excludeSubjects{i});
    end
    T = T(~ismember(T.subject, exSubs), :);

    if ~isempty(excludeSessions)
        keep = true(height(T), 1);
        for i = 1:size(excludeSessions, 1)
            exSub = normalize_subject_id(excludeSessions{i, 1});
            exStim = string(excludeSessions{i, 2});
            keep(T.subject == exSub & T.stim == exStim) = false;
        end
        T = T(keep, :);
    end

    keep = strlength(T.target_group) > 0;
    for i = 1:height(T)
        keep(i) = keep(i) && exist(char(T.func_dir(i)), 'dir') == 7;
    end

    T = T(keep, :);
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

function onsetFile = find_new_onset_file(funcDir, subj)

    onsetFile = fullfile(funcDir, sprintf('onset_%s.xlsx', subj));

    if exist(onsetFile, 'file') == 2
        return;
    end

    d = dir(fullfile(funcDir, 'onset_sub*.xlsx'));

    if isempty(d)
        onsetFile = '';
    else
        onsetFile = fullfile(d(1).folder, d(1).name);
    end
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

            rows(end+1, :) = [{char(subj), char(stim), char(tg), char(circuit), char(seedSystem), char(targetSystem), isPrimary, nEdges, status, reason}, vals]; %#ok<AGROW>
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

function statsTable = run_primary_circuit_delta_interaction_stats(circuitTable, measures, activeStims, offStim)

    rows = {};
    vars = {'family','measure','active_stim','circuit','seed_system','target_system','n_SN','n_STN','mean_delta_SN','mean_delta_STN','diff_SN_minus_STN','t','df','p','p_fdr_primary','p_perm_max'};

    if isempty(circuitTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    circuits = unique(circuitTable(circuitTable.is_primary == true, {'circuit','seed_system','target_system'}));

    for m = 1:numel(measures)
        measure = measures{m};
        for c = 1:height(circuits)
            circuit = string(circuits.circuit(c));
            seedSystem = string(circuits.seed_system(c));
            targetSystem = string(circuits.target_system(c));

            for a = 1:numel(activeStims)
                active = string(activeStims{a});

                dSN = paired_circuit_delta(circuitTable, measure, circuit, "SN", active, string(offStim));
                dSTN = paired_circuit_delta(circuitTable, measure, circuit, "STN", active, string(offStim));

                [tval, df, pval] = welch_ttest2(dSN, dSTN);

                rows(end+1, :) = {
                    'primary_delta_interaction', measure, char(active), char(circuit), char(seedSystem), char(targetSystem), ...
                    sum(isfinite(dSN)), sum(isfinite(dSTN)), ...
                    mean(dSN, 'omitnan'), mean(dSTN, 'omitnan'), mean(dSN, 'omitnan') - mean(dSTN, 'omitnan'), ...
                    tval, df, pval, NaN, NaN
                }; %#ok<AGROW>
            end
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function statsTable = run_circuit_simple_delta_stats(circuitTable, measures, activeStims, offStim, targetGroups)

    rows = {};
    vars = {'family','measure','active_stim','target_group','circuit','n','mean_delta','sd_delta','t','df','p','p_fdr_by_measure'};

    if isempty(circuitTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    circuits = unique(string(circuitTable.circuit));

    for m = 1:numel(measures)
        measure = measures{m};
        for c = 1:numel(circuits)
            circuit = circuits(c);
            for tgCell = targetGroups
                tg = string(tgCell{1});
                for a = 1:numel(activeStims)
                    active = string(activeStims{a});
                    vals = paired_circuit_delta(circuitTable, measure, circuit, tg, active, string(offStim));
                    [tval, df, pval] = one_sample_ttest(vals);

                    rows(end+1, :) = {
                        'circuit_simple_delta', measure, char(active), char(tg), char(circuit), ...
                        sum(isfinite(vals)), mean(vals, 'omitnan'), std(vals, 'omitnan'), tval, df, pval, NaN
                    }; %#ok<AGROW>
                end
            end
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function statsTable = run_circuit_frequency_stats(circuitTable, measures, targetGroups)

    rows = {};
    vars = {'family','measure','contrast','target_group','circuit','n','mean_delta','sd_delta','t','df','p','p_fdr_by_measure'};

    if isempty(circuitTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    circuits = unique(string(circuitTable.circuit));

    for m = 1:numel(measures)
        measure = measures{m};
        for c = 1:numel(circuits)
            circuit = circuits(c);
            for tgCell = targetGroups
                tg = string(tgCell{1});
                vals = paired_circuit_delta(circuitTable, measure, circuit, tg, "stim10", "stim130");
                [tval, df, pval] = one_sample_ttest(vals);

                rows(end+1, :) = {
                    'circuit_frequency', measure, 'stim10_minus_stim130', char(tg), char(circuit), ...
                    sum(isfinite(vals)), mean(vals, 'omitnan'), std(vals, 'omitnan'), tval, df, pval, NaN
                }; %#ok<AGROW>
            end
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function vals = paired_circuit_delta(circuitTable, measure, circuit, targetGroup, stimA, stimB)

    rows = circuitTable(string(circuitTable.circuit) == circuit & ...
        string(circuitTable.target_group) == targetGroup & ...
        string(circuitTable.status) == "success", :);

    subsA = unique(string(rows.subject(string(rows.stim) == stimA)));
    subsB = unique(string(rows.subject(string(rows.stim) == stimB)));
    subjects = intersect(subsA, subsB);

    vals = NaN(numel(subjects), 1);
    for i = 1:numel(subjects)
        sub = subjects(i);
        rowA = rows(string(rows.subject) == sub & string(rows.stim) == stimA, :);
        rowB = rows(string(rows.subject) == sub & string(rows.stim) == stimB, :);

        if ~isempty(rowA) && ~isempty(rowB)
            vals(i) = double(rowA.(measure)(1)) - double(rowB.(measure)(1));
        end
    end
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

function [tval, df, pval] = welch_ttest2(x, y)

    x = x(isfinite(x));
    y = y(isfinite(y));

    if numel(x) < 3 || numel(y) < 3
        tval = NaN; df = NaN; pval = NaN;
        return;
    end

    [~, pval, ~, st] = ttest2(x, y, 'Vartype', 'unequal');
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

function T = add_fdr_by_behavior_measure_family(T, pName, qName)

    if isempty(T) || ~ismember(pName, T.Properties.VariableNames)
        return;
    end

    T.(qName) = NaN(height(T), 1);
    keys = unique(strcat(string(T.behavior_metric), "___", string(T.measure)));

    for k = 1:numel(keys)
        parts = split(keys(k), "___");
        idx = string(T.behavior_metric) == parts(1) & string(T.measure) == parts(2) & isfinite(T.(pName));
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

function T = add_permutation_maxstat_delta_interaction(T, circuitTable, measures, activeStims, offStim, nPerm)

    if isempty(T)
        return;
    end

    obs = abs(T.t);
    obs(~isfinite(obs)) = NaN;
    maxStats = NaN(nPerm, 1);

    fprintf('Permutation max-stat for primary delta interaction: %d permutations\n', nPerm);

    for b = 1:nPerm
        maxT = 0;

        for r = 1:height(T)
            measure = string(T.measure(r));
            active = string(T.active_stim(r));
            circuit = string(T.circuit(r));

            dSN = paired_circuit_delta(circuitTable, char(measure), circuit, "SN", active, string(offStim));
            dSTN = paired_circuit_delta(circuitTable, char(measure), circuit, "STN", active, string(offStim));

            x = [dSN(:); dSTN(:)];
            g = [ones(numel(dSN), 1); zeros(numel(dSTN), 1)];
            good = isfinite(x) & isfinite(g);
            x = x(good);
            g = g(good);

            if numel(x) < 6 || numel(unique(g)) < 2
                continue;
            end

            gPerm = g(randperm(numel(g)));
            xpSN = x(gPerm == 1);
            xpSTN = x(gPerm == 0);
            [tval, ~, ~] = welch_ttest2(xpSN, xpSTN);

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

    measures = measures;
    activeStims = activeStims;
end

function statsTable = run_primary_behavior_interaction_stats(circuitTable, behTable, measures, activeStims, offStim, behaviorMetric)

    rows = {};
    vars = {'family','measure','active_stim','behavior_metric','circuit','n','n_SN','n_STN', ...
        'slope_STN','slope_SN','beta_interaction_SN_minus_STN','t_interaction','p_interaction', ...
        'r_SN','p_SN','r_STN','p_STN','loo_min_abs_r_SN','loo_min_abs_r_STN','p_fdr_primary','p_perm_max'};

    if isempty(circuitTable) || isempty(behTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    primaryCircuits = unique(string(circuitTable.circuit(circuitTable.is_primary == true)));

    for m = 1:numel(measures)
        measure = measures{m};
        for c = 1:numel(primaryCircuits)
            circuit = primaryCircuits(c);
            for a = 1:numel(activeStims)
                active = string(activeStims{a});

                [x, y, g] = paired_circuit_behavior_delta(circuitTable, behTable, measure, behaviorMetric, circuit, active, string(offStim));
                [lm, ok] = fit_interaction_lm(x, y, g);
                [rSN, pSN] = corr_by_group(x, y, g, 1);
                [rSTN, pSTN] = corr_by_group(x, y, g, 0);
                looSN = loo_min_abs_corr(x(g == 1), y(g == 1));
                looSTN = loo_min_abs_corr(x(g == 0), y(g == 0));

                if ok
                    slopeSTN = lm.beta(2);
                    slopeSN = lm.beta(2) + lm.beta(4);
                    betaInt = lm.beta(4);
                    tInt = lm.t(4);
                    pInt = lm.p(4);
                else
                    slopeSTN = NaN; slopeSN = NaN; betaInt = NaN; tInt = NaN; pInt = NaN;
                end

                rows(end+1, :) = {
                    'primary_behavior_interaction', measure, char(active), behaviorMetric, char(circuit), ...
                    numel(y), sum(g == 1), sum(g == 0), slopeSTN, slopeSN, betaInt, tInt, pInt, ...
                    rSN, pSN, rSTN, pSTN, looSN, looSTN, NaN, NaN
                }; %#ok<AGROW>
            end
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function statsTable = run_secondary_behavior_stats(circuitTable, behTable, measures, activeStims, offStim, behaviorMetrics)

    rows = {};
    vars = {'family','measure','active_stim','behavior_metric','circuit','n','n_SN','n_STN', ...
        'slope_STN','slope_SN','beta_interaction_SN_minus_STN','t_interaction','p_interaction','p_fdr_by_behavior_measure'};

    if isempty(circuitTable) || isempty(behTable)
        statsTable = cell2table_or_empty(rows, vars);
        return;
    end

    circuits = unique(string(circuitTable.circuit));

    for bm = 1:numel(behaviorMetrics)
        behaviorMetric = behaviorMetrics{bm};
        for m = 1:numel(measures)
            measure = measures{m};
            for c = 1:numel(circuits)
                circuit = circuits(c);
                for a = 1:numel(activeStims)
                    active = string(activeStims{a});

                    [x, y, g] = paired_circuit_behavior_delta(circuitTable, behTable, measure, behaviorMetric, circuit, active, string(offStim));
                    [lm, ok] = fit_interaction_lm(x, y, g);

                    if ok
                        slopeSTN = lm.beta(2);
                        slopeSN = lm.beta(2) + lm.beta(4);
                        betaInt = lm.beta(4);
                        tInt = lm.t(4);
                        pInt = lm.p(4);
                    else
                        slopeSTN = NaN; slopeSN = NaN; betaInt = NaN; tInt = NaN; pInt = NaN;
                    end

                    rows(end+1, :) = {
                        'secondary_behavior_interaction', measure, char(active), behaviorMetric, char(circuit), ...
                        numel(y), sum(g == 1), sum(g == 0), slopeSTN, slopeSN, betaInt, tInt, pInt, NaN
                    }; %#ok<AGROW>
                end
            end
        end
    end

    statsTable = cell2table_or_empty(rows, vars);
end

function [x, y, g] = paired_circuit_behavior_delta(circuitTable, behTable, measure, behaviorMetric, circuit, activeStim, offStim)

    rows = circuitTable(string(circuitTable.circuit) == circuit & string(circuitTable.status) == "success", :);

    subsA = unique(string(rows.subject(string(rows.stim) == activeStim)));
    subsB = unique(string(rows.subject(string(rows.stim) == offStim)));
    subjects = intersect(subsA, subsB);

    x = [];
    y = [];
    g = [];

    for i = 1:numel(subjects)
        sub = subjects(i);
        rowA = rows(string(rows.subject) == sub & string(rows.stim) == activeStim, :);
        rowB = rows(string(rows.subject) == sub & string(rows.stim) == offStim, :);

        if isempty(rowA) || isempty(rowB)
            continue;
        end

        behA = behTable(string(behTable.subject) == sub & string(behTable.stim) == activeStim, :);
        behB = behTable(string(behTable.subject) == sub & string(behTable.stim) == offStim, :);

        if isempty(behA) || isempty(behB)
            continue;
        end

        xVal = double(rowA.(measure)(1)) - double(rowB.(measure)(1));
        yVal = compute_behavior_delta(behA(1, :), behB(1, :), behaviorMetric);

        tg = string(rowA.target_group(1));
        if tg == "SN"
            gVal = 1;
        elseif tg == "STN"
            gVal = 0;
        else
            continue;
        end

        if isfinite(xVal) && isfinite(yVal)
            x(end+1, 1) = xVal; %#ok<AGROW>
            y(end+1, 1) = yVal; %#ok<AGROW>
            g(end+1, 1) = gVal; %#ok<AGROW>
        end
    end
end

function val = compute_behavior_delta(behA, behB, behaviorMetric)

    behaviorMetric = string(behaviorMetric);

    switch behaviorMetric
        case "false_alarm_reduction"
            if ismember('false_alarm_rate', behA.Properties.VariableNames)
                val = double(behB.false_alarm_rate(1)) - double(behA.false_alarm_rate(1));
            else
                val = NaN;
            end

        case "d_prime_improvement"
            if ismember('d_prime', behA.Properties.VariableNames)
                val = double(behA.d_prime(1)) - double(behB.d_prime(1));
            else
                val = NaN;
            end

        case "hit_rate_change"
            if ismember('hit_rate', behA.Properties.VariableNames)
                val = double(behA.hit_rate(1)) - double(behB.hit_rate(1));
            else
                val = NaN;
            end

        case "Pr_improvement"
            if ismember('Pr', behA.Properties.VariableNames)
                val = double(behA.Pr(1)) - double(behB.Pr(1));
            else
                val = NaN;
            end

        case "bias_c_change"
            if ismember('bias_c', behA.Properties.VariableNames)
                val = double(behA.bias_c(1)) - double(behB.bias_c(1));
            else
                val = NaN;
            end

        otherwise
            val = NaN;
    end
end

function [lm, ok] = fit_interaction_lm(x, y, g)

    good = isfinite(x) & isfinite(y) & isfinite(g);
    x = x(good);
    y = y(good);
    g = g(good);

    lm = struct('beta', NaN(4,1), 'se', NaN(4,1), 't', NaN(4,1), 'p', NaN(4,1), 'df', NaN);
    ok = false;

    if numel(y) < 8 || numel(unique(g)) < 2 || sum(g == 1) < 3 || sum(g == 0) < 3
        return;
    end

    xz = safe_zscore(x);
    X = [ones(numel(y), 1), xz, g(:), xz .* g(:)];

    if rank(X) < size(X, 2)
        return;
    end

    beta = X \ y;
    resid = y - X * beta;
    n = numel(y);
    p = size(X, 2);
    df = n - p;

    if df <= 0
        return;
    end

    sigma2 = sum(resid.^2) / df;
    covb = sigma2 * inv(X' * X);
    se = sqrt(diag(covb));
    tvals = beta ./ se;
    pvals = 2 * tcdf(-abs(tvals), df);

    lm.beta = beta;
    lm.se = se;
    lm.t = tvals;
    lm.p = pvals;
    lm.df = df;
    ok = true;
end

function [r, p] = corr_by_group(x, y, g, groupValue)

    idx = isfinite(x) & isfinite(y) & g == groupValue;

    if sum(idx) < 5
        r = NaN; p = NaN;
        return;
    end

    [r, p] = corr(x(idx), y(idx), 'Type', 'Spearman', 'Rows', 'complete');
end

function v = loo_min_abs_corr(x, y)

    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    if numel(x) < 6
        v = NaN;
        return;
    end

    rs = NaN(numel(x), 1);
    for i = 1:numel(x)
        keep = true(numel(x), 1);
        keep(i) = false;
        rs(i) = corr(x(keep), y(keep), 'Type', 'Spearman', 'Rows', 'complete');
    end

    v = min(abs(rs), [], 'omitnan');
end

function T = add_permutation_maxstat_behavior_interaction(T, circuitTable, behTable, measures, activeStims, offStim, behaviorMetric, nPerm)

    if isempty(T)
        return;
    end

    obs = abs(T.t_interaction);
    obs(~isfinite(obs)) = NaN;
    maxStats = NaN(nPerm, 1);

    fprintf('Permutation max-stat for primary behavior interaction: %d permutations\n', nPerm);

    data = cell(height(T), 1);
    for r = 1:height(T)
        [x, y, g] = paired_circuit_behavior_delta(circuitTable, behTable, char(T.measure(r)), behaviorMetric, string(T.circuit(r)), string(T.active_stim(r)), string(offStim));
        data{r} = struct('x', x, 'y', y, 'g', g);
    end

    for b = 1:nPerm
        maxT = 0;

        for r = 1:height(T)
            x = data{r}.x;
            y = data{r}.y;
            g = data{r}.g;

            if numel(y) < 8
                continue;
            end

            yPerm = y;
            for gv = [0 1]
                idx = find(g == gv);
                if numel(idx) > 1
                    yPerm(idx) = y(idx(randperm(numel(idx))));
                end
            end

            [lm, ok] = fit_interaction_lm(x, yPerm, g);
            if ok && isfinite(lm.t(4))
                maxT = max(maxT, abs(lm.t(4)));
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

    measures = measures;
    activeStims = activeStims;
end

function edgeStats = run_edge_followup_delta(edgeTable, circuitDefTable, seedMembership, targetMembership, measures, activeStims, offStim)

    rows = {};
    vars = {'family','measure','active_stim','circuit','seed','target_roi','n_SN','n_STN','mean_delta_SN','mean_delta_STN','diff_SN_minus_STN','t','df','p','p_fdr_within_circuit_measure'};

    if isempty(edgeTable) || isempty(circuitDefTable)
        edgeStats = cell2table_or_empty(rows, vars);
        return;
    end

    for c = 1:height(circuitDefTable)
        circuit = string(circuitDefTable.circuit(c));
        seedNames = system_members(seedMembership, string(circuitDefTable.seed_system(c)));
        targetNames = system_members(targetMembership, string(circuitDefTable.target_system(c)));

        for m = 1:numel(measures)
            measure = measures{m};
            for si = 1:numel(seedNames)
                for ti = 1:numel(targetNames)
                    seed = seedNames(si);
                    targetROI = targetNames(ti);

                    for a = 1:numel(activeStims)
                        active = string(activeStims{a});
                        dSN = paired_edge_delta(edgeTable, measure, seed, targetROI, "SN", active, string(offStim));
                        dSTN = paired_edge_delta(edgeTable, measure, seed, targetROI, "STN", active, string(offStim));
                        [tval, df, pval] = welch_ttest2(dSN, dSTN);

                        rows(end+1, :) = {
                            'edge_followup_delta', measure, char(active), char(circuit), char(seed), char(targetROI), ...
                            sum(isfinite(dSN)), sum(isfinite(dSTN)), mean(dSN, 'omitnan'), mean(dSTN, 'omitnan'), ...
                            mean(dSN, 'omitnan') - mean(dSTN, 'omitnan'), tval, df, pval, NaN
                        }; %#ok<AGROW>
                    end
                end
            end
        end
    end

    edgeStats = cell2table_or_empty(rows, vars);
end

function vals = paired_edge_delta(edgeTable, measure, seed, targetROI, targetGroup, stimA, stimB)

    rows = edgeTable(string(edgeTable.seed) == seed & ...
        string(edgeTable.target_roi) == targetROI & ...
        string(edgeTable.target_group) == targetGroup & ...
        string(edgeTable.status) == "success", :);

    subsA = unique(string(rows.subject(string(rows.stim) == stimA)));
    subsB = unique(string(rows.subject(string(rows.stim) == stimB)));
    subjects = intersect(subsA, subsB);

    vals = NaN(numel(subjects), 1);
    for i = 1:numel(subjects)
        sub = subjects(i);
        rowA = rows(string(rows.subject) == sub & string(rows.stim) == stimA, :);
        rowB = rows(string(rows.subject) == sub & string(rows.stim) == stimB, :);

        if ~isempty(rowA) && ~isempty(rowB)
            vals(i) = double(rowA.(measure)(1)) - double(rowB.(measure)(1));
        end
    end
end

function edgeStats = run_edge_followup_behavior(edgeTable, behTable, circuitDefTable, seedMembership, targetMembership, measures, activeStims, offStim, behaviorMetric)

    rows = {};
    vars = {'family','measure','active_stim','behavior_metric','circuit','seed','target_roi','n','n_SN','n_STN', ...
        'slope_STN','slope_SN','beta_interaction_SN_minus_STN','t_interaction','p_interaction','p_fdr_within_circuit_measure'};

    if isempty(edgeTable) || isempty(behTable) || isempty(circuitDefTable)
        edgeStats = cell2table_or_empty(rows, vars);
        return;
    end

    for c = 1:height(circuitDefTable)
        circuit = string(circuitDefTable.circuit(c));
        seedNames = system_members(seedMembership, string(circuitDefTable.seed_system(c)));
        targetNames = system_members(targetMembership, string(circuitDefTable.target_system(c)));

        for m = 1:numel(measures)
            measure = measures{m};
            for si = 1:numel(seedNames)
                for ti = 1:numel(targetNames)
                    seed = seedNames(si);
                    targetROI = targetNames(ti);

                    for a = 1:numel(activeStims)
                        active = string(activeStims{a});
                        [x, y, g] = paired_edge_behavior_delta(edgeTable, behTable, measure, behaviorMetric, seed, targetROI, active, string(offStim));
                        [lm, ok] = fit_interaction_lm(x, y, g);

                        if ok
                            slopeSTN = lm.beta(2);
                            slopeSN = lm.beta(2) + lm.beta(4);
                            betaInt = lm.beta(4);
                            tInt = lm.t(4);
                            pInt = lm.p(4);
                        else
                            slopeSTN = NaN; slopeSN = NaN; betaInt = NaN; tInt = NaN; pInt = NaN;
                        end

                        rows(end+1, :) = {
                            'edge_followup_behavior', measure, char(active), behaviorMetric, char(circuit), char(seed), char(targetROI), ...
                            numel(y), sum(g == 1), sum(g == 0), slopeSTN, slopeSN, betaInt, tInt, pInt, NaN
                        }; %#ok<AGROW>
                    end
                end
            end
        end
    end

    edgeStats = cell2table_or_empty(rows, vars);
end

function [x, y, g] = paired_edge_behavior_delta(edgeTable, behTable, measure, behaviorMetric, seed, targetROI, activeStim, offStim)

    rows = edgeTable(string(edgeTable.seed) == seed & ...
        string(edgeTable.target_roi) == targetROI & ...
        string(edgeTable.status) == "success", :);

    subsA = unique(string(rows.subject(string(rows.stim) == activeStim)));
    subsB = unique(string(rows.subject(string(rows.stim) == offStim)));
    subjects = intersect(subsA, subsB);

    x = [];
    y = [];
    g = [];

    for i = 1:numel(subjects)
        sub = subjects(i);
        rowA = rows(string(rows.subject) == sub & string(rows.stim) == activeStim, :);
        rowB = rows(string(rows.subject) == sub & string(rows.stim) == offStim, :);

        if isempty(rowA) || isempty(rowB)
            continue;
        end

        behA = behTable(string(behTable.subject) == sub & string(behTable.stim) == activeStim, :);
        behB = behTable(string(behTable.subject) == sub & string(behTable.stim) == offStim, :);

        if isempty(behA) || isempty(behB)
            continue;
        end

        xVal = double(rowA.(measure)(1)) - double(rowB.(measure)(1));
        yVal = compute_behavior_delta(behA(1, :), behB(1, :), behaviorMetric);

        tg = string(rowA.target_group(1));
        if tg == "SN"
            gVal = 1;
        elseif tg == "STN"
            gVal = 0;
        else
            continue;
        end

        if isfinite(xVal) && isfinite(yVal)
            x(end+1, 1) = xVal; %#ok<AGROW>
            y(end+1, 1) = yVal; %#ok<AGROW>
            g(end+1, 1) = gVal; %#ok<AGROW>
        end
    end
end

function pathValue = request_existing_directory(promptText)
    pathValue = clean_user_path(input(promptText, 's'));
    if isempty(pathValue) || exist(pathValue, 'dir') ~= 7
        error('Directory not found: %s', pathValue);
    end
end

function pathValue = request_existing_file(promptText)
    pathValue = clean_user_path(input(promptText, 's'));
    if isempty(pathValue) || exist(pathValue, 'file') ~= 2
        error('File not found: %s', pathValue);
    end
end

function pathValue = request_output_directory(promptText)
    pathValue = clean_user_path(input(promptText, 's'));
    if isempty(pathValue)
        error('An output directory is required.');
    end
    if exist(pathValue, 'dir') ~= 7
        [ok, message] = mkdir(pathValue);
        if ~ok
            error('Cannot create output directory: %s', message);
        end
    end
end

function pathValue = clean_user_path(pathValue)
    pathValue = strtrim(pathValue);
    if numel(pathValue) >= 2
        quotedDouble = pathValue(1) == '"' && pathValue(end) == '"';
        quotedSingle = pathValue(1) == char(39) && pathValue(end) == char(39);
        if quotedDouble || quotedSingle
            pathValue = pathValue(2:end-1);
        end
    end
end
