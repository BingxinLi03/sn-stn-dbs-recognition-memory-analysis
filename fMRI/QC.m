%% firstlevel_oldnew_qc_reviewer.m
% Post-analysis quality control for firstlevel_oldnew_reviewer.m.
% This script checks output completeness, manuscript-defined motion limits,
% explicit-mask coverage, and first-level contrast-image integrity.

clear; clc;

SPM_DIR = prompt_existing_dir('Enter the full path to the SPM12 directory: ');
RAW_ROOT = prompt_existing_dir('Enter the full path to the first-level RawData directory: ');
RESULT_ROOT = prompt_existing_dir('Enter the full path to the first-level output directory: ');

addpath(SPM_DIR);
spm('defaults', 'FMRI');
spm_jobman('initcfg');

ANALYSIS_NAME = 'Simple_Old_vs_New';
LOG_FILE = fullfile(RESULT_ROOT, 'firstlevel_oldnew_QC_log.xlsx');
OUTPUT_FILE = fullfile(RESULT_ROOT, 'firstlevel_oldnew_postQC.xlsx');

assert(exist(LOG_FILE, 'file') == 2, 'First-level QC log not found: %s', LOG_FILE);

MAX_TRANSLATION_MM = 2;
MAX_ROTATION_DEG = 2;
MASK_LOW_RELATIVE_TO_MEDIAN = 0.50;

T = readtable(LOG_FILE, 'Sheet', 'session_log', 'VariableNamingRule', 'preserve');
T.analysis = string(T.analysis);
T.stim = string(T.stim);
T.subject = string(T.subject);
T.status = string(T.status);
T.reason = string(T.reason);
T.output_dir = string(T.output_dir);

usableStatus = T.status == "success" | T.status == "skipped_existing";
usable = T(T.analysis == ANALYSIS_NAME & usableStatus, :);

summaryRows = {
    'total_log_rows', height(T);
    'usable_sessions', height(usable);
    'success', sum(T.status == "success");
    'skipped_existing', sum(T.status == "skipped_existing");
    'skipped_QC', sum(T.status == "skipped_QC");
    'failed', sum(T.status == "failed")
};
summaryTable = cell2table(summaryRows, 'VariableNames', {'item', 'n'});

integrityRows = {};
motionRows = {};
maskRows = {};
contrastRows = {};

for r = 1:height(usable)
    stim = char(usable.stim(r));
    subject = char(usable.subject(r));
    outDir = char(usable.output_dir(r));
    funcDir = fullfile(RAW_ROOT, stim, subject);

    spmmat = fullfile(outDir, 'SPM.mat');
    con1 = fullfile(outDir, 'con_0001.nii');
    con2 = fullfile(outDir, 'con_0002.nii');
    con3 = fullfile(outDir, 'con_0003.nii');
    spmT1 = fullfile(outDir, 'spmT_0001.nii');
    spmT2 = fullfile(outDir, 'spmT_0002.nii');
    spmT3 = fullfile(outDir, 'spmT_0003.nii');
    resms = fullfile(outDir, 'ResMS.nii');
    implicitMask = fullfile(outDir, 'mask.nii');
    explicitMask = fullfile(outDir, 'explicit_mask_GM_AND_ARTIFACT.nii');

    existsSPM = exist(spmmat, 'file') == 2;
    existsCon1 = exist(con1, 'file') == 2;
    existsCon2 = exist(con2, 'file') == 2;
    existsCon3 = exist(con3, 'file') == 2;
    existsT1 = exist(spmT1, 'file') == 2;
    existsT2 = exist(spmT2, 'file') == 2;
    existsT3 = exist(spmT3, 'file') == 2;
    existsResMS = exist(resms, 'file') == 2;
    existsImplicitMask = exist(implicitMask, 'file') == 2;
    existsExplicitMask = exist(explicitMask, 'file') == 2;
    outputsComplete = existsSPM && existsCon1 && existsCon2 && existsCon3 && ...
        existsT1 && existsT2 && existsT3 && existsResMS;

    integrityRows(end+1, :) = { ...
        stim, subject, outDir, existsSPM, existsCon1, existsCon2, existsCon3, ...
        existsT1, existsT2, existsT3, existsResMS, existsImplicitMask, ...
        existsExplicitMask, outputsComplete}; %#ok<AGROW>

    rpFile = find_latest_rp_file(funcDir);
    if isempty(rpFile)
        motionRows(end+1, :) = {stim, subject, '', NaN, NaN, true, 'Motion file missing'}; %#ok<AGROW>
        motionExceeds = true;
    else
        try
            motion = compute_motion_qc(rpFile);
            motionExceeds = motion.maxTranslationMM > MAX_TRANSLATION_MM || ...
                motion.maxRotationDeg > MAX_ROTATION_DEG;
            reason = '';
            if motion.maxTranslationMM > MAX_TRANSLATION_MM
                reason = append_reason(reason, sprintf('Translation > %.1f mm.', MAX_TRANSLATION_MM));
            end
            if motion.maxRotationDeg > MAX_ROTATION_DEG
                reason = append_reason(reason, sprintf('Rotation > %.1f degrees.', MAX_ROTATION_DEG));
            end
            motionRows(end+1, :) = {stim, subject, rpFile, ...
                motion.maxTranslationMM, motion.maxRotationDeg, motionExceeds, reason}; %#ok<AGROW>
        catch ME
            motionRows(end+1, :) = {stim, subject, rpFile, NaN, NaN, true, ME.message}; %#ok<AGROW>
            motionExceeds = true;
        end
    end

    selectedMask = '';
    if existsExplicitMask
        selectedMask = explicitMask;
    elseif existsImplicitMask
        selectedMask = implicitMask;
    end

    if isempty(selectedMask)
        maskRows(end+1, :) = {stim, subject, '', NaN, NaN, true, 'Mask missing'}; %#ok<AGROW>
    else
        try
            [nVoxels, volumeMl] = compute_mask_size(selectedMask);
            maskRows(end+1, :) = {stim, subject, selectedMask, nVoxels, volumeMl, false, ''}; %#ok<AGROW>
        catch ME
            maskRows(end+1, :) = {stim, subject, selectedMask, NaN, NaN, true, ME.message}; %#ok<AGROW>
        end
    end

    contrastFiles = {con1, con2, con3};
    contrastNames = {'OldValid_gt_NewValid', 'NewValid_gt_OldValid', 'AllValidRetrieval_gt_Baseline'};
    for c = 1:numel(contrastFiles)
        if exist(contrastFiles{c}, 'file') ~= 2
            contrastRows(end+1, :) = {stim, subject, contrastNames{c}, contrastFiles{c}, ...
                NaN, NaN, NaN, NaN, 0, 0, true, 'Contrast image missing'}; %#ok<AGROW>
            continue;
        end
        try
            stats = compute_image_stats(contrastFiles{c}, selectedMask);
            warningFlag = stats.nFinite == 0 || stats.nNonZero == 0;
            reason = '';
            if stats.nFinite == 0
                reason = 'No finite voxels in the selected mask.';
            elseif stats.nNonZero == 0
                reason = 'All voxels are zero in the selected mask.';
            end
            contrastRows(end+1, :) = {stim, subject, contrastNames{c}, contrastFiles{c}, ...
                stats.meanValue, stats.stdValue, stats.minValue, stats.maxValue, ...
                stats.nFinite, stats.nNonZero, warningFlag, reason}; %#ok<AGROW>
        catch ME
            contrastRows(end+1, :) = {stim, subject, contrastNames{c}, contrastFiles{c}, ...
                NaN, NaN, NaN, NaN, 0, 0, true, ME.message}; %#ok<AGROW>
        end
    end
end

integrityTable = cell2table(integrityRows, 'VariableNames', {
    'stim', 'subject', 'output_dir', 'SPM_mat', 'con_0001', 'con_0002', 'con_0003', ...
    'spmT_0001', 'spmT_0002', 'spmT_0003', 'ResMS', 'mask_nii', 'explicit_mask', ...
    'outputs_complete'
});

motionTable = cell2table(motionRows, 'VariableNames', {
    'stim', 'subject', 'rp_file', 'max_translation_mm', 'max_rotation_deg', ...
    'exceeds_manuscript_motion_limit', 'reason'
});

maskTable = cell2table(maskRows, 'VariableNames', {
    'stim', 'subject', 'mask_file', 'n_mask_voxels', 'mask_volume_ml', ...
    'mask_warning', 'reason'
});

contrastTable = cell2table(contrastRows, 'VariableNames', {
    'stim', 'subject', 'contrast', 'contrast_file', 'mean_value', 'std_value', ...
    'min_value', 'max_value', 'n_finite_voxels', 'n_nonzero_voxels', ...
    'contrast_warning', 'reason'
});

if ~isempty(maskTable)
    validMask = ~isnan(maskTable.n_mask_voxels);
    medianMask = median(maskTable.n_mask_voxels(validMask), 'omitnan');
    lowMask = validMask & maskTable.n_mask_voxels < medianMask * MASK_LOW_RELATIVE_TO_MEDIAN;
    maskTable.mask_warning(lowMask) = true;
    for i = find(lowMask)'
        maskTable.reason{i} = append_reason(maskTable.reason{i}, ...
            sprintf('Mask contains fewer than %.0f%% of the median voxel count.', ...
            MASK_LOW_RELATIVE_TO_MEDIAN * 100));
    end
end

flaggedRows = outerjoin(integrityTable(:, {'stim', 'subject', 'outputs_complete'}), ...
    motionTable(:, {'stim', 'subject', 'exceeds_manuscript_motion_limit'}), ...
    'Keys', {'stim', 'subject'}, 'MergeKeys', true, 'Type', 'full');
flaggedRows = outerjoin(flaggedRows, ...
    maskTable(:, {'stim', 'subject', 'mask_warning'}), ...
    'Keys', {'stim', 'subject'}, 'MergeKeys', true, 'Type', 'full');
flaggedRows.flagged = ~flaggedRows.outputs_complete | ...
    flaggedRows.exceeds_manuscript_motion_limit | flaggedRows.mask_warning;
flaggedRows = flaggedRows(flaggedRows.flagged, :);

if exist(OUTPUT_FILE, 'file') == 2
    delete(OUTPUT_FILE);
end
writetable(summaryTable, OUTPUT_FILE, 'Sheet', 'summary');
writetable(T, OUTPUT_FILE, 'Sheet', 'firstlevel_log');
writetable(integrityTable, OUTPUT_FILE, 'Sheet', 'output_integrity');
writetable(motionTable, OUTPUT_FILE, 'Sheet', 'motion_qc');
writetable(maskTable, OUTPUT_FILE, 'Sheet', 'mask_coverage');
writetable(contrastTable, OUTPUT_FILE, 'Sheet', 'contrast_integrity');
writetable(flaggedRows, OUTPUT_FILE, 'Sheet', 'flagged_sessions');

notes = {
    'Motion exclusion rule', sprintf('Translation > %.1f mm or rotation > %.1f degrees in any direction', MAX_TRANSLATION_MM, MAX_ROTATION_DEG);
    'Mask warning rule', sprintf('Mask voxel count < %.0f%% of the median across usable sessions', MASK_LOW_RELATIVE_TO_MEDIAN * 100);
    'Formal second-level analysis', 'Not performed by this script';
};
writetable(cell2table(notes, 'VariableNames', {'item', 'value'}), OUTPUT_FILE, 'Sheet', 'notes');

fprintf('\nPost-analysis QC completed.\nReport: %s\n', OUTPUT_FILE);

function pathOut = prompt_existing_dir(message)
    pathOut = clean_path_input(input(message, 's'));
    assert(exist(pathOut, 'dir') == 7, 'Directory not found: %s', pathOut);
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

function rpFile = find_latest_rp_file(funcDir)
    rpFile = '';
    files = dir(fullfile(funcDir, 'rp*.txt'));
    if isempty(files)
        return;
    end
    [~, idx] = max([files.datenum]);
    rpFile = fullfile(files(idx).folder, files(idx).name);
end

function motion = compute_motion_qc(rpFile)
    parameters = load(rpFile);
    assert(size(parameters, 2) >= 6, 'Motion file has fewer than six columns.');
    translation = parameters(:, 1:3);
    rotationDeg = parameters(:, 4:6) * 180 / pi;
    motion.maxTranslationMM = max(abs(translation(:)));
    motion.maxRotationDeg = max(abs(rotationDeg(:)));
end

function [nVoxels, volumeMl] = compute_mask_size(maskFile)
    V = spm_vol(maskFile);
    Y = spm_read_vols(V);
    nVoxels = sum(Y(:) > 0);
    voxelVolumeMm3 = abs(det(V.mat(1:3, 1:3)));
    volumeMl = nVoxels * voxelVolumeMm3 / 1000;
end

function stats = compute_image_stats(imageFile, maskFile)
    V = spm_vol(imageFile);
    Y = spm_read_vols(V);

    if ~isempty(maskFile) && exist(maskFile, 'file') == 2
        Vm = spm_vol(maskFile);
        M = spm_read_vols(Vm) > 0;
        assert(isequal(size(Y), size(M)), 'Image and mask dimensions differ.');
    else
        M = true(size(Y));
    end

    values = Y(M & isfinite(Y));
    if isempty(values)
        stats.meanValue = NaN;
        stats.stdValue = NaN;
        stats.minValue = NaN;
        stats.maxValue = NaN;
        stats.nFinite = 0;
        stats.nNonZero = 0;
        return;
    end

    stats.meanValue = mean(values);
    stats.stdValue = std(values);
    stats.minValue = min(values);
    stats.maxValue = max(values);
    stats.nFinite = numel(values);
    stats.nNonZero = sum(values ~= 0);
end

function reason = append_reason(reason, addition)
    if isempty(reason)
        reason = addition;
    else
        reason = [reason ' ' addition];
    end
end
