%% Healthy-control New-versus-Old contrasts
% The first-level GLM must already be specified and estimated in SPM12.
% Expected conditions: Old_valid and New_valid.

clear; clc;

SPM_DIR = select_directory('Select the SPM12 directory');
ROOT_DIR = select_directory('Select the HC first-level root directory');

FIRSTLEVEL_FOLDER = 'GLM_NewOld';
SUBJECTS = {};
OVERWRITE_EXISTING_CONTRASTS = true;
OLD_CONDITION_ALIASES = {'Old_valid', 'OldValid'};
NEW_CONDITION_ALIASES = {'New_valid', 'NewValid'};

addpath(SPM_DIR);
assert(exist('spm', 'file') == 2, 'SPM12 was not found in the selected directory.');
spm('Defaults', 'fMRI');
spm_get_defaults('cmdline', true);

if isempty(SUBJECTS)
    d = dir(fullfile(ROOT_DIR, 'sub*'));
    d = d([d.isdir]);
    SUBJECTS = {d.name};
end
SUBJECTS = sort_subject_folders(SUBJECTS);

report = struct('subject', {}, 'status', {}, 'firstlevel_dir', {}, ...
    'n_columns', {}, 'old_columns', {}, 'new_columns', {}, 'message', {});

for iSub = 1:numel(SUBJECTS)
    subName = SUBJECTS{iSub};
    subDir = fullfile(ROOT_DIR, subName);

    report(end+1).subject = subName; %#ok<SAGROW>
    report(end).status = 'NOT_RUN';
    report(end).firstlevel_dir = '';
    report(end).n_columns = NaN;
    report(end).old_columns = '';
    report(end).new_columns = '';
    report(end).message = '';

    try
        spmMat = locate_firstlevel_spm(subDir, FIRSTLEVEL_FOLDER);
        [SPM, oldIdx, newIdx] = create_new_old_contrasts( ...
            spmMat, OLD_CONDITION_ALIASES, NEW_CONDITION_ALIASES, ...
            OVERWRITE_EXISTING_CONTRASTS);

        report(end).status = 'OK';
        report(end).firstlevel_dir = fileparts(spmMat);
        report(end).n_columns = numel(SPM.xX.name);
        report(end).old_columns = strtrim(sprintf('%d ', oldIdx));
        report(end).new_columns = strtrim(sprintf('%d ', newIdx));
        report(end).message = 'Contrasts created';

        fprintf('[OK] %s\n', subName);
    catch ME
        report(end).status = 'FAILED';
        report(end).message = ME.message;
        fprintf('[FAILED] %s: %s\n', subName, ME.message);

        try
            spmMat = locate_firstlevel_spm(subDir, FIRSTLEVEL_FOLDER);
            dump_design_names(spmMat);
        catch
        end
    end
end

outCsv = fullfile(ROOT_DIR, 'hc_new_old_contrast_report.csv');
writetable(struct2table(report), outCsv);
fprintf('Report saved to: %s\n', outCsv);

function spmMat = locate_firstlevel_spm(subDir, firstlevelFolder)
    candidates = {
        fullfile(subDir, firstlevelFolder, 'SPM.mat')
        fullfile(subDir, 'SPM.mat')
    };

    spmMat = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            spmMat = candidates{i};
            return;
        end
    end

    error('SPM.mat was not found for %s.', subDir);
end

function [SPM, oldIdx, newIdx] = create_new_old_contrasts( ...
    spmMat, oldAliases, newAliases, overwriteExisting)

    load(spmMat, 'SPM');
    firstDir = fileparts(spmMat);
    SPM.swd = firstDir;

    names = SPM.xX.name(:);
    nCols = numel(names);
    oldIdx = find_condition_columns(names, oldAliases);
    newIdx = find_condition_columns(names, newAliases);

    if isempty(oldIdx) || isempty(newIdx)
        dump_design_names(spmMat);
        error('Old_valid and/or New_valid regressors were not identified.');
    end

    if overwriteExisting
        delete_existing_contrast_files(firstDir);
    end

    wNewGtOld = zeros(1, nCols);
    wOldGtNew = zeros(1, nCols);
    wNewGtOld(newIdx) = 1 / numel(newIdx);
    wNewGtOld(oldIdx) = -1 / numel(oldIdx);
    wOldGtNew(newIdx) = -1 / numel(newIdx);
    wOldGtNew(oldIdx) = 1 / numel(oldIdx);

    xCon = spm_FcUtil('Set', 'New_valid > Old_valid', 'T', 'c', ...
        wNewGtOld', SPM.xX.xKXs);
    xCon(2) = spm_FcUtil('Set', 'Old_valid > New_valid', 'T', 'c', ...
        wOldGtNew', SPM.xX.xKXs);
    SPM.xCon = xCon;

    save(spmMat, 'SPM');

    oldPwd = pwd;
    cleanupObj = onCleanup(@() cd(oldPwd)); %#ok<NASGU>
    cd(firstDir);
    SPM = spm_contrasts(SPM, 1:numel(SPM.xCon));
    save(spmMat, 'SPM');

    assert(contrast_file_exists(firstDir, 1), 'con_0001 was not created.');
    assert(contrast_file_exists(firstDir, 2), 'con_0002 was not created.');
end

function idx = find_condition_columns(names, aliases)
    aliasNorm = cellfun(@normalize_label, aliases, 'UniformOutput', false);
    idx = [];

    for i = 1:numel(names)
        name = names{i};
        if isempty(regexp(name, '\*bf\(1\)', 'once'))
            continue;
        end
        if ~isempty(regexp(name, '\^\d+\*bf\(1\)', 'once'))
            continue;
        end

        condition = extract_condition_label(name);
        if any(strcmp(normalize_label(condition), aliasNorm))
            idx(end+1) = i; %#ok<AGROW>
        end
    end
end

function condition = extract_condition_label(name)
    token = regexp(name, 'Sn\(\d+\)\s+(.+?)\*bf\(1\)', 'tokens', 'once');
    if isempty(token)
        token = regexp(name, '(.+?)\*bf\(1\)', 'tokens', 'once');
    end

    if isempty(token)
        condition = name;
    else
        condition = strtrim(token{1});
    end
end

function value = normalize_label(value)
    value = lower(char(value));
    value = regexprep(value, '[^a-z0-9]', '');
end

function delete_existing_contrast_files(firstDir)
    patterns = {'con_*.nii', 'con_*.img', 'con_*.hdr', ...
        'spmT_*.nii', 'spmT_*.img', 'spmT_*.hdr'};

    for i = 1:numel(patterns)
        files = dir(fullfile(firstDir, patterns{i}));
        for j = 1:numel(files)
            delete(fullfile(files(j).folder, files(j).name));
        end
    end
end

function tf = contrast_file_exists(firstDir, index)
    base = sprintf('con_%04d', index);
    tf = exist(fullfile(firstDir, [base '.nii']), 'file') == 2 || ...
        exist(fullfile(firstDir, [base '.img']), 'file') == 2;
end

function dump_design_names(spmMat)
    load(spmMat, 'SPM');
    outFile = fullfile(fileparts(spmMat), 'design_regressor_names.txt');
    fid = fopen(outFile, 'w');
    if fid < 0
        return;
    end

    cleanupObj = onCleanup(@() fclose(fid)); 
    for i = 1:numel(SPM.xX.name)
        fprintf(fid, '%03d\t%s\n', i, SPM.xX.name{i});
    end
end

function sorted = sort_subject_folders(subjects)
    numbers = nan(size(subjects));
    for i = 1:numel(subjects)
        token = regexp(subjects{i}, '\d+', 'match', 'once');
        if ~isempty(token)
            numbers(i) = str2double(token);
        end
    end
    [~, order] = sortrows([numbers(:), (1:numel(subjects))']);
    sorted = subjects(order);
end

function folder = select_directory(prompt)
    folder = uigetdir(pwd, prompt);
    if isequal(folder, 0)
        error('Directory selection was cancelled.');
    end
end
