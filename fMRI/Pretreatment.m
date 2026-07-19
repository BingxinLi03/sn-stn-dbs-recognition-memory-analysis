%% DBS task-fMRI preprocessing
% Input: slice-time-corrected functional NIfTI files matching A*.nii.
% Pipeline: realignment, T1 segmentation, artifact-mask preparation,
% optional enantiomorphic filling for coregistration estimation,
% EPI-to-T1 coregistration, normalization, 6-mm smoothing, and GLM mask creation.

clear; clc;

STIMS = {'stim0','stim10','stim130'};
FUNC_FOLDER = 'FunImgHA';
T1_STIM = 'stim0';
T1_FOLDER = 'T1ImgH';
MASK_FOLDER = 'mask';

RUN_STAGE = 2;                  % 1: prepare files for manual masks; 2: preprocess
FORCE_REALIGN = false;
USE_ENANTIOMORPHIC_FILL = true;
USE_COREG_WEIGHTING = true;    % Requires an SPM build that accepts eoptions.weight
CLEAR_WORKDIR_EACH_RUN = true;
RESTORE_R_HEADERS = true;

NORM_VOX = [3 3 3];
SMOOTH_FWHM = [6 6 6];
BB = NaN(2,3);

SPM_DIR = request_directory('Enter the full path to the SPM12 directory: ');
ROOT_DIR = request_directory('Enter the full path to the preprocessing data root: ');

CODE_DIR = '';
if USE_ENANTIOMORPHIC_FILL
    CODE_DIR = request_directory(...
        'Enter the directory containing entiamorphicSub1.m: ');
end

assert(exist(SPM_DIR,'dir') == 7, 'SPM directory not found: %s', SPM_DIR);
assert(exist(ROOT_DIR,'dir') == 7, 'Data root not found: %s', ROOT_DIR);

addpath(SPM_DIR);
if ~isempty(CODE_DIR)
    addpath(CODE_DIR);
    assert(exist('entiamorphicSub1','file') == 2, ...
        'entiamorphicSub1.m was not found in: %s', CODE_DIR);
end

spm('defaults','fmri');
spm_jobman('initcfg');

subRoot = fullfile(ROOT_DIR, STIMS{1}, FUNC_FOLDER);
assert(exist(subRoot,'dir') == 7, ...
    'Functional-data directory not found: %s', subRoot);

subDirs = dir(fullfile(subRoot,'sub*'));
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name},{'.','..'}));
assert(~isempty(subDirs), 'No sub* directories were found under: %s', subRoot);

fprintf('Found %d subjects under %s\n', numel(subDirs), subRoot);

for i = 1:numel(subDirs)
    subName = subDirs(i).name;
    subNum = parse_sub_number(subName);
    if isempty(subNum)
        fprintf('[SKIP] Could not parse subject number: %s\n', subName);
        continue;
    end

    fprintf('\n================ %s ================\n', subName);

    t1Dir = fullfile(ROOT_DIR, T1_STIM, T1_FOLDER, subName);
    t1File = find_t1_file_robust(t1Dir);

    if RUN_STAGE == 1
        for s = 1:numel(STIMS)
            stim = STIMS{s};
            funcDir = fullfile(ROOT_DIR, stim, FUNC_FOLDER, subName);
            if exist(funcDir,'dir') ~= 7
                fprintf('[MISS] %s | directory not found: %s\n', subName, funcDir);
                continue;
            end

            fprintf('--- Stage 1 | %s | %s ---\n', subName, stim);
            [meanFile, ~] = ensure_realign(funcDir, FORCE_REALIGN);

            maskDir = fullfile(ROOT_DIR, stim, MASK_FOLDER, subName);
            if exist(maskDir,'dir') ~= 7
                mkdir(maskDir);
            end

            meanCopyForMask = fullfile(maskDir, ...
                sprintf('mean_for_mask_sub%02d_%s.nii', subNum, stim));
            if exist(meanCopyForMask,'file') ~= 2
                copyfile(meanFile, meanCopyForMask);
                fprintf('  Created drawing reference: %s\n', meanCopyForMask);
            else
                fprintf('  Drawing reference exists: %s\n', meanCopyForMask);
            end

            templateMask = fullfile(maskDir, ...
                sprintf('mask_template_sub%02d_%s.nii', subNum, stim));
            if exist(templateMask,'file') ~= 2
                V = spm_vol(meanFile);
                Z = zeros(V.dim);
                V.fname = templateMask;
                spm_write_vol(V, Z);
                fprintf('  Created mask template: %s\n', templateMask);
            else
                fprintf('  Mask template exists: %s\n', templateMask);
            end

            fprintf('  Draw and save the artifact mask as: %s\n', ...
                fullfile(maskDir, ...
                sprintf('mask_sub%02d_%s.nii', subNum, stim)));
        end

        fprintf('Stage 1 completed for %s\n', subName);
        continue;
    end

    yFile = ensure_segmentation_y(t1File);

    for s = 1:numel(STIMS)
        stim = STIMS{s};
        funcDir = fullfile(ROOT_DIR, stim, FUNC_FOLDER, subName);
        if exist(funcDir,'dir') ~= 7
            fprintf('[MISS] %s | directory not found: %s\n', subName, funcDir);
            continue;
        end

        maskDir = fullfile(ROOT_DIR, stim, MASK_FOLDER, subName);
        maskFile = fullfile(maskDir, ...
            sprintf('mask_sub%02d_%s.nii', subNum, stim));

        if exist(maskFile,'file') ~= 2 && strcmpi(stim,'stim0')
            legacyMask = fullfile(maskDir, sprintf('mask_sub%02d.nii', subNum));
            if exist(legacyMask,'file') == 2
                maskFile = legacyMask;
                fprintf('  Using legacy mask name: %s\n', maskFile);
            end
        end

        if exist(maskFile,'file') ~= 2
            fprintf('[SKIP] %s | %s | artifact mask not found: %s\n', ...
                subName, stim, maskFile);
            continue;
        end

        fprintf('\n--- %s | %s ---\n', subName, stim);
        [meanFile, rScans] = ensure_realign(funcDir, FORCE_REALIGN);

        workDir = fullfile(funcDir, '_dbs_work');
        if CLEAR_WORKDIR_EACH_RUN && exist(workDir,'dir') == 7
            rmdir(workDir,'s');
        end
        if exist(workDir,'dir') ~= 7
            mkdir(workDir);
        end

        if RESTORE_R_HEADERS
            backupFile = fullfile(funcDir, 'precoreg_realign_headers.mat');
            restore_or_backup_realign_headers(meanFile, rScans, backupFile);
        end

        maskCopy = fullfile(workDir, ...
            sprintf('mask_sub%02d_%s_src.nii', subNum, stim));
        copyfile(maskFile, maskCopy);

        P = char([meanFile ',1'], [maskCopy ',1']);
        flags = struct('mask',0,'mean',0,'interp',0,'which',1, ...
            'wrap',[0 0 0],'prefix','r');
        spm_reslice(P, flags);

        [~, maskBase, maskExt] = fileparts(maskCopy);
        rmaskCopy = fullfile(workDir, ['r' maskBase maskExt]);
        assert(exist(rmaskCopy,'file') == 2, ...
            'Resliced artifact mask was not created: %s', rmaskCopy);

        Vref = spm_vol(meanFile);
        Vrm = spm_vol(rmaskCopy);
        artifactVoxels = spm_read_vols(Vrm) > 0;

        maskEPI = fullfile(workDir, ...
            sprintf('maskEPI_sub%02d_%s.nii', subNum, stim));
        Vout = Vref;
        Vout.fname = maskEPI;
        spm_write_vol(Vout, double(artifactVoxels));

        weightFile = fullfile(workDir, ...
            sprintf('weight_sub%02d_%s.nii', subNum, stim));
        Vweight = Vref;
        Vweight.fname = weightFile;
        spm_write_vol(Vweight, double(~artifactVoxels));

        maskForCoreg = fullfile(workDir, ...
            sprintf('maskForCoreg_sub%02d_%s.nii', subNum, stim));
        copyfile(maskEPI, maskForCoreg);

        if USE_ENANTIOMORPHIC_FILL
            meanCopy = fullfile(workDir, ...
                sprintf('meanCopy_sub%02d_%s.nii', subNum, stim));
            copyfile(meanFile, meanCopy);

            maskForFill = fullfile(workDir, ...
                sprintf('maskForFill_sub%02d_%s.nii', subNum, stim));
            copyfile(maskEPI, maskForFill);

            intactImg = entiamorphicSub1(meanCopy, maskForFill);
            Vf = spm_vol(intactImg);
            filledData = spm_read_vols(Vf);

            sourceForCoreg = fullfile(workDir, ...
                sprintf('filledMeanForCoreg_sub%02d_%s.nii', subNum, stim));
            Vsource = Vref;
            Vsource.fname = sourceForCoreg;
            spm_write_vol(Vsource, filledData);
        else
            sourceForCoreg = meanFile;
        end

        fprintf('  Coregistering EPI to T1...\n');
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {t1File};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {sourceForCoreg};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = ...
            [rScans; {maskForCoreg}];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.0200 0.0200 0.0200 0.0010 0.0010 0.0010 ...
             0.0100 0.0100 0.0100 0.0010 0.0010 0.0010];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        if USE_COREG_WEIGHTING
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.weight = ...
                {weightFile};
        end
        spm_jobman('run', matlabbatch);

        fprintf('  Normalizing to MNI space at [%d %d %d] mm...\n', NORM_VOX);
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {yFile};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = rScans;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = NORM_VOX;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        spm_jobman('run', matlabbatch);

        fprintf('  Smoothing with [%d %d %d] mm FWHM...\n', SMOOTH_FWHM);
        wScans = cellstr(spm_select('ExtFPList', funcDir, '^wr.*\.nii$', Inf));
        assert(~isempty(wScans), 'No normalized wr*.nii files found in: %s', funcDir);

        matlabbatch = [];
        matlabbatch{1}.spm.spatial.smooth.data = wScans;
        matlabbatch{1}.spm.spatial.smooth.fwhm = SMOOTH_FWHM;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run', matlabbatch);

        fprintf('  Warping the artifact mask and creating a GLM mask...\n');
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {yFile};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ...
            {maskForCoreg};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = NORM_VOX;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        spm_jobman('run', matlabbatch);

        wMaskTmp = spm_file(maskForCoreg, 'prefix', 'w');
        if exist(wMaskTmp,'file') == 2
            if exist(maskDir,'dir') ~= 7
                mkdir(maskDir);
            end

            outWMask = fullfile(maskDir, ...
                sprintf('wmask_sub%02d_%s.nii', subNum, stim));
            copyfile(wMaskTmp, outWMask);

            Vwm = spm_vol(outWMask);
            normalizedArtifactVoxels = spm_read_vols(Vwm) > 0;

            Vglm = Vwm;
            Vglm.fname = fullfile(maskDir, ...
                sprintf('glm_mask_sub%02d_%s.nii', subNum, stim));
            spm_write_vol(Vglm, double(~normalizedArtifactVoxels));

            fprintf('    Saved normalized artifact mask: %s\n', outWMask);
            fprintf('    Saved GLM explicit mask: %s\n', Vglm.fname);
        else
            warning('Normalized artifact mask was not created: %s', wMaskTmp);
        end

        fprintf('Completed %s | %s\n', subName, stim);
    end
end

fprintf('\nAll preprocessing completed.\n');
fprintf('GLM images: swr*.nii in each functional directory.\n');
fprintf('Explicit masks: glm_mask_subXX_stimX.nii in each mask directory.\n');

function directoryPath = request_directory(promptText)
directoryPath = strtrim(input(promptText, 's'));
assert(~isempty(directoryPath), 'A directory path is required.');
assert(exist(directoryPath,'dir') == 7, ...
    'Directory not found: %s', directoryPath);
end

function subNum = parse_sub_number(subName)
token = regexp(subName,'^sub(\d+)$','tokens','once');
if isempty(token)
    subNum = [];
    return;
end
subNum = str2double(token{1});
if isnan(subNum)
    subNum = [];
end
end

function t1File = find_t1_file_robust(t1Dir)
assert(exist(t1Dir,'dir') == 7, 'T1 directory not found: %s', t1Dir);

nii = dir(fullfile(t1Dir,'*.nii'));
assert(~isempty(nii), 'No .nii files found in T1 directory: %s', t1Dir);

keep = false(size(nii));
score = zeros(size(nii));

for k = 1:numel(nii)
    name = nii(k).name;
    lowerName = lower(name);

    if contains(lowerName,'crop')
        continue;
    end
    if startsWith(lowerName,'y_') || startsWith(lowerName,'iy_')
        continue;
    end
    if ~isempty(regexp(lowerName,'^c[1-6]','once'))
        continue;
    end
    if startsWith(lowerName,'m') && contains(lowerName,'t1')
        continue;
    end
    if startsWith(lowerName,'w') || startsWith(lowerName,'sw') || ...
            startsWith(lowerName,'r')
        continue;
    end

    try
        V = spm_vol(fullfile(t1Dir,name));
        if numel(V) ~= 1
            continue;
        end
    catch
        continue;
    end

    keep(k) = true;
    if contains(lowerName,'t1')
        score(k) = score(k) + 2;
    end
    if contains(lowerName,'mprage') || contains(lowerName,'spgr')
        score(k) = score(k) + 2;
    end
    if contains(lowerName,'anat')
        score(k) = score(k) + 1;
    end
end

candidates = nii(keep);
scores = score(keep);
assert(~isempty(candidates), 'No eligible original T1 image found in: %s', t1Dir);

bytes = reshape([candidates.bytes], [], 1);
rankTable = [-scores(:), -bytes];
[~, order] = sortrows(rankTable, [1 2]);
candidates = candidates(order);

t1File = fullfile(t1Dir, candidates(1).name);
fprintf('  Using T1: %s\n', t1File);
end

function restore_or_backup_realign_headers(meanFileWithIdx, rScansWithIdx, backupFile)
meanFile = regexprep(meanFileWithIdx, ',\d+$', '');
rFiles = cellfun(@(x) regexprep(x, ',\d+$', ''), rScansWithIdx, ...
    'UniformOutput', false);
rFiles = unique(rFiles);

if exist(backupFile,'file') == 2
    S = load(backupFile);
    mismatch = ~isfield(S,'meanFile') || ~isfield(S,'meanMat') || ...
        ~isfield(S,'rFiles') || ~isfield(S,'rMats') || ...
        ~strcmpi(S.meanFile, meanFile) || numel(S.rFiles) ~= numel(rFiles);

    if mismatch
        fprintf('  Header backup mismatch; creating a new backup.\n');
        S = collect_header_data(meanFile, rFiles);
        save(backupFile, '-struct', 'S');
        return;
    end

    spm_get_space(meanFile, S.meanMat);
    for i = 1:numel(S.rFiles)
        if exist(S.rFiles{i},'file') == 2
            spm_get_space(S.rFiles{i}, S.rMats{i});
        end
    end
    fprintf('  Restored realignment headers from: %s\n', backupFile);
else
    S = collect_header_data(meanFile, rFiles);
    save(backupFile, '-struct', 'S');
    fprintf('  Saved realignment-header backup: %s\n', backupFile);
end
end

function S = collect_header_data(meanFile, rFiles)
S.meanFile = meanFile;
S.meanMat = spm_get_space(meanFile);
S.rFiles = rFiles;
S.rMats = cell(size(rFiles));
for i = 1:numel(rFiles)
    S.rMats{i} = spm_get_space(rFiles{i});
end
end

function yFile = ensure_segmentation_y(t1File)
[t1Dir, t1Base, t1Ext] = fileparts(t1File);
yFile = fullfile(t1Dir, ['y_' t1Base t1Ext]);

if exist(yFile,'file') == 2
    fprintf('  Found deformation field: %s\n', yFile);
    return;
end

fprintf('  Segmenting T1 and estimating deformation fields...\n');
tpm = fullfile(spm('Dir'),'tpm','TPM.nii');
assert(exist(tpm,'file') == 2, 'SPM TPM file not found: %s', tpm);

matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols = {t1File};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];

ngaus = [1 1 2 3 4 2];
for t = 1:6
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = ...
        {sprintf('%s,%d', tpm, t)};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = ngaus(t);
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end

matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('run', matlabbatch);

assert(exist(yFile,'file') == 2, ...
    'T1 segmentation did not create the deformation field: %s', yFile);
fprintf('  Created deformation field: %s\n', yFile);
end

function [meanFile, rScans] = ensure_realign(funcDir, forceRealign)
inputScans = spm_select('ExtFPList', funcDir, '^[Aa].*\.nii$', Inf);
if isempty(inputScans)
    error('No slice-time-corrected A*.nii files found in: %s', funcDir);
end
inputScans = cellstr(inputScans);

meanList = dir(fullfile(funcDir,'mean*.nii'));
rList = dir(fullfile(funcDir,'r*.nii'));

if ~forceRealign && ~isempty(meanList) && ~isempty(rList)
    [~, idx] = max([meanList.datenum]);
    meanFile = fullfile(funcDir, meanList(idx).name);
    rScans = cellstr(spm_select('ExtFPList', funcDir, '^r.*\.nii$', Inf));
    if isempty(rScans)
        error('r*.nii files exist but could not be selected in: %s', funcDir);
    end
    fprintf('  Reusing existing realignment outputs. Mean image: %s\n', meanFile);
    return;
end

fprintf('  Running realignment (estimate and reslice)...\n');
matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = inputScans;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run', matlabbatch);

meanList = dir(fullfile(funcDir,'mean*.nii'));
assert(~isempty(meanList), 'Realignment did not create mean*.nii in: %s', funcDir);
[~, idx] = max([meanList.datenum]);
meanFile = fullfile(funcDir, meanList(idx).name);

rScans = cellstr(spm_select('ExtFPList', funcDir, '^r.*\.nii$', Inf));
assert(~isempty(rScans), 'Realignment did not create r*.nii in: %s', funcDir);
fprintf('  Realignment completed. Mean image: %s | Volumes: %d\n', ...
    meanFile, numel(rScans));
end
