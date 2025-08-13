% % main.m 
% %
% % Edited by Minseung Kim, 2025-08-13 % %
%
% % This code manage the simulation environment of TMM (Thelen Muscle Model)
% % We adjusted the dynamic equation from the literature of D. G. Thelen (2003)
%
% % Entry point for MuscleFreqResp repository.
% % Allows selection of model (TMM or MMM) and mode (active or passive).
% % For MMM, prompts for a specific variant if not provided.

% [Getting started]

function main(varargin)

p = inputParser;

addParameter(p, 'Model', '',          @(s)ischar(s)||isstring(s));
addParameter(p, 'Mode',  '',          @(s)ischar(s)||isstring(s));
addParameter(p, 'MMMVariant', '',     @(s)ischar(s)||isstring(s));
addParameter(p, 'ForcePreview', true, @(x)islogical(x)||isnumeric(x));

parse(p, varargin{:});
opt = p.Results;

opt.Model           = upper(char(opt.Model));
opt.Mode            = lower(char(opt.Mode));
opt.MMMVariant      = char(opt.MMMVariant);
opt.ForcePreview    = logical(opt.ForcePreview);

repoRoot        = fileparts(mfilename('fullpath'));

if isempty(opt.Model)
    fprintf('Select model:\n  1) TMM (Thelen)\n  2) MMM (Millard)\n');
    m           = safePick(2, 1, 'Enter choice [1-2]: ');
    opt.Model   = tern(m==2, 'MMM', 'TMM');
end

if isempty(opt.Mode)
    fprintf('Select simulation mode:\n  1) active\n  2) passive\n');
    c           = safePick(2, 1, 'Enter choice [1-2]: ');
    opt.Mode    = tern(c==2, 'passive', 'active');
end

switch opt.Model
    case 'TMM'
        modelDir    = fullfile(repoRoot, 'TMM test', 'TMM');
        if strcmp(opt.Mode, 'passive')
            entry   = 'TMM_passive_ver1.m';
        else
            entry   = 'TMM_ver1.m';
        end

    case 'MMM'
        modelDir    = fullfile(repoRoot, 'MMM test', 'MMM', 'src');
        listing     = dir(fullfile(modelDir, 'MMM*.m'));
        names       = {listing.name};

        if isempty(names)
            error('No MMM files found in %s', modelDir);
        end

        isPassive           = contains(names, 'passive', 'IgnoreCase', true);
        candidatesActive    = names(~isPassive);
        candidatesPassive   = names( isPassive);
        candidates          = iff(strcmp(opt.Mode,'passive'), candidatesPassive, candidatesActive);

        if isempty(candidates)
            error('No MMM %s variants found in %s', opt.Mode, modelDir);
        end

        fprintf('\n=== MMM %s variant list ===\n', upper(opt.Mode));
        for i = 1:numel(candidates)
            fprintf('  %d) %s\n', i, candidates{i});
        end

        if ~isempty(opt.MMMVariant)
            fprintf('\nPreset MMMVariant: %s\n', opt.MMMVariant);
        end

        if opt.ForcePreview
            if isempty(opt.MMMVariant)
                k       = safePick(numel(candidates), 1, 'Enter choice number: ');
                entry   = candidates{k};
            else
                inp     = input('Press ENTER to use preset, or enter number to override: ', 's');
                if isempty(inp)
                    entry = opt.MMMVariant;
                else
                    k       = str2double(inp);
                    assert(~isnan(k) && k>=1 && k<=numel(candidates), 'Invalid selection.');
                    entry   = candidates{k};
                end
            end
         
        else
            entry = ifempty(opt.MMMVariant, candidates{1});
        end

    otherwise
        error('Unknown Model: %s (use ''TMM'' or ''MMM'')', opt.Model);
end

if exist(modelDir, 'dir') ~= 7
    error('Model directory not found: %s', modelDir);
end
addpath(genpath(modelDir));

assert(exist(modelDir,'dir') == 7, 'Model dir not found: %s', modelDir);
addpath(genpath(modelDir));

entryPath = fullfile(modelDir, entry);
assert(exist(entryPath,'file') == 2, 'Entry file not found: %s', entryPath);

fprintf('\n>> Running %s | Model=%s | Mode=%s\n', entry, opt.Model, opt.Mode);
run(entryPath);

end

% % Helpers

function y = tern(cond, a, b),  if cond, y=a; else, y=b; end, end

function y = iff(cond, a, b),   if cond, y=a; else, y=b; end, end

function y = ifempty(x, d),     if isempty(x), y=d; else, y=x; end, end

function k = safePick(n, defaultIdx, prompt)
    if nargin<2 || isempty(defaultIdx), defaultIdx=1; end
    if nargin<3, prompt='Select: '; end
    s = input(prompt, 's');
    if isempty(s), k = defaultIdx; return; end
    k = str2double(s);
    assert(~isnan(k) && k>=1 && k<=n, 'Invalid selection.');
end