% % MFR_main.m 
% %
% % Edited by Minseung Kim, 2025-08-13 % %
% % Edited by Minseung Kim, 2025-11-26 % %

function MFR_main()

repoRoot = fileparts(mfilename('fullpath'));

S = struct();
S.repoRoot = repoRoot;

S.fig = figure( ...
    'Name', 'MuscleFreqResp – Frequency Response Simulator', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'off', ...
    'Position', [100 100 900 630]);   % [left bottom width height]
   
uicontrol(S.fig, ...
    'Style', 'text', ...
    'String', 'MuscleFreqResp – Simulation GUI', ...
    'FontSize', 19, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'Position', [10 540 860 40]);

S.p1 = uipanel(S.fig, ...
    'Title', 'Phase 1 - Choose Model / Mode / Variant', ...
    'TitlePosition','centertop', ...
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Units','pixels', ...
    'Position', [40 110 380 410]);   % left, bottom, width, height

S.p2 = uipanel(S.fig, ...
    'Title', 'Phase 2 - Simulation Parameters', ...
    'TitlePosition','centertop', ...
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Units','pixels', ...
    'Position', [440 110 420 410]);

% === Phase 1 controls (left) ===========================================
y0   = 250;
dy   = 50;
xL   = 25;
wLbl = 90;
wCtl = 220;
fs   = 10;

% === Config mode: Preset / Custom ======================================
uicontrol(S.p1, 'Style','text', ...
    'String','Config mode', ...
    'FontSize',fs, ...
    'HorizontalAlignment','left', ...
    'Position',[xL, y0+dy, wLbl, 22]);

S.popupConfigMode = uicontrol(S.p1, 'Style','popupmenu', ...
    'String',{'Paper preset (reproduction)','Custom (not yet implemented)'}, ...
    'FontSize',fs, ...
    'Position',[xL+wLbl, y0+dy, wCtl, 24]);

uicontrol(S.p1, 'Style','text', ...
    'String','Model', ...
    'FontSize',fs, ...
    'HorizontalAlignment','left', ...
    'Position',[xL y0 wLbl 22]);

S.popupModel = uicontrol(S.p1, 'Style','popupmenu', ...
    'String',{'TMM (Thelen)','MMM (Millard)'}, ...
    'FontSize',fs, ...
    'Position',[xL+wLbl y0 wCtl 26], ...
    'Callback',@(h,~)onModelOrModeChange(h));

uicontrol(S.p1, 'Style','text', ...
    'String','Mode', ...
    'FontSize',fs, ...
    'HorizontalAlignment','left', ...
    'Position',[xL y0-dy wLbl 22]);

S.popupMode = uicontrol(S.p1, 'Style','popupmenu', ...
    'String',{'active','passive'}, ...
    'FontSize',fs, ...
    'Position',[xL+wLbl y0-dy wCtl 26], ...
    'Callback',@(h,~)onModelOrModeChange(h));

uicontrol(S.p1, 'Style','text', ...
    'String','MMM variant', ...
    'FontSize',fs, ...
    'HorizontalAlignment','left', ...
    'Position',[xL y0-2*dy wLbl 22]);

S.popupMMM = uicontrol(S.p1, 'Style','popupmenu', ...
    'String',{'(not used for TMM)'}, ...
    'Enable','off', ...
    'FontSize',fs, ...
    'Position',[xL+wLbl y0-2*dy wCtl 26]);

S.textInfo = uicontrol(S.p1, 'Style','text', ...
    'String', 'Choose configuration and press "Choose" to unlock Phase 2.', ...
    'FontSize',fs, ...
    'HorizontalAlignment','left', ...
    'Position',[xL 95 320 35]);

S.btnChoose = uicontrol(S.p1, 'Style', 'pushbutton', ...
    'String', 'Choose', ...
    'FontSize', fs+1, ...
    'FontWeight', 'bold', ...
    'Position', [xL+20 20 130 32], ...
    'Callback', @(h,~)onChoosePressed(h));

S.btnClose  = uicontrol(S.p1, 'Style', 'pushbutton', ...
    'String', 'Close', ...
    'FontSize', fs+1, ...
    'FontWeight', 'bold', ...
    'Position', [xL+180 20 130 32], ...
    'Callback', @(h,~)onClosePressed(h));

% === Phase 2 controls (right) ==========================================
xR    = 30;
yR0   = 350;
dy2   = 32;
wLbl2 = 160;
wEdit = 200;
fs2   = 10;

uicontrol(S.p2,'Style','text', ...
    'String','L_mn0_values (row)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0 wLbl2 22]);

S.editLmn0 = uicontrol(S.p2,'Style','edit', ...
    'String','1.0', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0 wEdit 24]);

uicontrol(S.p2,'Style','text', ...
    'String','V_m0_values (row)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0-dy2 wLbl2 22]);

S.editVm0 = uicontrol(S.p2,'Style','edit', ...
    'String','0.0', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0-dy2 wEdit 24]);

uicontrol(S.p2,'Style','text', ...
    'String','u0 (activation)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0-2*dy2 wLbl2 22]);

S.editU0 = uicontrol(S.p2,'Style','edit', ...
    'String','0.5', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0-2*dy2 wEdit 24]);

uicontrol(S.p2,'Style','text', ...
    'String','const_b (excitation amp)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0-3*dy2 wLbl2 22]);

S.editB = uicontrol(S.p2,'Style','edit', ...
    'String','0.01', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0-3*dy2 wEdit 24]);

uicontrol(S.p2,'Style','text', ...
    'String','mass (kg or N·s²/m)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0-4*dy2 wLbl2 22]);

S.editMass = uicontrol(S.p2,'Style','edit', ...
    'String','3e6', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0-4*dy2 wEdit 24]);

uicontrol(S.p2,'Style','text', ...
    'String','damping (N·s/m)', ...
    'FontSize',fs2, ...
    'HorizontalAlignment','left', ...
    'Position',[xR yR0-5*dy2 wLbl2 22]);

S.editDamp = uicontrol(S.p2,'Style','edit', ...
    'String','0.0', ...
    'FontSize',fs2, ...
    'Enable','off', ...
    'Position',[xR+wLbl2 yR0-5*dy2 wEdit 24]);

uicontrol(S.p2, 'Style', 'text', ...
    'String', 'totalTime (s)', ...
    'FontSize', fs2, ...
    'HorizontalAlignment', 'left', ...
    'Position', [xR yR0-6*dy2 wLbl2 22]);

S.editTotalT = uicontrol(S.p2, 'Style', 'edit', ...
    'String', '120', ...
    'FontSize', fs2, ...
    'Enable', 'off', ...
    'Position', [xR+wLbl2 yR0-6*dy2 wEdit 24]);

uicontrol(S.p2, 'Style', 'text', ...
    'String', 'dt (s)', ...
    'FontSize', fs2, ...
    'HorizontalAlignment', 'left', ...
    'Position', [xR yR0-7*dy2 wLbl2 22]);

S.editDt = uicontrol(S.p2, 'Style', 'edit', ...
    'String', '0.001', ...
    'FontSize', fs2, ...
    'Enable', 'off', ...
    'Position', [xR+wLbl2 yR0-7*dy2 wEdit 24]);

uicontrol(S.p2, 'Style', 'text', ...
    'String', 'freq range [Hz] (min, max, N)', ...
    'FontSize', fs2, ...
    'HorizontalAlignment', 'left', ...
    'Position', [xR yR0-8*dy2 wLbl2+40 22]);

S.editFreq = uicontrol(S.p2, 'Style', 'edit', ...
    'String', '0.1, 100, 400', ...
    'FontSize', fs2, ...
    'Enable', 'off', ...
    'Position', [xR+wLbl2+40 yR0-8*dy2 wEdit-40 24]);

S.textPhase2 = uicontrol(S.p2, 'Style', 'text', ...
    'String', 'In passive mode, u0 and const_b are fixed to 0 and disabled.', ...
    'FontSize', fs2, ...
    'HorizontalAlignment', 'center', ...
    'Position', [xR-10 20 380 60]);

S.btnRun = uicontrol(S.p2, 'Style', 'pushbutton', ...
    'String', 'Run Simulation', ...
    'FontSize', fs2 + 1, ...
    'FontWeight', 'bold', ...
    'Enable', 'off', ...
    'Position', [xR+20 20 170 32], ...
    'Callback', @(h,~)onRunPressed(h));

S.btnStop = uicontrol(S.p2, 'Style', 'pushbutton', ...
    'String', 'Stop', ...
    'FontSize', fs2 + 1, ...
    'FontWeight', 'bold', ...
    'Enable', 'off', ...
    'Position', [xR+230 20 110 32], ...
    'Callback', @(h,~)onStopPressed(h));

S.textStatus = uicontrol(S.fig, 'Style', 'text', ...
    'Style', 'text', ...
    'String', 'Status: Ready', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'Position', [20 65 860 22]);

S.btnViewer = uicontrol(S.fig, ...
    'Style', 'pushbutton', ...
    'String', 'Open Bode viewer', ...
    'FontSize', 10, ...
    'Position', [25 15 150 25], ... 
    'Callback', @(h,~)onOpenViewerPressed(h));

S.textFooter = uicontrol(S.fig, ...
    'Style', 'text', ...
    'String', 'Last updated: 2025-12-04', ... 
    'FontSize', 8, ...
    'HorizontalAlignment', 'right', ...
    'Position', [20 40 860 12]);

S.model     = 'TMM';    
S.mode      = 'active';
S.entryPath = '';
S.modelDir  = '';
S.setupMode = 1;
S.mmmListActive  = {};
S.mmmListPassive = {};

S.isRunning = false;

guidata(S.fig, S);
refreshMMMVariantPopup(S.fig);

end

function onModelOrModeChange(h)
    fig = ancestor(h, 'figure');
    S   = guidata(fig);

    modelIdx = get(S.popupModel, 'Value');
    modeIdx  = get(S.popupMode,  'Value');

    S.model = tern(modelIdx==2, 'MMM', 'TMM');
    S.mode  = tern(modeIdx==2,  'passive', 'active');

    guidata(fig, S);

    if strcmp(S.model,'TMM')
        set(S.popupMMM, 'Enable','off', ...
            'String', {'(not used for TMM)'}, ...
            'Value', 1);

        S = guidata(fig);
        disablePhase2(S);
    else
        set(S.popupMMM, 'Enable','on');
        refreshMMMVariantPopup(fig);   

        S = guidata(fig);
        disablePhase2(S);         
    end
end

function refreshMMMVariantPopup(fig)
    S = guidata(fig);
    
    modelDir = fullfile(S.repoRoot, 'MMM test', 'MMM', 'src');
    names = {};
    
    if exist(modelDir,'dir') == 7
        listing = dir(fullfile(modelDir,'MMM*.m'));
        names   = {listing.name};
    end
    
    if isempty(names)
        set(S.popupMMM,'String',{'(no MMM*.m found)'},'Value',1);
        guidata(fig,S);
        return;
    end
    
    isPassive = contains(names,'passive','IgnoreCase',true);
    
    if strcmp(S.mode,'passive')
        candidates = names(isPassive);
    else
        candidates = names(~isPassive);
    end
    
    if isempty(candidates)
        if strcmp(S.mode,'passive')
            msg = '(no passive MMM*.m)';
        else
            msg = '(no active MMM*.m)';
        end
        set(S.popupMMM,'String',{msg},'Value',1);
    else
        set(S.popupMMM,'String',candidates,'Value',1);
    end
    
    guidata(fig,S);
end

function onChoosePressed(h)
    fig = ancestor(h,'figure');
    S = guidata(fig);
    
    cfgIdx      = get(S.popupConfigMode,'Value');   % 1: preset, 2: custom
    S.setupMode = cfgIdx;
    
    if cfgIdx == 2
        warndlg({'Custom mode is not implemented yet.', ...
                 'This version supports only "Paper preset (reproduction)" mode.'}, ...
                'Not yet implemented','modal');
        if isfield(S,'textStatus') && ishandle(S.textStatus)
            set(S.textStatus,'String', ...
                'Status: Custom mode is not implemented in this version.');
        end
        guidata(fig,S);
        return; 
    end
    
    % Model / Mode
    modelIdx = get(S.popupModel,'Value');
    modeIdx  = get(S.popupMode, 'Value');
    
    S.model = tern(modelIdx==2, 'MMM', 'TMM');
    S.mode  = tern(modeIdx==2,  'passive', 'active');
    
    repoRoot = S.repoRoot;
    
    switch S.model
        case 'TMM'
            modelDir = fullfile(repoRoot, 'TMM test', 'TMM');
            if strcmp(S.mode,'passive')
                entry = 'TMM_passive.m';
            else
                entry = 'TMM_active.m';
            end
    
        case 'MMM'
            modelDir = fullfile(repoRoot, 'MMM test', 'MMM', 'src');
            listStr  = get(S.popupMMM,'String');
            if ischar(listStr)
                candidates = {listStr};
            else
                candidates = listStr;
            end
            idx = get(S.popupMMM,'Value');
            if isempty(candidates) || idx<1 || idx>numel(candidates)
                errordlg('No valid MMM variant selected.','Error','modal');
                return;
            end
            entry = candidates{idx};
    
        otherwise
            errordlg('Unknown model selected.','Error','modal');
            return;
    end
    
    if exist(modelDir,'dir') ~= 7
        errordlg(sprintf('Model directory not found:\n%s', modelDir), ...
                 'Error','modal');
        return;
    end
    
    entryPath = fullfile(modelDir, entry);
    if exist(entryPath,'file') ~= 2
        errordlg(sprintf('Entry file not found:\n%s', entryPath), ...
                 'Error','modal');
        return;
    end
    
    S.modelDir  = modelDir;
    S.entryPath = entryPath;
    
    configurePhase2Defaults(S);
    enablePhase2(S);
    
    guidata(fig,S);
end

function configurePhase2Defaults(S)

    LmnStr = '1.0';
    VmStr  = '0.0';
    
    if strcmp(S.mode,'passive')
        u0Str = '0.0';
        bStr  = '0.0';
    else
        u0Str = '0.5';
        bStr  = '0.01';
    end
    
    if strcmp(S.mode,'passive')
        massStr = '0';       
        dampStr = '0';      
    else
        massStr = '3e6';      
        dampStr = '0.0';
    end
    
    totalTStr   = '120';
    dtStr       = '0.001';
    freqStr     = '0.1, 100, 400';
    
    set(S.editLmn0,     'String', LmnStr);
    set(S.editVm0,      'String', VmStr);
    set(S.editU0,       'String', u0Str);
    set(S.editB,        'String', bStr);
    set(S.editMass,     'String', massStr);
    set(S.editDamp,     'String', dampStr);
    set(S.editTotalT,   'String', totalTStr);
    set(S.editDt,       'String', dtStr);
    set(S.editFreq,     'String', freqStr);
end

function enablePhase2(S)
    set(S.editLmn0,     'Enable', 'on');
    set(S.editVm0,      'Enable', 'on');
    set(S.editMass,     'Enable', 'on');
    set(S.editDamp,     'Enable', 'on');
    set(S.editTotalT,   'Enable', 'on');
    set(S.editDt,       'Enable', 'on');
    set(S.editFreq,     'Enable', 'on');
    set(S.btnRun,       'Enable', 'on');
    set(S.btnStop,      'Enable', 'off');

    if strcmp(S.mode,   'passive')
        set(S.editU0,   'Enable', 'off');
        set(S.editB,    'Enable', 'off');
    else
        set(S.editU0,   'Enable', 'on');
        set(S.editB,    'Enable', 'on');
    end
end

function disablePhase2(S)
    set(S.editLmn0,     'Enable', 'off');
    set(S.editVm0,      'Enable', 'off');
    set(S.editU0,       'Enable', 'off');
    set(S.editB,        'Enable', 'off');
    set(S.editMass,     'Enable', 'off');
    set(S.editDamp,     'Enable', 'off');
    set(S.editTotalT,   'Enable', 'off');
    set(S.editDt,       'Enable', 'off');
    set(S.editFreq,     'Enable', 'off');
    set(S.btnRun,       'Enable', 'off');
    set(S.btnStop,      'Enable', 'off');
end

function onRunPressed(h)

    fig = ancestor(h,'figure');
    S = guidata(fig);
    
    if isfield(S,'isRunning') && S.isRunning
        set(S.textStatus,'String','Status: Simulation is already running.');
        return;
    end

    if isempty(S.entryPath) || isempty(S.modelDir)
        errordlg('No configuration selected in Phase 1. Press "Choose" first.', ...
                 'Error','modal');
        return;
    end

    LmnStr      = get(S.editLmn0,'String');
    VmStr       = get(S.editVm0, 'String');
    u0Str       = get(S.editU0,  'String');
    bStr        = get(S.editB,   'String');
    massStr     = get(S.editMass,'String');
    dampStr     = get(S.editDamp,'String');
    totStr      = get(S.editTotalT,'String');
    dtStr       = get(S.editDt,    'String');
    freqStr     = get(S.editFreq,  'String');
    
    LmnVals     = str2num(LmnStr);
    VmVals      = str2num(VmStr);  
    u0Val       = str2double(u0Str);
    bVal        = str2double(bStr);
    massVal     = str2double(massStr);
    dampVal     = str2double(dampStr);
    totVal      = str2double(totStr);
    dtVal       = str2double(dtStr);
    
    if isempty(freqStr)
        error('Frequency range is empty. Please specify as "min, max, steps".');
    end

    freqStrClean = strrep(freqStr, ',', ' ');
    freqTokens = sscanf(freqStrClean, '%f');
    
    if numel(freqTokens) ~= 3 || any(isnan(freqTokens))
        error('Invalid frequency range. Use "min, max, steps" (e.g., "0.1, 100, 400").');
    end

    freqlbVal = freqTokens(1);
    freqhbVal = freqTokens(2);
    stepsVal  = round(freqTokens(3));

    if isempty(LmnVals) || isempty(VmVals) || isnan(u0Val) || isnan(bVal) || ...
       isnan(massVal) || isnan(dampVal) || isnan(totVal) || isnan(dtVal)
        errordlg('One or more parameters are invalid. Check inputs.','Error','modal');
        return;
    end

    if strcmp(S.mode,'passive')
        u0Val = 0.0;
        bVal  = 0.0;
    end

    assignin('base', 'L_mn0_values', LmnVals);
    assignin('base', 'V_m0_values',  VmVals);
    assignin('base', 'u0',           u0Val);
    assignin('base', 'const_b',      bVal);
    assignin('base', 'mass',         massVal);
    assignin('base', 'damp',         dampVal);
    assignin('base', 'ext_damping',  dampVal);
    assignin('base', 'totalTime',    totVal);
    assignin('base', 'dt',           dtVal);
    assignin('base', 'freqlb',       freqlbVal);
    assignin('base', 'freqhb',       freqhbVal);
    assignin('base', 'steps',        stepsVal);
    assignin('base', 'MFR_STOP',     false);

    L_mn0_values    = LmnVals;
    V_m0_values     = VmVals;
    u0              = u0Val;
    const_b         = bVal;
    mass            = massVal;
    damp            = dampVal;
    ext_damping     = dampVal;
    totalTime       = totVal;
    dt              = dtVal;
    freqlb          = freqlbVal;
    freqhb          = freqhbVal;
    steps           = stepsVal;

    S.isRunning     = true;
    set(S.btnRun, 'Enable','off');
    set(S.btnStop,'Enable','on');

    [~,entryName,entryExt] = fileparts(S.entryPath);
    set(S.textStatus,'String',sprintf('Status: Running %s', [entryName entryExt]));
    guidata(fig,S);

    drawnow;
    
    fprintf('\n>> Running %s\n   Model=%s | Mode=%s\n', ...
            S.entryPath, S.model, S.mode);
    
    try
        run(S.entryPath);
        msg = sprintf('Status: Finished %s', [entryName entryExt]);
    catch ME
        msg = sprintf('Status: Error in %s (%s)', [entryName entryExt], ME.message);
    end
    
    S = guidata(fig);
    S.isRunning = false;
    set(S.btnRun, 'Enable','on');
    set(S.btnStop,'Enable','off');
    set(S.textStatus,'String',msg);
    
    guidata(fig,S);
end

function onStopPressed(h)
    fig = ancestor(h,'figure');
    S   = guidata(fig);
    
    disp('[DEBUG] Stop button pressed');

    if ~isfield(S,'isRunning') || ~S.isRunning
        set(S.textStatus,'String','Status: No simulation is running.');
        return;
    end

    set(S.textStatus,'String','Status: Stop requested (will stop after current step)...');
    drawnow;

    assignin('base', 'MFR_STOP', true);
end

function onClosePressed(h)
    fig = ancestor(h,'figure');
    delete(fig);
end

function onOpenViewerPressed(h)
    fig = ancestor(h,'figure');
    S   = guidata(fig);

    try
        fullpaths = {};

        keepSelecting = true;
        lastPath = S.repoRoot;

        while keepSelecting
            [files, pathName] = uigetfile( ...
                '*.mat', ...
                'Select Bode data files (afterFFT .mat)', ...
                lastPath, ...
                'MultiSelect','on');

            if isequal(files,0)
                break;
            end

            if ischar(files)
                files = {files};
            end

            for i = 1:numel(files)
                fullpaths{end+1,1} = fullfile(pathName, files{i});
            end

            lastPath = pathName;

            answer = questdlg( ...
                'Select more files from another folder?', ...
                'Add files', ...
                'Yes','No','No');

            if isempty(answer) || strcmp(answer,'No')
                keepSelecting = false;
            end
        end

        if isempty(fullpaths)
            disp('[INFO] Bode viewer: no files selected.');
            return;
        end

        assignin('base','GUI_Bode_files', fullpaths);

        viewerPath = fullfile(S.repoRoot, 'MFR_Bode_viewer.m');
        if exist(viewerPath,'file') ~= 2
            error('MFR_Bode_viewer.m not found at %s', viewerPath);
        end

        run(viewerPath);

    catch ME
        warning(sprintf('Failed to open Bode viewer: %s', ME.message));
    end
end

function y = tern(cond,a,b)
if cond, y = a; else, y = b; end
end
