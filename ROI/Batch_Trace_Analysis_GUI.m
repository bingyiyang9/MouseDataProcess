classdef Batch_Trace_Analysis_GUI < matlab.apps.AppBase
    % Batch Trace Analysis GUI
    % 基于 ReadImageJROI 工具箱进行批量 Tiff 灰度提取
    
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        GridLayout           matlab.ui.container.GridLayout
        LeftPanel            matlab.ui.container.Panel
        RightPanel           matlab.ui.container.Panel
        
        % UI 组件
        BtnSelectTiffs       matlab.ui.control.Button
        BtnSelectROI         matlab.ui.control.Button
        ListBoxFiles         matlab.ui.control.ListBox
        LblROIStatus         matlab.ui.control.Label
        SpinnerFPS           matlab.ui.control.Spinner
        LblFPS               matlab.ui.control.Label
        BtnRun               matlab.ui.control.Button
        TextAreaLog          matlab.ui.control.TextArea
        
        UIAxes               matlab.ui.control.UIAxes
    end
    
    properties (Access = private)
        TiffFiles = {};      % 存储选中的 TIFF 文件全路径
        TiffNames = {};      % 文件名用于显示
        ROIPath = '';        % ROI 文件路径
        ROIData = {};        % 解析后的 ROI 数据
        ROI_Masks = {};      % 转换后的逻辑掩膜
        ImgWidth = 0;
        ImgHeight = 0;
    end
    
    methods (Access = private)
        
        % === 1. 选择 TIFF 文件 (支持多选) ===
        function onSelectTiffs(app, ~, ~)
            [files, path] = uigetfile({'*.tif;*.tiff', 'Tiff Stack'}, ...
                                      'Select TIFF Files (Multi-Select)', 'MultiSelect', 'on');
            if isequal(files, 0), return; end
            
            % 处理单选和多选的情况
            if ischar(files)
                files = {files};
            end
            
            % --- 关键：自然排序 (Natural Sort) ---
            % 防止 Video_10 排在 Video_2 前面
            [~, idx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), files));
            files = files(idx);
            
            app.TiffFiles = fullfile(path, files);
            app.TiffNames = files;
            app.ListBoxFiles.Items = files;
            
            app.log(sprintf('已加载 %d 个 TIFF 文件。', length(files)));
            
            % 预读取第一个文件的尺寸，用于生成 Mask
            info = imfinfo(app.TiffFiles{1});
            app.ImgWidth = info(1).Width;
            app.ImgHeight = info(1).Height;
        end
        
        % === 2. 选择 ROI 文件 ===
        function onSelectROI(app, ~, ~)
            % 检查工具箱
            if ~exist('ReadImageJROI', 'file')
                uialert(app.UIFigure, '未找到 ReadImageJROI 工具箱！请先安装。', 'Error');
                return;
            end
            
            if app.ImgWidth == 0
                uialert(app.UIFigure, '请先选择 TIFF 文件，以便确定 ROI 尺寸。', 'Warning');
                return;
            end
            
            [f, p] = uigetfile({'*.zip;*.roi', 'ImageJ ROI (*.zip)'}, 'Select ROI File');
            if isequal(f, 0), return; end
            
            app.ROIPath = fullfile(p, f);
            app.log(['正在解析 ROI: ' f ' ...']);
            
            try
                % 调用工具箱
                rawROIs = ReadImageJROI(app.ROIPath);
                if isstruct(rawROIs), rawROIs = {rawROIs}; end % 统一格式
                
                app.ROIData = rawROIs;
                nROIs = length(rawROIs);
                
                % --- 预计算 Masks (这是加速的关键) ---
                % 提前把 ROI 转成 0/1 矩阵，处理视频时直接点乘，速度极快
                app.ROI_Masks = cell(1, nROIs);
                for i = 1:nROIs
                    coords = app.get_coords_from_roi(rawROIs{i});
                    if ~isempty(coords)
                        app.ROI_Masks{i} = poly2mask(coords(:,1), coords(:,2), app.ImgHeight, app.ImgWidth);
                    else
                        app.ROI_Masks{i} = false(app.ImgHeight, app.ImgWidth);
                    end
                end
                
                app.LblROIStatus.Text = sprintf('ROI: %s (%d cells)', f, nROIs);
                app.LblROIStatus.FontColor = [0 0.6 0];
                app.BtnRun.Enable = 'on';
                app.log(sprintf('成功解析 %d 个 ROI 并生成掩膜。', nROIs));
                
            catch ME
                uialert(app.UIFigure, ME.message, 'ROI Error');
            end
        end
        
        % === 3. 开始批处理 ===
        function onRun(app, ~, ~)
            if isempty(app.TiffFiles) || isempty(app.ROI_Masks), return; end
            
            nFiles = length(app.TiffFiles);
            nROIs = length(app.ROI_Masks);
            fps = app.SpinnerFPS.Value;
            
            FullData = []; % 存储所有数据的累加矩阵
            
            d = uiprogressdlg(app.UIFigure, 'Title', 'Batch Processing', 'Message', 'Starting...', 'Cancelable', 'on');
            
            try
                % --- 文件循环 ---
                for k = 1:nFiles
                    if d.CancelRequested, break; end
                    
                    fname = app.TiffNames{k};
                    fpath = app.TiffFiles{k};
                    d.Message = sprintf('Processing File %d/%d: %s', k, nFiles, fname);
                    d.Value = (k-1)/nFiles;
                    
                    % 1. 读取整个 Stack (使用 Tiff 库比 imread 快)
                    Stack = app.read_tiff_stack_fast(fpath);
                    [h, w, nFrames] = size(Stack);
                    
                    if h ~= app.ImgHeight || w ~= app.ImgWidth
                        app.log(['警告: ' fname ' 尺寸不匹配，跳过！']);
                        continue;
                    end
                    
                    % 2. 提取灰度 (矩阵化加速运算)
                    % 将 stack 变形为 [Pixels, Frames]
                    ReshapedStack = reshape(Stack, h*w, nFrames); 
                    ReshapedStack = double(ReshapedStack);
                    
                    FileTraces = zeros(nFrames, nROIs);
                    
                    for r = 1:nROIs
                        mask = app.ROI_Masks{r};
                        if ~any(mask(:)), continue; end
                        
                        % 利用逻辑索引快速提取均值
                        % mask(:) 是列向量，直接在 ReshapedStack 中提取对应行
                        pixel_vals = ReshapedStack(mask(:), :);
                        FileTraces(:, r) = mean(pixel_vals, 1)';
                    end
                    
                    % 3. 拼接到总数据
                    FullData = [FullData; FileTraces];
                    
                    app.log(sprintf('完成: %s (%d 帧)', fname, nFrames));
                end
                
                % --- 导出与绘图 ---
                if ~isempty(FullData)
                    d.Message = 'Saving CSV...';
                    d.Value = 1;
                    
                    totalFrames = size(FullData, 1);
                    TimeVec = (1:totalFrames)' / fps;
                    
                    % 1. 保存 CSV
                    outName = fullfile(fileparts(app.TiffFiles{1}), 'Batch_Results_Merged.csv');
                    
                    % 构建表头
                    varNames = {'Frame', 'Time'};
                    for i=1:nROIs, varNames{end+1} = sprintf('Mean_ROI_%03d', i); end
                    
                    T = array2table([(1:totalFrames)', TimeVec, FullData], 'VariableNames', varNames);
                    writetable(T, outName);
                    
                    % 2. 绘图预览 (只画前5个 ROI)
                    plot(app.UIAxes, TimeVec, FullData(:, 1:min(5, nROIs)));
                    title(app.UIAxes, 'Trace Preview (First 5 ROIs)');
                    xlabel(app.UIAxes, 'Time (s)'); ylabel(app.UIAxes, 'Raw Intensity');
                    
                    app.log('批处理完成！');
                    app.log(['结果已保存: ' outName]);
                    uialert(app.UIFigure, ['Saved to: ' outName], 'Success');
                end
                
            catch ME
                app.log(['Error: ' ME.message]);
                uialert(app.UIFigure, ME.message, 'Processing Error');
            end
            close(d);
        end
        
        % --- 辅助: ReadImageJROI 坐标提取 ---
        function coords = get_coords_from_roi(~, roi)
            coords = [];
            if isfield(roi, 'mnCoordinates') && ~isempty(roi.mnCoordinates)
                coords = roi.mnCoordinates;
            elseif isfield(roi, 'vnRectBounds')
                % 矩形或椭圆转换
                t = roi.vnRectBounds(1); l = roi.vnRectBounds(2);
                b = roi.vnRectBounds(3); r = roi.vnRectBounds(4);
                if strcmpi(roi.strType, 'Oval')
                    theta = linspace(0, 2*pi, 30)';
                    w=(r-l)/2; h=(b-t)/2; cx=l+w; cy=t+h;
                    coords = [cx+w*cos(theta), cy+h*sin(theta)];
                else % Rect
                    coords = [l, t; r, t; r, b; l, b; l, t];
                end
            end
        end
        
        % --- 辅助: 快速读取 Tiff ---
        function Stack = read_tiff_stack_fast(~, path)
            t = Tiff(path, 'r');
            info = imfinfo(path);
            h = info(1).Height; w = info(1).Width; nF = numel(info);
            Stack = zeros(h, w, nF, 'uint16'); % 假设16位，兼容性好
            for k = 1:nF
                t.setDirectory(k);
                Stack(:,:,k) = t.read();
            end
            t.close();
        end
        
        function log(app, msg)
            app.TextAreaLog.Value = [app.TextAreaLog.Value; {msg}];
            scroll(app.TextAreaLog, 'bottom');
            drawnow;
        end
    end
    
    % === 构造函数: 创建 UI ===
    methods (Access = public)
        function app = Batch_Trace_Analysis_GUI
            % Create UIFigure and Components
            app.UIFigure = uifigure('Position', [100 100 900 600], 'Name', 'Batch Trace Analysis (ImageJ ROI)');
            app.GridLayout = uigridlayout(app.UIFigure, [1 2]);
            app.GridLayout.ColumnWidth = {300, '1x'};
            
            % Left Panel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Column = 1;
            
            app.BtnSelectTiffs = uibutton(app.LeftPanel, 'Text', '1. Select TIFFs (Multi)', ...
                'Position', [20 520 260 40], 'ButtonPushedFcn', @app.onSelectTiffs);
            
            app.ListBoxFiles = uilistbox(app.LeftPanel, 'Position', [20 350 260 160]);
            
            app.BtnSelectROI = uibutton(app.LeftPanel, 'Text', '2. Select ROI Set (.zip)', ...
                'Position', [20 300 260 40], 'ButtonPushedFcn', @app.onSelectROI);
            
            app.LblROIStatus = uilabel(app.LeftPanel, 'Text', 'ROI: Not Loaded', ...
                'Position', [20 270 260 20], 'FontColor', [0.5 0.5 0.5]);
            
            app.LblFPS = uilabel(app.LeftPanel, 'Text', 'FPS:', 'Position', [20 230 40 20]);
            app.SpinnerFPS = uispinner(app.LeftPanel, 'Limits', [1 200], 'Value', 30, ...
                'Position', [70 230 100 20]);
            
            app.BtnRun = uibutton(app.LeftPanel, 'Text', '3. Run Batch Processing', ...
                'Position', [20 160 260 50], 'BackgroundColor', [0.4 0.6 1], ...
                'FontWeight', 'bold', 'Enable', 'off', 'ButtonPushedFcn', @app.onRun);
            
            app.TextAreaLog = uitextarea(app.LeftPanel, 'Position', [20 20 260 120]);
            
            % Right Panel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Column = 2;
            
            app.UIAxes = uiaxes(app.RightPanel, 'Position', [20 20 520 540]);
            title(app.UIAxes, 'Preview');
        end
    end
end