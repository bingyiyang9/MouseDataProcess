function PlotROI_TriColor_Logic()
    %% 0. 环境检查
    if ~exist('ReadImageJROI', 'file')
        errordlg('请确保 ReadImageJROI 工具箱已在路径中。');
        return;
    end
    clc; clear; close all;

    % 过滤参数
    MaxAspectRatio = 6.0; 
    MinArea = 10;         

    %% 1. 用户输入阈值
    prompt = {'输入重叠阈值 (0.1-1.0，按较小ROI面积计算):'};
    dlgtitle = '共定位判定参数';
    definput = {'0.5'}; 
    answer = inputdlg(prompt, dlgtitle, [1 50], definput);
    if isempty(answer), return; end
    OverlapThreshold = str2double(answer{1});

    %% 2. 文件选择
    [f9, p9] = uigetfile({'*.zip;*.roi'}, '1. 选择 920nm (Green) ROI 文件');
    if isequal(f9,0), return; end
    
    [f10, p10] = uigetfile({'*.zip;*.roi'}, '2. 选择 1030nm (Red) ROI 文件');
    if isequal(f10,0), return; end

    %% 3. 解析与预处理
    fprintf('正在解析 ROI 数据...\n');
    roi9_raw = ReadImageJROI(fullfile(p9, f9));
    roi10_raw = ReadImageJROI(fullfile(p10, f10));
    if isstruct(roi9_raw), roi9_raw = {roi9_raw}; end
    if isstruct(roi10_raw), roi10_raw = {roi10_raw}; end

    [list9, max9] = extract_coords_helper(roi9_raw, MaxAspectRatio, MinArea);
    [list10, max10] = extract_coords_helper(roi10_raw, MaxAspectRatio, MinArea);
    
    % 确定画布尺寸
    canvasWH = max([max9; max10]) + 50; 
    W = ceil(canvasWH(1)); H = ceil(canvasWH(2));

    num9 = length(list9);
    num10 = length(list10);
    
    % 预计算 Mask 和 面积
    masks9 = cell(1, num9); areas9 = zeros(1, num9);
    for i = 1:num9
        masks9{i} = poly2mask(list9{i}(:,1), list9{i}(1:end,2), H, W);
        areas9(i) = sum(masks9{i}(:));
    end
    
    masks10 = cell(1, num10); areas10 = zeros(1, num10);
    for j = 1:num10
        masks10{j} = poly2mask(list10{j}(:,1), list10{j}(:,2), H, W);
        areas10(j) = sum(masks10{j}(:));
    end

    %% 4. 核心三色分类逻辑
    fprintf('正在进行空间拓扑分析...\n');
    
    % 标记哪些 ROI 已经归为"黄色共定位"组
    matched9 = false(1, num9);
    matched10 = false(1, num10);
    yellow_list = {}; 

    for i = 1:num9
        m9 = masks9{i};
        a9 = areas9(i);
        
        for j = 1:num10
            if matched10(j), continue; end % 如果1030已被占用，跳过
            
            m10 = masks10{j};
            a10 = areas10(j);
            
            % 计算交集
            intersectArea = sum(m9(:) & m10(:));
            
            if intersectArea > 0
                % 计算比例：交集 / 较小的那个ROI面积
                overlapRatio = intersectArea / min(a9, a10);
                
                if overlapRatio >= OverlapThreshold
                    % 判定为共定位
                    matched9(i) = true;
                    matched10(j) = true;
                    
                    % 选取两者中面积较大的形状作为黄色代表
                    if a9 >= a10
                        yellow_list{end+1} = list9{i};
                    else
                        yellow_list{end+1} = list10{j};
                    end
                    break; % 920(i) 已经处理完毕
                end
            end
        end
    end

    %% 5. 绘图 (白色背景，互斥绘制)
    figure('Name', 'Neuron Map: Exclusive Tri-Color', 'Color', 'w', 'Position', [100, 100, 850, 850]);
    hold on; axis image;
    % 设置坐标轴：白色背景，黑轴，Y轴反向
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse');
    box on;

    % 颜色配置
    Color920Only  = [0, 0.5, 0];  % 深绿色
    Color1030Only = [1, 0, 0];    % 红色
    ColorColoc    = [1, 0.9, 0];  % 黄色

    % A. 绘制 920 Only (深绿)
    for i = find(~matched9)
        patch(list9{i}(:,1), list9{i}(:,2), Color920Only, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end

    % B. 绘制 1030 Only (红色) - 此时已自动排除了与920重叠的部分
    for j = find(~matched10)
        patch(list10{j}(:,1), list10{j}(:,2), Color1030Only, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end

    % C. 绘制 共定位 (黄色)
    for k = 1:length(yellow_list)
        patch(yellow_list{k}(:,1), yellow_list{k}(:,2), ColorColoc, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    end

    % 界面美化
    title(['Neuron Dist: DeepGreen(920 Only), Red(1030 Only), Yellow(Co-localized, n=', ...
           num2str(length(yellow_list)), ')'], 'FontSize', 12);
    xlabel('X (px)'); ylabel('Y (px)');
    grid on;
    
    fprintf('分析完成：\n - 纯920神经元: %d\n - 纯1030神经元: %d\n - 共定位神经元: %d\n', ...
            sum(~matched9), sum(~matched10), length(yellow_list));
end

%% --- 坐标提取辅助函数 ---
function [coordsList, maxXY] = extract_coords_helper(roiRaw, maxAR, minArea)
    coordsList = {};
    maxXY = [0 0];
    for i = 1:length(roiRaw)
        curr = roiRaw{i};
        coords = [];
        if isfield(curr, 'mnCoordinates') && ~isempty(curr.mnCoordinates)
            coords = curr.mnCoordinates;
        elseif isfield(curr, 'vnRectBounds')
            t = curr.vnRectBounds(1); l = curr.vnRectBounds(2);
            b = curr.vnRectBounds(3); r = curr.vnRectBounds(4);
            if strcmpi(curr.strType, 'Oval')
                tt = linspace(0, 2*pi, 40)';
                w = (r-l)/2; h = (b-t)/2;
                coords = [(l+w) + w*cos(tt), (t+h) + h*sin(tt)];
            else 
                coords = [l, t; r, t; r, b; l, b; l, t];
            end
        end
        if isempty(coords), continue; end
        w_roi = max(coords(:,1)) - min(coords(:,1));
        h_roi = max(coords(:,2)) - min(coords(:,2));
        if w_roi <= 0, w_roi = 1; end; if h_roi <= 0, h_roi = 1; end
        if (max(w_roi, h_roi)/min(w_roi, h_roi) > maxAR) || (w_roi*h_roi < minArea)
            continue; 
        end
        coordsList{end+1} = coords;
        maxXY = max([maxXY; coords]);
    end
end