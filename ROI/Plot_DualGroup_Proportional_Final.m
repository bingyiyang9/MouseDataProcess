function Plot_DualGroup_Proportional_Final
    %% 1. 初始化与参数设置
    clc; clear; close all;
    
    % --- 核心参数 ---
    FixedBgVal = 26500;       % 背景扣除值
    CLim = [0, 4];            % 热图色阶范围
    NumTraces = 6;            % 右侧随机展示的波形数量
    % ----------------
    
    % 交互输入参数
    prompt = {'视频帧率 (FPS):', '基线滑动窗口 (秒):'};
    definput = {'8.4', '60'}; % 默认 8.4 FPS, 120秒滑动窗口
    answer = inputdlg(prompt, '参数设置', [1 40], definput);
    if isempty(answer), return; end
    
    FPS = str2double(answer{1});
    WindowSeconds = str2double(answer{2});
    WinSize = round(WindowSeconds * FPS); % 转换为帧数

    %% 2. 加载数据
    % Group 1: 920 Only (Green) - 放下面
    disp('Step 1: 请选择 [920 only] 的 CSV (Green)...');
    [f1, p1] = uigetfile('*.csv', 'CSV: 920 Only');
    if isequal(f1, 0), return; end
    [dFF1, Time1] = process_csv(fullfile(p1, f1), FixedBgVal, FPS, WinSize);
    
    % Group 2: 920 & 1030 (Dual) - 放上面
    disp('Step 2: 请选择 [920 & 1030] 的 CSV (Dual)...');
    [f2, p2] = uigetfile('*.csv', 'CSV: 920 & 1030');
    if isequal(f2, 0), return; end
    [dFF2, Time2] = process_csv(fullfile(p2, f2), FixedBgVal, FPS, WinSize);

    %% 3. 计算动态高度比例
    n1 = size(dFF1, 2); % Group 1 细胞数
    n2 = size(dFF2, 2); % Group 2 细胞数
    TotalN = n1 + n2;
    
    % 定义绘图区域的上下边距 (归一化坐标 0-1)
    MarginBottom = 0.12; % 留出空间给 X 轴标签
    MarginTop = 0.08;
    Gap = 0.02; 
    AvailableHeight = 1 - MarginBottom - MarginTop - Gap;
    
    % 动态分配高度 (保证最小高度，防止压扁)
    MinHeight = 0.05; 
    H2 = max(MinHeight, AvailableHeight * (n2 / TotalN)); % 上图高度
    H1 = AvailableHeight - H2;                            % 下图高度
    
    Y1 = MarginBottom;              % 下图起点
    Y2 = MarginBottom + H1 + Gap;   % 上图起点

    %% 4. 绘图
    hFig = figure('Name', 'Dual Group Analysis (Minutes)', 'Color', 'w', 'Position', [50, 50, 1400, 900]);
    
    % 左侧布局
    LeftMargin = 0.08;
    WidthHeatmap = 0.60;
    
    % --- 下图: Group 1 (920 Only) ---
    ax1 = axes('Position', [LeftMargin, Y1, WidthHeatmap, H1]);
    imagesc(Time1, 1:n1, dFF1');
    colormap(ax1, 'parula'); clim(ax1, CLim);
    ylabel('Neuron #', 'FontWeight', 'bold');
    xlabel('Time (min)', 'FontWeight', 'bold', 'FontSize', 11); % 横轴改为分钟
    title(sprintf('Group 1: 920 Only (n=%d)', n1));
    
    % --- 上图: Group 2 (Dual) ---
    ax2 = axes('Position', [LeftMargin, Y2, WidthHeatmap, H2]);
    imagesc(Time2, 1:n2, dFF2');
    colormap(ax2, 'parula'); clim(ax2, CLim);
    ylabel('Neuron #', 'FontWeight', 'bold');
    title(sprintf('Group 2: 920 & 1030 (n=%d)', n2));
    set(gca, 'XTickLabel', []); % 隐藏上图 X 轴标尺
    
    % Colorbar
    cb = colorbar;
    cb.Position = [LeftMargin+WidthHeatmap+0.01, Y1, 0.02, H1+H2+Gap];
    cb.Label.String = '\DeltaF/F';
    
    linkaxes([ax1, ax2], 'x'); % 缩放同步
    
    %% 5. 右侧：波形图 (随机抽取)
    
    % --- [修改点] 随机抽取索引 ---
    % Group 1 随机
    if n1 <= NumTraces
        rand1 = 1:n1;
    else
        rand1 = sort(randperm(n1, NumTraces));
    end
    
    % Group 2 随机
    if n2 <= NumTraces
        rand2 = 1:n2;
    else
        rand2 = sort(randperm(n2, NumTraces));
    end
    
    % 坐标设置
    ax3 = axes('Position', [0.78, MarginBottom, 0.20, 1-MarginBottom-MarginTop]);
    hold on;
    
    offset = 0;
    StackSpace = 5.0; % 波形间距
    
    % 画 Group 1 (绿色)
    c1 = [0 0.6 0.2];
    for i = 1:length(rand1)
        plot(Time1, dFF1(:, rand1(i)) + offset, 'Color', c1, 'LineWidth', 0.8);
        text(Time1(end), offset+0.5, sprintf(' G1-%d', rand1(i)), 'Color', c1, 'FontSize',8);
        offset = offset + StackSpace;
    end
    
    % 分隔线
    yline(offset + StackSpace/2, 'k:', 'LineWidth', 1.5);
    offset = offset + StackSpace * 1.5;
    
    % 画 Group 2 (紫红)
    c2 = [0.8 0 0.4];
    for i = 1:length(rand2)
        plot(Time2, dFF2(:, rand2(i)) + offset, 'Color', c2, 'LineWidth', 0.8);
        text(Time2(end), offset+0.5, sprintf(' G2-%d', rand2(i)), 'Color', c2, 'FontSize',8);
        offset = offset + StackSpace;
    end
    
    title('Random Sample Traces');
    xlabel('Time (min)', 'FontWeight', 'bold'); % 横轴为分钟
    axis tight; axis off;
    ylim([-1, offset + 2]);
    
    fprintf('绘图完成！横轴单位：分钟。右侧为随机抽样。\n');
end

%% === 数据处理函数 ===
function [dFF, Time] = process_csv(filepath, bg_val, fps, win_size)
    opts = detectImportOptions(filepath);
    opts.VariableNamingRule = 'preserve';
    T = readtable(filepath, opts);
    
    dataCols = [];
    varNames = T.Properties.VariableNames;
    for i = 1:length(varNames)
        if (startsWith(varNames{i}, 'Mean') || startsWith(varNames{i}, 'ROI')) && isnumeric(T.(varNames{i}))
            dataCols = [dataCols, i];
        end
    end
    RawF = table2array(T(:, dataCols));
    [nFrames, nROIs] = size(RawF);
    
    % [修改点] 时间轴转换为分钟
    Time = (1:nFrames)' / fps / 60; 
    
    dFF = zeros(size(RawF));
    
    for i = 1:nROIs
        F = RawF(:, i);
        F_corr = F - bg_val; 
        F_smooth = smoothdata(F_corr, 'movmean', 5);
        
        % [滑动窗口说明]
        % 这里的 win_size 就是滑动窗口的长度 (以帧为单位)
        % movmin 会在当前点的前后各取半个窗口，找最小值作为 F0
        if exist('movprctile', 'file')
            F0 = movprctile(F_smooth, win_size, 10);
        else
            F0 = movmin(F_smooth, win_size);
            F0 = smoothdata(F0, 'gaussian', fps*10); % 对基线做一点平滑
        end
        
        F0(F0 < 1) = 1; 
        dFF(:, i) = (F_smooth - F0) ./ F0;
    end
end