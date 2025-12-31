function CSV_Trajectory_Cleaner()
%% CSV_TRAJECTORY_CLEANER
%  功能：读取行为学 CSV 数据 -> 剔除剧烈震荡/跳变 -> 平滑处理 -> 导出 .fig
%
%  针对痛点：
%  1. 轨迹中有"飞线" (瞬移到远处的错误点)。
%  2. 轨迹线条锯齿严重，不平滑。
%  3. 需要生成可编辑的 .fig 文件用于论文排版。
%
%  Author: Gemini
%  Date: 2025-12-29

    clc; clear; close all;

    %% 1. 参数设置 (根据震荡程度调整)
    % ============================================================
    % --- 强力去震荡参数 ---
    cleanParams.maxJump = 5;          % [关键] 最大跳变阈值 (像素)。
                                       % 任何两帧间距离超过此值，会被视为"严重震荡"并直接删除。
                                       % 建议值：30 - 80。越小删得越狠。

    cleanParams.smoothWin = 100;        % [关键] 平滑窗口大小。
                                       % 值越大线条越圆润，但可能丢失微小转弯细节。
                                       % 建议值：10 - 20。
    % ============================================================

    %% 2. 读取 CSV 文件
    [fileName, pathName] = uigetfile({'*.csv', 'CSV Files (*.csv)'}, '选择包含轨迹的 CSV 文件');
    if isequal(fileName, 0), return; end
    fullPath = fullfile(pathName, fileName);
    
    % 读取数据
    T = readtable(fullPath);
    
    % 自动识别 X 和 Y 列 (假设在第2、3列，或者根据表头查找)
    % 尝试常见的列名
    if ismember('X', T.Properties.VariableNames) && ismember('Y', T.Properties.VariableNames)
        rawX = T.X; rawY = T.Y;
    elseif ismember('X_px', T.Properties.VariableNames)
        rawX = T.X_px; rawY = T.Y_px;
    else
        % 如果找不到名字，默认取第2和第3列
        fprintf('警告：未找到标准列名，默认使用第2列作为X，第3列作为Y。\n');
        rawX = T{:, 2};
        rawY = T{:, 3};
    end
    
    fprintf('成功读取数据: %d 行\n', length(rawX));

    %% 3. 数据清洗 (去震荡核心)
    cleanX = rawX;
    cleanY = rawY;

    % --- 步骤 A: 剔除剧烈跳变 (Outlier Removal) ---
    % 计算帧间距离
    diffX = diff(cleanX);
    diffY = diff(cleanY);
    dist = [0; sqrt(diffX.^2 + diffY.^2)];
    
    % 找到震荡点
    outlierIdx = dist > cleanParams.maxJump;
    numOutliers = sum(outlierIdx);
    
    % 将震荡点设为 NaN (删除)
    cleanX(outlierIdx) = NaN;
    cleanY(outlierIdx) = NaN;
    
    fprintf('  -> 已剔除 %d 个严重震荡点 (阈值: %d px)\n', numOutliers, cleanParams.maxJump);
    
    % --- 步骤 B: 线性插值 (补全断点) ---
    cleanX = fillmissing(cleanX, 'linear');
    cleanY = fillmissing(cleanY, 'linear');
    
    % --- 步骤 C: 平滑处理 (Smoothing) ---
    finalX = smoothdata(cleanX, 'movmean', cleanParams.smoothWin);
    finalY = smoothdata(cleanY, 'movmean', cleanParams.smoothWin);
    
    fprintf('  -> 已完成平滑处理 (窗口: %d)\n', cleanParams.smoothWin);

    %% 4. 绘图与导出
    [~, nameBody, ~] = fileparts(fileName);
    outDir = fullfile(pathName, [nameBody '_Pinghua']);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    
    % --- 图 1: 对比图 (原始 vs 清洗后) ---
    hFig1 = figure('Name', 'Comparison', 'Color', 'w');
    subplot(1,2,1);
    plot(rawX, rawY, 'Color', [0.7 0.7 0.7]); 
    axis ij; axis equal; title('原始数据 (含震荡)');
    
    subplot(1,2,2);
    plot(finalX, finalY, 'r-', 'LineWidth', 1.5);
    axis ij; axis equal; title('清洗后数据');
    
    % --- 图 2: Nature 风格最终图 (红线白底) ---
    hFig2 = figure('Name', 'Final_Nature', 'Color', 'w');
    plot(finalX, finalY, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5);
    axis ij; axis equal; axis off; % 极简风格
    % 自动调整范围
    margin = 50;
    xlim([min(finalX)-margin, max(finalX)+margin]);
    ylim([min(finalY)-margin, max(finalY)+margin]);
    title('Final Trajectory');

    % --- 图 3: 时间渐变图 (Time Gradient) ---
    hFig3 = figure('Name', 'Final_Gradient', 'Color', 'w');
    z = zeros(size(finalX));
    col = 1:length(finalX);
    patch([finalX;NaN], [finalY;NaN], [z;NaN], [col(:);NaN], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2);
    colormap(jet); colorbar;
    axis ij; axis equal; axis off;
    title('Time-Coded Trajectory');

    %% 5. 保存 .fig 和 .png
    fprintf('正在保存文件...\n');
    
    % 保存 Nature 风格
    savefig(hFig2, fullfile(outDir, [nameBody '_Nature.fig']));
    print(hFig2, fullfile(outDir, [nameBody '_Nature.png']), '-dpng', '-r300');
    
    % 保存 渐变风格
    savefig(hFig3, fullfile(outDir, [nameBody '_Gradient.fig']));
    print(hFig3, fullfile(outDir, [nameBody '_Gradient.png']), '-dpng', '-r300');
    
    % 导出一份清洗后的新 CSV
    T_new = table((1:length(finalX))', finalX, finalY, ...
        'VariableNames', {'Frame', 'X_Clean', 'Y_Clean'});
    writetable(T_new, fullfile(outDir, [nameBody '_CleanedData.csv']));
    
    fprintf('✅ 全部完成！结果保存在: %s\n', outDir);
    
    msgbox(['处理完成！已剔除 ' num2str(numOutliers) ' 个震荡点。'], '成功');
end