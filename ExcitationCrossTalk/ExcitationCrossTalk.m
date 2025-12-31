%% 双光子激发串扰定标分析 (v7: 多Ref平均 + ID追踪 + 纯最小值基线)
% Author: AI Expert
% Date: 2025-12-22

clear; clc; close all;

%% 1. 参数与路径设置
base_path = 'C:\Users\yangbingyi\Downloads\crosstalk_R_2X\crosstalk_R_2X'; 

% 功率点列表
power_list = [1, 1.5, 2, 2.5, 3, 3.5,4]; 

% 固定的子路径
sub_path = fullfile('CellVideo2', 'CellVideo', 'CellVideo 1.tif');

fprintf('=== 实验配置 ===\n');

% --- (A) 寻找并处理 Dark (背景) ---
dark_struct = dir(fullfile(base_path, '*Crosstalk_R_dark*'));
if isempty(dark_struct), error('找不到 Dark 文件夹'); end
path_dark = fullfile(dark_struct(1).folder, dark_struct(1).name, sub_path);
fprintf('Dark File: %s\n', dark_struct(1).name);

% 计算基线: 帧平均 -> 全局最小值 (无中值滤波)
bg_raw_img = load_tiff_avg(path_dark);
bg_val = min(bg_raw_img(:)); 
fprintf('Background Baseline (Global Min): %.2f\n', bg_val);

% --- (B) 寻找并平均多个 Ref (1030nm) ---
ref_struct = dir(fullfile(base_path, '*Crosstalk_R_Ref*'));
if isempty(ref_struct), error('找不到 Ref 1030nm 文件夹'); end

fprintf('Found %d Ref folders. Averaging them...\n', length(ref_struct));
img_ref_accum = 0;
for i = 1:length(ref_struct)
    p_ref = fullfile(ref_struct(i).folder, ref_struct(i).name, sub_path);
    fprintf('  - Ref %d: %s\n', i, ref_struct(i).name);
    img_ref_accum = img_ref_accum + load_tiff_avg(p_ref);
end
% 计算 Ref 平均图并扣除背景
img_ref_master = (img_ref_accum / length(ref_struct)) - bg_val;
img_ref_master(img_ref_master < 0) = 0; % 归零

%% 2. [交互式] 定义 Master ROIs 并生成 ID 图
fprintf('\nStep 1: 设定 Reference 母板 (请框选 Beads)...\n');

% 弹出选圆 GUI
[ref_centers, ref_radii] = gui_roi_selector_interactive(img_ref_master, 'Master Ref: Averaged Image');

if isempty(ref_centers)
    error('未选中任何 Beads，程序终止');
end

num_beads = length(ref_radii);
bead_ids = (1:num_beads)'; % 生成 ID: 1, 2, 3...

% 计算 Reference 的强度 (分母)
ref_intensities = extract_intensity(img_ref_master, ref_centers, ref_radii);

% --- 输出图1: 整合的 ID 对比图 ---
f_map = figure('Name', 'Bead ID Map', 'Color', 'w');
imshow(img_ref_master, []); 
title(['Master Reference ID Map (n=' num2str(num_beads) ')']); hold on;
viscircles(ref_centers, ref_radii, 'Color', 'g', 'LineWidth', 1);
% 标注 ID
for i = 1:num_beads
    text(ref_centers(i,1), ref_centers(i,2), num2str(bead_ids(i)), ...
        'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hold off;
saveas(f_map, 'Crosstalk_Bead_ID_Map.png');
fprintf('已保存 ID 映射图: Crosstalk_Bead_ID_Map.png\n');

%% 3. 循环处理所有功率点 (Data Extraction)
% 初始化数据容器
% Columns: Method(1=DUET, 2=Trad), Power, BeadID, Val_920, Val_Ref, Ratio
final_data_matrix = []; 

for p = power_list
    if mod(p, 1) == 0, p_str = num2str(p); else, p_str = num2str(p); end
    
    fprintf('\n----------------------------------------\n');
    fprintf('正在分析功率点: %s mW\n', p_str);
    
    pattern_duet = ['*Crosstalk_DUET_R_' p_str 'mW_*'];
    pattern_trad = ['*Crosstalk_R_' p_str 'mW_*'];
    
    % === 处理 DUET ===
    d_list = dir(fullfile(base_path, pattern_duet));
    if ~isempty(d_list)
        f_name = d_list(1).name;
        fprintf('  > DUET: %s\n', f_name);
        path_img = fullfile(d_list(1).folder, f_name, sub_path);
        
        [img_curr, ~] = load_preprocess(path_img, bg_val);
        
        % 配准
        shifted_c = gui_registration_cross(img_curr, ref_centers, ref_radii, bead_ids, ['Align DUET: ' p_str 'mW']);
        
        % 提取数据
        vals = extract_intensity(img_curr, shifted_c, ref_radii);
        ratios = vals ./ ref_intensities;
        
        % 存入矩阵: [Method=1, Power, ID, Val920, ValRef, Ratio]
        current_data = [ones(num_beads,1), repmat(p, num_beads,1), bead_ids, vals, ref_intensities, ratios];
        final_data_matrix = [final_data_matrix; current_data];
    else
        warning('Missing DUET data for %s mW', p_str);
    end
    
    % === 处理 Traditional ===
    t_list = dir(fullfile(base_path, pattern_trad));
    % 排除 DUET
    valid_idx = 0;
    for k=1:length(t_list)
        if ~contains(t_list(k).name, 'DUET'), valid_idx = k; break; end
    end
    
    if valid_idx > 0
        f_name = t_list(valid_idx).name;
        fprintf('  > Trad: %s\n', f_name);
        path_img = fullfile(t_list(valid_idx).folder, f_name, sub_path);
        
        [img_curr, ~] = load_preprocess(path_img, bg_val);
        
        % 配准
        shifted_c = gui_registration_cross(img_curr, ref_centers, ref_radii, bead_ids, ['Align Trad: ' p_str 'mW']);
        
        % 提取数据
        vals = extract_intensity(img_curr, shifted_c, ref_radii);
        ratios = vals ./ ref_intensities;
        
        % 存入矩阵: [Method=2, Power, ID, Val920, ValRef, Ratio]
        current_data = [ones(num_beads,1)*2, repmat(p, num_beads,1), bead_ids, vals, ref_intensities, ratios];
        final_data_matrix = [final_data_matrix; current_data];
    else
        warning('Missing Trad data for %s mW', p_str);
    end
end

%% 4. 生成 CSV 表格
fprintf('\nStep 3: 导出 CSV 数据...\n');

Method_Strs = cell(size(final_data_matrix, 1), 1);
for i = 1:length(Method_Strs)
    if final_data_matrix(i, 1) == 1
        Method_Strs{i} = 'With DUET';
    else
        Method_Strs{i} = 'Traditional';
    end
end

T = table(Method_Strs, final_data_matrix(:,2), final_data_matrix(:,3), ...
    final_data_matrix(:,4), final_data_matrix(:,5), final_data_matrix(:,6), ...
    'VariableNames', {'Method', 'Power_mW', 'Bead_ID', 'Intensity_920', 'Intensity_Ref_1030', 'Crosstalk_Ratio'});

disp(head(T));
writetable(T, 'Crosstalk_Excitation_Detail.csv');
fprintf('表格已保存: Crosstalk_Excitation_Detail.csv\n');

%% 5. 绘图与拟合
fprintf('Step 4: 生成趋势图...\n');
visualize_final_trend(final_data_matrix);


%% ================== 函数库 ==================

% --- 1. 交互式 Ref 选圆 (保留v6逻辑) ---
function [final_centers, final_radii] = gui_roi_selector_interactive(img, win_title)
    hFig = figure('Name', win_title, 'NumberTitle', 'off', 'MenuBar', 'none', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
    hAx = axes('Parent', hFig, 'Position', [0.05, 0.1, 0.7, 0.85]);
    imshow(img, [], 'Parent', hAx); title(hAx, '请调整滑块以圈中 Ref 图上所有的 Beads'); hold(hAx, 'on');
    
    uicontrol('Style','text','Units','normalized','Position',[0.8,0.85,0.15,0.03],'String','Min Radius:','FontWeight','bold');
    hSMin = uicontrol('Style','slider','Units','normalized','Position',[0.8,0.82,0.15,0.03],'Min',3,'Max',50,'Value',6);
    hEMin = uicontrol('Style','edit','Units','normalized','Position',[0.955,0.82,0.04,0.03],'String','6');
    uicontrol('Style','text','Units','normalized','Position',[0.8,0.75,0.15,0.03],'String','Max Radius:','FontWeight','bold');
    hSMax = uicontrol('Style','slider','Units','normalized','Position',[0.8,0.72,0.15,0.03],'Min',5,'Max',100,'Value',30);
    hEMax = uicontrol('Style','edit','Units','normalized','Position',[0.955,0.72,0.04,0.03],'String','30');
    uicontrol('Style','text','Units','normalized','Position',[0.8,0.65,0.15,0.03],'String','Sensitivity:','FontWeight','bold');
    hSSens = uicontrol('Style','slider','Units','normalized','Position',[0.8,0.62,0.15,0.03],'Min',0.6,'Max',0.99,'Value',0.90);
    hESens = uicontrol('Style','edit','Units','normalized','Position',[0.955,0.62,0.04,0.03],'String','0.90');
    hText = uicontrol('Style','text','Units','normalized','Position',[0.8,0.5,0.15,0.05],'String','Found: 0','FontSize',12,'ForegroundColor','b');
    uicontrol('Style','pushbutton','Units','normalized','Position',[0.8,0.1,0.15,0.08],'String','Confirm & Save','FontSize',14,'BackgroundColor','g','Callback',@(s,e) uiresume(hFig));
    
    % 简单的归一化辅助显示
    img_norm = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    gui_data.img = img_norm; gui_data.hAx = hAx; gui_data.hText = hText;
    
    cb = {@update_circles, hSMin, hEMin, hSMax, hEMax, hSSens, hESens, gui_data};
    set(hSMin,'Callback',cb); set(hSMax,'Callback',cb); set(hSSens,'Callback',cb);
    set(hEMin,'Callback',cb); set(hEMax,'Callback',cb); set(hESens,'Callback',cb);
    update_circles([],[],hSMin, hEMin, hSMax, hEMax, hSSens, hESens, gui_data);
    uiwait(hFig);
    rmin=round(get(hSMin,'Value')); rmax=round(get(hSMax,'Value')); sens=get(hSSens,'Value'); if rmin>=rmax, rmax=rmin+1; end
    [final_centers, final_radii] = imfindcircles(gui_data.img, [rmin rmax], 'ObjectPolarity','bright','Sensitivity',sens);
    close(hFig);
    
    function update_circles(~,~,hSMin, hEMin, hSMax, hEMax, hSSens, hESens, gdata)
        vmin=round(get(hSMin,'Value')); set(hEMin,'String',num2str(vmin)); vmax=round(get(hSMax,'Value')); set(hEMax,'String',num2str(vmax));
        vsens=get(hSSens,'Value'); set(hESens,'String',num2str(vsens,'%.2f'));
        if vmin>=vmax, vmax=vmin+1; set(hSMax,'Value',vmax); set(hEMax,'String',num2str(vmax)); end
        try delete(findobj(gdata.hAx,'Type','line')); [c,r] = imfindcircles(gdata.img, [vmin vmax], 'ObjectPolarity','bright','Sensitivity',vsens);
            if ~isempty(c), viscircles(gdata.hAx, c, r, 'Color','g'); end; set(gdata.hText, 'String', sprintf('Found: %d', length(r)));
        catch; end
    end
end

% --- 2. 十字滑块配准 GUI (含ID显示) ---
function new_centers = gui_registration_cross(img, centers_orig, radii, ids, title_str)
    hFig = figure('Name', title_str, 'NumberTitle', 'off', 'MenuBar', 'none', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.7]);
    hAx = axes('Parent', hFig, 'Position', [0.1, 0.25, 0.75, 0.7]);
    imshow(img, [], 'Parent', hAx); title(['[配准] ' title_str]); hold on;
    
    % 绘制初始圆和ID
    draw_circles(hAx, centers_orig, radii, ids);
    
    uicontrol('Style','text','Units','normalized','Position',[0.1, 0.12, 0.75, 0.03],'String','X Shift', 'FontWeight','bold');
    hSlX = uicontrol('Style','slider','Units','normalized','Position',[0.1, 0.08, 0.75, 0.04],'Min',-50, 'Max',50, 'Value',0, 'Callback', @update_shift);
    uicontrol('Style','text','Units','normalized','Position',[0.88, 0.25, 0.1, 0.03],'String','Y Shift', 'FontWeight','bold');
    hSlY = uicontrol('Style','slider','Units','normalized','Position',[0.90, 0.25, 0.04, 0.7],'Min',-50, 'Max',50, 'Value',0, 'Callback', @update_shift);
    uicontrol('Style','pushbutton','Units','normalized','Position',[0.4, 0.01, 0.2, 0.05],'String','Confirm','FontSize',12, 'BackgroundColor','g','Callback', @(s,e) uiresume(hFig));
    
    uiwait(hFig);
    sx = get(hSlX, 'Value'); sy = get(hSlY, 'Value');
    new_centers = centers_orig + [sx, -sy];
    close(hFig);
    
    function update_shift(~,~)
        sx = get(hSlX, 'Value'); sy = get(hSlY, 'Value');
        cla(hAx); imshow(img, [], 'Parent', hAx); hold(hAx, 'on');
        draw_circles(hAx, centers_orig + [sx, -sy], radii, ids);
    end

    function draw_circles(ax, centers, r, id_list)
        viscircles(ax, centers, r, 'Color', 'g', 'LineWidth', 1);
        for i=1:length(id_list)
            text(ax, centers(i,1), centers(i,2), num2str(id_list(i)), 'Color', 'y', 'FontWeight', 'bold', 'HorizontalAlignment','center');
        end
    end
end

% --- 辅助函数 ---
function avg_img = load_tiff_avg(filepath)
    if ~isfile(filepath), error('File not found: %s', filepath); end
    info = imfinfo(filepath);
    sum_img = zeros(info(1).Height, info(1).Width);
    for k = 1:numel(info), sum_img = sum_img + double(imread(filepath, k)); end
    avg_img = sum_img / numel(info);
end

function [img_out, avg_raw] = load_preprocess(filepath, bg_val)
    avg_raw = load_tiff_avg(filepath);
    img_sub = avg_raw - bg_val;
    img_sub(img_sub < 0) = 0; 
    img_out = img_sub;
end

function intensities = extract_intensity(img, centers, radii)
    intensities = zeros(length(radii), 1);
    [rows, cols] = size(img);
    [X, Y] = meshgrid(1:cols, 1:rows);
    for i = 1:length(radii)
        mask = (X - centers(i, 1)).^2 + (Y - centers(i, 2)).^2 <= radii(i)^2;
        intensities(i) = mean(img(mask));
    end
end
function visualize_final_trend(data)
    % Data cols: [Method, Power, ID, Val920, ValRef, Ratio]
    f = figure('Color', 'w', 'Position', [100, 100, 700, 500]); hold on;
    
    % 分离数据
    mask_duet = (data(:,1) == 1);
    mask_trad = (data(:,1) == 2);
    
    duet_pow = data(mask_duet, 2); duet_rat = data(mask_duet, 6);
    trad_pow = data(mask_trad, 2); trad_rat = data(mask_trad, 6);
    
    % 绘制散点
    s1 = scatter(trad_pow, trad_rat, 50, [0.8 0.4 0.1], 'filled', 'MarkerFaceAlpha', 0.6);
    s2 = scatter(duet_pow, duet_rat, 50, [0.1 0.4 0.6], 'filled', 'MarkerFaceAlpha', 0.6);
    
    x_fit = linspace(min(data(:,2)), max(data(:,2)), 100);
    
    fprintf('\n========== 拟合结果报告 (Polynomial Degree = 2) ==========\n');
    
    % --- Traditional 拟合 ---
    if ~isempty(trad_pow)
        [p1, S1] = polyfit(trad_pow, trad_rat, 2); % p1 是系数向量 [a, b, c]
        [y_fit_trad, delta1] = polyval(p1, trad_pow, S1);
        
        % 计算 R-square
        yresid = trad_rat - y_fit_trad;
        SSresid = sum(yresid.^2);
        SStotal = (length(trad_rat)-1) * var(trad_rat);
        rsq_trad = 1 - SSresid/SStotal;
        
        % 画线
        plot(x_fit, polyval(p1, x_fit), '--', 'Color', [0.8 0.4 0.1], 'LineWidth', 2);
        
        % 输出公式
        fprintf('>> Traditional Method:\n');
        fprintf('   Equation: y = (%.4f) * x^2 + (%.4f) * x + (%.4f)\n', p1(1), p1(2), p1(3));
        fprintf('   R-squared: %.4f\n', rsq_trad);
    end
    
    % --- DUET 拟合 ---
    if ~isempty(duet_pow)
        [p2, S2] = polyfit(duet_pow, duet_rat, 2);
        [y_fit_duet, delta2] = polyval(p2, duet_pow, S2);
        
        % 计算 R-square
        yresid = duet_rat - y_fit_duet;
        SSresid = sum(yresid.^2);
        SStotal = (length(duet_rat)-1) * var(duet_rat);
        rsq_duet = 1 - SSresid/SStotal;
        
        % 画线
        plot(x_fit, polyval(p2, x_fit), '--', 'Color', [0.1 0.4 0.6], 'LineWidth', 2);
        
        % 输出公式
        fprintf('>> With DUET:\n');
        fprintf('   Equation: y = (%.4f) * x^2 + (%.4f) * x + (%.4f)\n', p2(1), p2(2), p2(3));
        fprintf('   R-squared: %.4f\n', rsq_duet);
    end
    fprintf('==========================================================\n');
    
    xlabel('920nm Power (mW)', 'FontSize', 12);
    ylabel('Crosstalk Ratio (Signal_{920} / Signal_{Ref})', 'FontSize', 12);
    legend([s2, s1], {'With DUET', 'Traditional'}, 'Location', 'northwest');
    title('Excitation Crosstalk Calibration');
    grid on; box on;
    hold off;
    
    saveas(f, 'Crosstalk_Trend_Plot.png');
    fprintf('已保存趋势图: Crosstalk_Trend_Plot.png\n');
end