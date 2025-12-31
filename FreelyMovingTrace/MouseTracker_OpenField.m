function MouseTracker_Classic_Pro()
%% MOUSETRACKER_CLASSIC_PRO
%  【回归经典+全能导出版】行为学轨迹追踪系统
%  
%  核心逻辑：
%  1. [算法回滚] 使用"背景减除法" (Max Projection)，识别最稳健，不丢失老鼠细节。
%  2. [参数优化] 适度过滤线缆 (半径4px)，允许快跑 (跳变100px)。
%  3. [全能导出] 保留 Batch/GIF/Fig源文件/Nature图/渐变图/汇总图。
%
%  Author: Gemini
%  Date: 2025-12-29

    clc; clear; close all;

    %% 1. 参数设置 (回归稳健参数)
    % ============================================================
    % --- 图像识别 (背景减除法) ---
    sysParams.diffThresh = 100/255;     % 差异阈值 (越小越灵敏)。
                                       % 只要比背景黑一点点就会被识别。
    sysParams.erodeRadius = 6;         % 【回调】半径4px (直径8px)。
                                       % 足够擦除普通数据线，但绝对不会擦掉老鼠。
    sysParams.minArea = 180;           % 最小面积
    
    % --- 运动学 (放宽限制) ---
    sysParams.maxJump = 5;           % 【回调】允许帧间移动5px。
                                       % 防止老鼠快跑时被当成坏点误删。
    sysParams.smoothWin = 15;          
    
    % --- 导出设置 (保持不变) ---
    exportParams.saveGIF = false;       % 存 GIF
    exportParams.gifFrameSkip = 5;     % GIF 抽帧
    
    visParams.showRealTime = true;     
    visParams.drawInterval = 5;        
    % ============================================================

    %% 2. 批量文件选择
    [fileNames, pathName] = uigetfile({'*.avi;*.mp4;*.mov', 'Video Files'}, ...
                                      '选择视频 (支持多选)', 'MultiSelect', 'on');
    if isequal(fileNames, 0), return; end
    if ischar(fileNames), fileNames = {fileNames}; end
    
    numFiles = length(fileNames);
    fprintf('共选中 %d 个视频，准备开始...\n', numFiles);

    allTracksData = cell(numFiles, 1);
    allVideoNames = cell(numFiles, 1);
    commonRoiPos = [];
    commonBgImg = []; % 存储通用的背景图(如果场景不变)

    %% 3. 批处理循环
    for fIdx = 1:numFiles
        currFileName = fileNames{fIdx};
        fullPath = fullfile(pathName, currFileName);
        
        fprintf('\n[%d/%d] 处理中: %s\n', fIdx, numFiles, currFileName);
        
        try
            % --- 核心处理 ---
            [smoothTraj, finalRoi, bgImg] = processSingleVideo(fullPath, sysParams, visParams, exportParams, fIdx, numFiles, commonBgImg);
            
            % --- 记录数据 ---
            allTracksData{fIdx} = smoothTraj;
            allVideoNames{fIdx} = currFileName;
            
            if fIdx == 1
                commonRoiPos = finalRoi; 
                % 如果摄像机不动，可以沿用计算好的背景，大大加速
                % commonBgImg = bgImg; 
            end
            
        catch ME
            fprintf('  ❌ 错误: %s\n', ME.message);
        end
    end
    
    %% 4. 生成汇总图 (Summary Plot)
    if ~isempty(commonRoiPos)
        fprintf('\n正在生成汇总图 (All-in-One)...\n');
        hSum = figure('Name', 'Summary', 'Color', 'w', 'Visible', 'off');
        hold on;
        
        % 画 ROI 框
        rectangle('Position', commonRoiPos, 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
        
        % 叠加所有轨迹
        colors = jet(numFiles);
        for k = 1:numFiles
            traj = allTracksData{k};
            if ~isempty(traj)
                plot(traj(:,1), traj(:,2), 'Color', [colors(k,:) 0.6], 'LineWidth', 1.5, 'DisplayName', allVideoNames{k});
            end
        end
        
        axis ij; axis equal; axis off;
        xlim([commonRoiPos(1)-50, commonRoiPos(1)+commonRoiPos(3)+50]);
        ylim([commonRoiPos(2)-50, commonRoiPos(2)+commonRoiPos(4)+50]);
        title('Summary of All Tracks');
        
        % 保存
        saveas(hSum, fullfile(pathName, 'All_Trajectories_Summary.png'));
        savefig(hSum, fullfile(pathName, 'All_Trajectories_Summary.fig'));
        fprintf('🎉 汇总图已保存\n');
        close(hSum);
    end
    
    msgbox('所有处理完成！', 'Success');


    %% --- 内部核心处理函数 ---
    function [finalTraj, roiOut, bgImage] = processSingleVideo(videoPath, sys, vis, expP, currentIdx, totalFiles, inputBg)
        [saveDir, fName, ~] = fileparts(videoPath);
        outDir = fullfile(saveDir, [fName '_Results']);
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        
        videoObj = VideoReader(videoPath);
        nFrames = videoObj.NumFrames;
        
        % --- 步骤A: 自动背景建模 (回到老方案的核心) ---
        if isempty(inputBg)
            fprintf('  -> 正在计算背景模型...\n');
            % 随机采50帧取最大值，生成纯净背景
            sampleIdx = floor(linspace(1, nFrames, 50));
            stack = [];
            for i = sampleIdx
                videoObj.CurrentTime = (i-1)/videoObj.FrameRate;
                f = readFrame(videoObj);
                if size(f,3)==3, f=rgb2gray(f); end
                stack = cat(3, stack, f);
            end
            bgImage = max(stack, [], 3); % 白底黑鼠用max，黑底白鼠用min
        else
            bgImage = inputBg;
            fprintf('  -> 沿用已有背景模型\n');
        end
        
        videoObj.CurrentTime = 0;
        frame1 = readFrame(videoObj);
        
        % --- 步骤B: ROI 智能逻辑 ---
        persistent lastRoi;
        if currentIdx == 1
            hFig = figure('Name', 'ROI', 'NumberTitle', 'off');
            imshow(frame1);
            title(sprintf('[%d/%d] 框选实验箱区域', currentIdx, totalFiles), 'Color', 'r', 'FontSize', 12);
            roiRect = drawrectangle('Color', 'r', 'Label', 'Arena');
            customWait(roiRect); 
            roiOut = roiRect.Position;
            lastRoi = roiOut;
            close(hFig);
            
            if totalFiles > 1
                choice = questdlg('后续视频是否沿用此方框?', '批处理', '是', '否', '是');
                if strcmp(choice, '否'), lastRoi = []; end
            end
        else
            if ~isempty(lastRoi)
                roiOut = lastRoi;
            else
                hFig = figure('Name', 'ROI', 'NumberTitle', 'off');
                imshow(frame1);
                roiRect = drawrectangle('Color', 'r');
                customWait(roiRect); 
                roiOut = roiRect.Position; close(hFig);
            end
        end

        % --- 步骤C: 追踪循环 (背景减除法) ---
        videoObj.CurrentTime = 0;
        rawTraj = nan(nFrames, 2);
        seDisk = strel('disk', sys.erodeRadius);
        
        % 预裁剪背景
        bgROI = imcrop(bgImage, roiOut);
        
        % GIF 准备
        gifFileName = fullfile(outDir, [fName '_Tracking.gif']);
        firstGifFrame = true;
        
        % 可视化
        hVis = figure('Name', ['Tracking: ' fName], 'NumberTitle', 'off', 'Visible', 'on');
        ax = axes('Parent', hVis);
        hImg = imshow(frame1, 'Parent', ax); hold on;
        rectangle('Position', roiOut, 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--');
        hPt = plot(NaN, NaN, 'g.', 'MarkerSize', 20);
        hLine = plot(NaN, NaN, 'r-', 'LineWidth', 1);
        
        wb = waitbar(0, 'Processing...');
        
        frameCount = 0;
        while hasFrame(videoObj)
            frameCount = frameCount + 1;
            frameRGB = readFrame(videoObj);
            
            % 1. 裁剪
            frameROI = imcrop(frameRGB, roiOut);
            if size(frameROI,3)==3, grayROI = rgb2gray(frameROI); else, grayROI=frameROI; end
            
            % 2. 核心算法：背景减除 (比形态学更保真)
            % 背景(白) - 当前帧(黑) = 差异(正值)
            diffImg = imsubtract(bgROI, grayROI);
            
            % 3. 二值化
            bw = imbinarize(diffImg, sys.diffThresh);
            
            % 4. 适度去噪 (只去细线，不伤老鼠)
            bwClean = imopen(bw, seDisk);
            bwClean = imfill(bwClean, 'holes');
            
            % 5. 提取质心
            stats = regionprops(bwClean, 'Centroid', 'Area');
            validIdx = find([stats.Area] > sys.minArea);
            
            if ~isempty(validIdx)
                [~, maxI] = max([stats(validIdx).Area]);
                c = stats(validIdx(maxI)).Centroid;
                gX = c(1)+roiOut(1); gY = c(2)+roiOut(2);
                rawTraj(frameCount, :) = [gX, gY];
            end
            
            % 6. GIF 与显示
            if mod(frameCount, vis.drawInterval) == 0 && isvalid(hVis)
                set(hImg, 'CData', frameRGB);
                if ~isnan(rawTraj(frameCount,1))
                    set(hPt, 'XData', rawTraj(frameCount,1), 'YData', rawTraj(frameCount,2));
                    startP = max(1, frameCount-50);
                    set(hLine, 'XData', rawTraj(startP:frameCount,1), 'YData', rawTraj(startP:frameCount,2));
                end
                drawnow limitrate;
                
                if expP.saveGIF && mod(frameCount, expP.gifFrameSkip) == 0
                    frame = getframe(hVis); 
                    im = frame2im(frame); 
                    [imind, cm] = rgb2ind(im, 256); 
                    if firstGifFrame
                        imwrite(imind, cm, gifFileName, 'gif', 'Loopcount', inf, 'DelayTime', 0.1); 
                        firstGifFrame = false; 
                    else 
                        imwrite(imind, cm, gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1); 
                    end 
                end
            end
            if mod(frameCount, 100) == 0, waitbar(frameCount/nFrames, wb); end
        end
        delete(wb); if isvalid(hVis), close(hVis); end
        
        % --- 数据处理 ---
        cleanTraj = rawTraj;
        % 1. 跳变剔除 (100px)
        diffXY = diff(cleanTraj);
        dist = [0; sqrt(sum(diffXY.^2, 2))];
        cleanTraj(dist > sys.maxJump, :) = NaN; 
        
        % 2. 插值
        cleanTraj(:,1) = fillmissing(cleanTraj(:,1), 'linear');
        cleanTraj(:,2) = fillmissing(cleanTraj(:,2), 'linear');
        
        % 3. 平滑
        finalTraj = zeros(size(cleanTraj));
        finalTraj(:,1) = smoothdata(cleanTraj(:,1), 'movmean', sys.smoothWin);
        finalTraj(:,2) = smoothdata(cleanTraj(:,2), 'movmean', sys.smoothWin);
        
        % --- 导出 1: Nature 风格 (红线白底) ---
        hF1 = figure('Visible', 'off', 'Color', 'w');
        rectangle('Position', roiOut, 'EdgeColor', 'k', 'LineWidth', 2); hold on;
        plot(finalTraj(:,1), finalTraj(:,2), 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5);
        axis ij; axis equal; axis off;
        xlim([roiOut(1)-50, roiOut(1)+roiOut(3)+50]);
        ylim([roiOut(2)-50, roiOut(2)+roiOut(4)+50]);
        
        print(hF1, fullfile(outDir, [fName '_NatureStyle.png']), '-dpng', '-r300');
        savefig(hF1, fullfile(outDir, [fName '_NatureStyle.fig'])); 
        close(hF1);
        
        % --- 导出 2: 时间渐变图 (Gradient) ---
        hF2 = figure('Visible', 'off', 'Color', 'w');
        rectangle('Position', roiOut, 'EdgeColor', 'k', 'LineWidth', 2); hold on;
        
        x = finalTraj(:,1); y = finalTraj(:,2);
        z = zeros(size(x)); col = 1:length(x);
        patch([x;NaN], [y;NaN], [z;NaN], [col(:);NaN], ...
              'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2);
        
        colormap(jet); c = colorbar; c.Label.String = 'Time (Frame)';
        axis ij; axis equal; axis off;
        xlim([roiOut(1)-50, roiOut(1)+roiOut(3)+50]);
        ylim([roiOut(2)-50, roiOut(2)+roiOut(4)+50]);
        title('Time-Coded Trajectory');
        
        print(hF2, fullfile(outDir, [fName '_Gradient.png']), '-dpng', '-r300');
        savefig(hF2, fullfile(outDir, [fName '_Gradient.fig']));
        close(hF2);
        
        % --- 导出 CSV ---
        T = table((1:nFrames)', finalTraj(:,1), finalTraj(:,2), ...
             'VariableNames', {'Frame', 'X', 'Y'});
        writetable(T, fullfile(outDir, [fName '_Data.csv']));
        
        fprintf('  -> 结果已保存: %s\n', outDir);
    end

    function customWait(hROI)
        l = addlistener(hROI,'ROIClicked',@(src, evt) uiresume);
        uiwait; delete(l);
    end
end