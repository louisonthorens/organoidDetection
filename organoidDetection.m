% ------------------------------------------------------------------------------
% Automated Organoids Detection (Template Matching & Bottleneck Splitting)
% Author: Louison Thorens - Institute of Mechanobiology - Northeastern University
% Contact: l.thorens@northeastern.edu
% Date: 2026-02-23
% Version: 2.0
%
% Description:
% This script automates the detection and segmentation of organoids in stitched 
% grayscale images using 2D normalized cross-correlation (NCC) against a filter 
% bank of disk templates. It incorporates adaptive thresholding and morphological 
% bottleneck analysis to computationally bisect overlapping clusters.
%
% License:
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/
% ------------------------------------------------------------------------------

close all; clearvars; clc;

%%%
% PARAMETERS
%%%
folder             = '';
plateLetters       = {'B'}; 
plateNumbers       = 1;

% Template Filter Bank Parameters
radii              = 9:1:40;       % Radii range for templates (px)
rBoundary          = 4;            % Edge thickness for hollow templates (px)

% Image Filtering & Thresholding Parameters
blurSigma          = 5;            % Sigma for initial Gaussian blur
bgSmoothSigma      = 100;          % Sigma for local background estimation
localContrastThresh= 2.5;          % SNR threshold for initial binarization
BWthreshold        = 500;          % Image background intensity threshold

% Morphological Analysis Parameters
minAreaFilter      = 100;          % Minimum area to keep initial regions (px^2)
maxEccentricity    = 0.99;         % 0=circle, ->1=line (filters out artifacts)
eccentricityCut    = 0.7;          % Threshold to trigger bottleneck splitting
finalAreaThreshold = 200;          % Final minimum area for valid organoids (px^2)

%%%
% MAIN LOOP
%%%
for iLetter = 1:length(plateLetters)
    letter = plateLetters{iLetter};
    
    for number = plateNumbers
        close all;
        fprintf('Processing %s - %d...\n', letter, number);
        
        %%% 
        % LOADING AND PRE-PROCESSING OF IMAGE
        %%%
        imgFile = fullfile(folder, [letter ' - ' num2str(number) '_stitched_gray.tif']);
        if ~isfile(imgFile)
            warning('File %s not found. Skipping...', imgFile);
            continue;
        end
        img = double(imread(imgFile));
        
        %%% 
        % GENERATE TEMPLATE FILTER BANK
        %%%
        nR = numel(radii);
        templates = cell(1, nR * 2); % Allocate for solid and edge templates
        
        % 1. Solid disk templates
        for ii = 1:nR
            r = radii(ii);
            t = fspecial('disk', r);          % filled disk
            t = t ./ sum(t(:));               % unit energy
            templates{ii} = t;
        end
        
        % 2. Edge-emphasized (dim center) templates
        for ii = 1:nR
            r = radii(ii);
            t = double(fspecial('disk', r));  
            t(t > 0) = 1;
            [x, y] = meshgrid(1:size(t,1), 1:size(t,2));
            % Dim the interior slightly
            t(((x-mean(x(1,:))).^2 + (y-mean(y(:,1))).^2) < (r-rBoundary).^2) = 0.5;
            t = t ./ sum(t(:));               
            templates{ii+nR} = t;
        end
        
        %%% 
        % CROSS-CORRELATION (TEMPLATE MATCHING)
        %%%
        I = img;
        I = I - mean(I(:));
        I = I ./ sqrt(sum(I(:).^2) + eps);     % zero‑mean, unit‑norm
        Iblurred = imgaussfilt(I, blurSigma);
        
        Rmax   = -Inf(size(I));                % best score so far
        Rscale = zeros(size(I));               % radius that produced Rmax
        
        % Initialize command-line progress tracker
        fprintf('Cross-correlating %d templates\n', nR);
        
        % Loop through all templates and find best matches
        for ii = 1:nR

            fprintf('%d / %d \n',ii,nR);
            
            t = templates{ii};
            t = t - mean(t(:));                % zero‑mean, unit‑norm template
            t = t ./ sqrt(sum(t(:).^2) + eps);
            
            r_temp = floor(size(t,1)/2);
            R = normxcorr2(t, Iblurred);
            R = R(r_temp:end-r_temp-1, r_temp:end-r_temp-1); % Crop to original size
            
            % Update the “winner‑takes‑all” maps
            mask         = R > Rmax;
            Rmax(mask)   = R(mask);
            Rscale(mask) = radii(ii);        
        end
        fprintf('\n');
        
        %%% 
        % ADAPTIVE THRESHOLDING & INITIAL MASKS
        %%%
        newImage = Rmax .* img; % Amplify original image with correlation score
        newImage(newImage <= 1e-4) = 0;

        BG = img >= BWthreshold;
        
        % Local background correction
        A = newImage ./ imgaussfilt(newImage, bgSmoothSigma);
        A(A > 10) = 10; % Cap extreme outliers

        A(BG==0) = 0;
        
        % Binarization and morphological cleanup
        BW = A > localContrastThresh;
        BW = imfill(BW, 'holes');
        I  = newImage;
        
        % Extract region properties
        CC    = bwconncomp(BW);
        stats = regionprops(CC, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', ...
            'Orientation', 'Area', 'Eccentricity', 'PixelIdxList');
        
        % Apply initial physical filters
        keep = [stats.Area] >= minAreaFilter & ...
               [stats.Eccentricity] <= maxEccentricity & ...
               [stats.MinorAxisLength] > 0;
        stats = stats(keep);
        
        %%% 
        % MORPHOLOGICAL SPLITTING (BOTTLENECK ANALYSIS)
        %%%
        allBoundaries = {}; 
        areas = [];
        t_param = linspace(0, 2*pi, 200); % For parametric ellipses

        % Initialize command-line progress tracker
        fprintf('Analyzing %d objects\n', numel(stats));          
        
        for k = 1:numel(stats)

            fprintf('%d / %d \n',k,numel(stats));

            cx = stats(k).Centroid(1);
            cy = stats(k).Centroid(2);
            a  = stats(k).MajorAxisLength/2;   
            b  = stats(k).MinorAxisLength/2;   
            th = -deg2rad(stats(k).Orientation);
            e  = stats(k).Eccentricity;
            
            % Isolate the current component mask
            maskHere = I*0; 
            maskHere(stats(k).PixelIdxList) = 1;
            maskHere = imclose(maskHere, strel('disk', 5));
            
            % If highly eccentric, attempt to cut
            if e > eccentricityCut
                x1 = cx + a*cos(0)*cos(th) - b*sin(0)*sin(th); 
                y1 = cy + a*cos(0)*sin(th) + b*sin(0)*cos(th);
                x2 = cx + a*cos(pi)*cos(th) - b*sin(pi)*sin(th); 
                y2 = cy + a*cos(pi)*sin(th) + b*sin(pi)*cos(th);
                
                halfWidth = 50;
                [H, W] = size(I);
                xm = linspace(x1, x2, 1e3);                     
                ym = linspace(y1, y2, 1e3);
                dx = x2 - x1;                         
                dy = y2 - y1;
                len = hypot(dx, dy);
                
                if len < 1e-6
                    dx = cos(th);                    
                    dy = -sin(th);
                    len = hypot(dx, dy);
                end
                
                nx = -dy/len; 
                ny = dx/len;            
                Npts = max(3, 2*round(halfWidth)+1);
                s    = linspace(-halfWidth, halfWidth, Npts);
                xN   = xm' + s*nx;
                yN   = ym' + s*ny;
                
                xN = min(max(xN, 1), W);
                yN = min(max(yN, 1), H);
                
                % Sample profile along the normal
                profile = interp2(medfilt2(I, [5 5]), xN, yN);
                profile = mean(profile, 2);
                profile = profile / max(profile);
                
                firstPos = find(profile >= 0.35, 1, 'first'); 
                lastPos  = find(profile >= 0.35, 1, 'last');
                profile(1:firstPos) = 1; 
                profile(lastPos:end) = 1;
                
                toRemove = find(profile < 0.35); 
                toRemove = toRemove(1:min(5, length(toRemove)));
                toRemovej = round(xN(toRemove, :)); 
                toRemovei = round(yN(toRemove, :));
                
                toRemoveIndex = sub2ind(size(maskHere), toRemovei, toRemovej);
                newMask = maskHere;
                newMask(toRemoveIndex) = 0;
                newMask = bwareaopen(newMask, 50);
                
                statsHere = regionprops(bwmorph(newMask, 'thicken', 10), 'Centroid');
                newMask = newMask*0;
                centers = cat(1, statsHere.Centroid);
                for ic = 1:size(centers, 1)
                    newMask(round(centers(ic, 2)), round(centers(ic, 1))) = 1;
                end
                
                newMask = bwmorph(newMask, 'thicken', 1e2);
                boundaries = bwboundaries(newMask .* maskHere);
                
                % Check boundaries for neck points
                for iBoundary = 1:length(boundaries)
                    boundary = boundaries{iBoundary};
                    xB = boundary(:, 2); 
                    yB = boundary(:, 1);
                    
                    [doCut, i1, i2, width, score] = neck_decide_from_boundary(xB, yB);
                    
                    if ~doCut
                        allBoundaries = cat(1, allBoundaries, boundary);
                        areas(end+1) = polyarea(xB, yB);
                    else
                        % Bisect the boundary
                        sorted = sort([i1, i2]);
                        i1 = sorted(1); i2 = sorted(2);
                        
                        xB1 = xB(i1:i2); yB1 = yB(i1:i2);
                        xB2 = xB; xB2(i1:i2) = []; 
                        yB2 = yB; yB2(i1:i2) = [];
                        
                        allBoundaries = cat(1, allBoundaries, {[yB1, xB1]});
                        areas(end+1) = polyarea(xB1, yB1);
                        allBoundaries = cat(1, allBoundaries, {[yB2, xB2]});
                        areas(end+1) = polyarea(xB2, yB2);
                    end
                end
            else
                % Normal cell, keep boundary as is
                boundaries = bwboundaries(maskHere);
                allBoundaries = cat(1, allBoundaries, boundaries);
                for iBoundary = 1:length(boundaries)
                    boundary = boundaries{iBoundary};
                    xB = boundary(:, 2); 
                    yB = boundary(:, 1);
                    areas(end+1) = polyarea(xB, yB);
                end
            end
        end
        fprintf('\n');
        
        %%%
        % SAVING OUTPUTS
        %%%
        close all;
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        imagesc(img);
        hold on; axis equal; colormap gray; axis off;
        clim([0 5e3]);
        
        % Filter out final small noise objects
        goodOnes = find(areas >= finalAreaThreshold);
        for k = goodOnes
            plot(allBoundaries{k}(:, 2), allBoundaries{k}(:, 1), 'LineWidth', 2);
        end
        title(sprintf('Detected %d organoids', length(goodOnes)));
        
        % Save MAT and TIF
        matName = fullfile(folder, [letter ' - ' num2str(number) '_stitched_gray.mat']);
        tifName = fullfile(folder, [letter ' - ' num2str(number) '_stitched_gray_detection.tif']);
        
        save(matName, 'allBoundaries', 'areas');
        saveas(gcf, tifName);
        
        fprintf('Saved results for %s - %d\n', letter, number);
    end
end

%%%
% HELPER FUNCTIONS
%%%
function [doCut, i1, i2, width, score] = neck_decide_from_boundary(x, y, params)
    % Decide if a component should be split at a "neck" using only its boundary.
    if nargin < 3
        params.arcFracMin  = 0.15;
        params.maxWidthPx  = 75;
        params.scoreThr    = 0.3;
        params.sgolWinFrac = 0.02;
    end
    
    % Ensure open polygon
    if x(1)==x(end) && y(1)==y(end)
        x(end)=[]; y(end)=[];
    end
    x = double(x(:)); y = double(y(:));
    N = numel(x);
    
    if N < 10
        doCut=false; i1=[]; i2=[]; width=Inf; score=Inf; return; 
    end
    
    % Smooth boundary
    win = max(5, 2*floor(params.sgolWinFrac*N)+1); 
    if mod(win, 2) == 0, win = win + 1; end
    
    try
        x = sgolayfilt(x, 3, win); 
        y = sgolayfilt(y, 3, win);
    catch
        x = movmean([x(end-win+1:end); x; x(1:win)], win); x = x(win+1:win+N);
        y = movmean([y(end-win+1:end); y; y(1:win)], win); y = y(win+1:win+N);
    end
    
    % Perimeter via arc lengths
    perim = sum(hypot(diff([x; x(1)]), diff([y; y(1)])));
    
    idx  = (1:N)';
    dIdx = abs(idx - idx');                        
    dIdx = min(dIdx, N - dIdx);                    
    arcDist = dIdx * (perim / N);
    
    % Pairwise Euclidean distance
    Dx = x - x'; Dy = y - y';
    D  = hypot(Dx, Dy);
    
    % Candidate mask
    amin  = params.arcFracMin * perim;
    cand  = arcDist >= amin;
    D(~cand) = Inf;         
    D(D < 2) = Inf;         
    
    % Score and Decision
    Score = D ./ (arcDist + eps);
    [score, idxLin] = min(Score(:));
    
    if ~isfinite(score)
        doCut=false; i1=[]; i2=[]; width=Inf; return;
    end
    
    [i1, i2] = ind2sub([N N], idxLin);
    width   = D(i1, i2);
    doCut   = (width <= params.maxWidthPx) && (score <= params.scoreThr);
end