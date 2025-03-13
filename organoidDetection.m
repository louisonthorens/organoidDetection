% ------------------------------------------------------------------------------
% Automated Organoids Detection & Image Stitching
% Author: Louison Thorens - Institute of Mechanobiology - Northeastern University
% Contact: l.thorens@northeastern.edu
% Date: 2025-03-13
% Version: 1.0
%
% Description:
% This script automates image stitching and detection of organoids.
%
% License:
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/
%
% If you use this code, please cite:
% - Paper: TBD
% - GitHub: Louison Thorens. (2025). Automated Cell Detection on Plates (v1.0)
%           Available at: https://github.com/louisonthorens/organoidDetection
% ------------------------------------------------------------------------------

close all; clearvars; 

% Initialize a table to store results
results = table('Size', [0 4], ...
               'VariableTypes', {'string', 'double', 'double', 'double'}, ...
               'VariableNames', {'Label', 'Num_Cells', 'Average_Area_um2', 'Median_Area_um2'});

folder = 'dataExample';

%%%
% PARAMETERS
%%%
plateLetters = {'B'}; % for running through all plate positions
plateNumbers = 1:2;

radiusCirclemm = 6.5; % radius of the container in mm

overlap_percentagex = 0.03; % 3% overlap in x
overlap_percentagey = 0.02; % 2% overlap in y

rows = 3; cols = 3; % arrangement of images for stitching

contrast = [100 3000];

binaryThreshold = 0.4;

thresholdArea = 2000; % in px2

%%%
% MAIN LOOP
%%%

for iLetter = 1:length(plateLetters)
    letter = plateLetters{iLetter};
    for iAB = plateNumbers
        %%% 
        % LOADING AND PRE-PROCESSING OF IMAGE
        %%%
        fprintf('Processing %s%d...\n', letter, iAB);
        
        xdceFile = dir([folder '/*.xdce']);
        imageInfo = BioformatsImage([folder '/' xdceFile.name]);
        scale = imageInfo.pxSize(1);
       
        files = ['*' letter ' - ' num2str(iAB) '(fld * wv Green - dsRed).tif'];

        allFiles = dir(fullfile(folder, files));

        if isempty(allFiles)
            warning('No files found for %s. Skipping...', files);
            continue;
        end

        sample_img = imread(fullfile(folder, allFiles(1).name));
        [height, width] = size(sample_img);
        image_size = [height width];

        % Compute overlap in pixels
        overlap_x = round(width * overlap_percentagex);
        overlap_y = round(height * overlap_percentagey);

        % Calculate new dimensions after considering overlap
        step_width = width - overlap_x;
        step_height = height - overlap_y;

        % Initialize the large image canvas and weight matrices
        stitched_width = step_width * cols + overlap_x;
        stitched_height = step_height * rows + overlap_y;
        stitched_image = zeros(stitched_height, stitched_width);
        weight_image = zeros(stitched_height, stitched_width); % For averaging

        % Process and stitch the images
        % Straight-forward stitching using known overlap
        for r = 1:rows
            for c = 1:cols
                % Compute linear index
                idx = (r - 1) * cols + c;

                % Check if idx exceeds number of available files
                if idx > length(allFiles)
                    warning('Index %d exceeds number of available files. Skipping...', idx);
                    continue;
                end

                % Read and crop the current image
                img = imread(fullfile(folder, allFiles(idx).name));
                cropped_img = img(1:step_height + overlap_y, 1:step_width + overlap_x); % Include overlap
                % Compute position in the large image
                start_x = (c - 1) * step_width + 1;
                start_y = (r - 1) * step_height + 1;

                end_x = start_x + size(cropped_img, 2) - 1;
                end_y = start_y + size(cropped_img, 1) - 1;

                % Ensure the stitched image is large enough
                if end_x > stitched_width || end_y > stitched_height
                    error('Stitched image size is too small. Increase stitched_width or stitched_height.');
                end

                % Accumulate pixel values for averaging
                stitched_image(start_y:end_y, start_x:end_x) = stitched_image(start_y:end_y, start_x:end_x) + double(cropped_img);
                weight_image(start_y:end_y, start_x:end_x) = weight_image(start_y:end_y, start_x:end_x) + 1;
            end
        end

        % Compute the average in overlapping regions
        stitched_image = double(stitched_image) ./ weight_image;                     

        % Normalize the stitched image
        stitched_image = (stitched_image - min(contrast)) / (max(contrast) - min(contrast));
        stitched_image(stitched_image < 0) = 0; 
        stitched_image(stitched_image > 1) = 1;

        % Apply Non-Local Means Filtering - time consuming step
        % stitched_image = imnlmfilt(stitched_image);
        % disp('NLM done')
        
        % Binarization of the image using user threshold
        BG = stitched_image>0.01;
        stats = regionprops(BG); areas = cat(1,stats.Area); [~,idMax] = max(areas);
        center = stats(idMax).Centroid;
        mask = zeros(size(stitched_image)); [X,Y] = meshgrid(1:size(mask,2),1:size(mask,1));
        radiusCircle = 6.5*1e3/scale/2;
        mask(sqrt((X-center(1)).^2 + (Y-center(2)).^2) <= radiusCircle) = 1;
        imagesc(mask.*stitched_image)
        
        stitched_image = stitched_image.*mask;

        % Display the image (optional)
        figure();
        imagesc(stitched_image);
        clim([0 1]);
        drawnow

        % Create a red colormap
        red_colormap = [linspace(0, 1, 256)', zeros(256, 1), zeros(256, 1)];
        colormap(red_colormap);

        colorbar_handle = colorbar;
        ylabel(colorbar_handle, 'Normalized Fluorescence Intensity', 'FontSize', 12); % Add label to colorbar
        axis equal; axis off;

        title(files);
        
        % Save the stitched image
        clean_filename = strrep(files, '*', ''); % Remove '*' if present
        output_image_path = fullfile(folder, '/results', [clean_filename]);
        imwrite(uint16((stitched_image)*2^16), output_image_path);
        
        % Binarization of the image using user threshold
        bw = imbinarize(stitched_image, binaryThreshold);
        bw = imclose(bw, strel('disk', 5));
        bw = imfill(bw, 'holes');

        stats = regionprops(bw, 'Centroid', 'Area','PixelList');

        areas = [stats.Area];
        stats = stats(areas >= thresholdArea);

        if isempty(stats)
            warning('No regions found for %s%d with area >= %d. Skipping...', letter, iAB, thresholdArea);
            continue;
        end

        centers = cat(1, stats.Centroid);

        % Initialize a cell array to store boundaries
        num_regions = numel(stats); % Number of regions
        boundaries = cell(num_regions, 1);

        % Loop through each region to extract boundaries
        for k = 1:num_regions
            % Create a blank binary image for the current region
            region_image = false(stitched_height, stitched_width); % Updated to stitched image size

            % Get the PixelList for the current region
            pixel_list = stats(k).PixelList; % Nx2 array [x, y]

            % Convert PixelList to linear indices and mark in binary image
            linear_indices = sub2ind([stitched_height, stitched_width], round(pixel_list(:,2)), round(pixel_list(:,1)));
            region_image(linear_indices) = true;

            % Extract boundaries using bwboundaries
            boundary = bwboundaries(region_image, 'noholes'); % Ignore holes for clean boundary

            if ~isempty(boundary)
                boundaries{k} = boundary{1};
            else
                boundaries{k} = [];
            end
        end

        % Display the stitched image with boundaries
        figure()
        imagesc(stitched_image);
        stitched_imageColor = zeros(size(stitched_image,1),size(stitched_image,2),3);
        stitched_imageColor(:,:,1) = stitched_image;
        overlayBoundaries = zeros(size(stitched_imageColor));
        clim([0 1]);
        colormap(red_colormap);
        axis equal; axis off;
        hold on;
        for k = 1:num_regions
            boundary = boundaries{k};
            for k = 1:length(boundary)
                overlayBoundaries(boundary(k,1), boundary(k,2),:) = 1;
            end
            if ~isempty(boundary)
                plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5); % (x, y) -> column, row
            end
        end

        finalImageColor = stitched_imageColor + imdilate(overlayBoundaries,strel('disk',5));
        
        % Compute metrics
        num_cells = length(stats);
        average_area = mean([stats.Area]);
        median_area = median([stats.Area]);

        % Save the annotated image
        annotated_filename = strrep(clean_filename, '.tif', '.jpg'); % Change extension to .jpg
        annotated_image_path = fullfile(folder, '/results', annotated_filename);
        imwrite(finalImageColor,annotated_image_path)

        % Close the stitched image figure
        close(gcf);

        % Create label 
        label = sprintf('%s%d', letter, iAB);

        % Compute scaled areas
        average_area_um2 = average_area * scale^2;
        median_area_um2 = median_area * scale^2;

        % Append the data to the results table
        newRow = {label, num_cells, average_area_um2, median_area_um2};
        results = [results; newRow];
    end
end

%%%
% SAVING
%%%

% Define the output CSV file path
output_csv = fullfile(folder, '/results', 'summary_results.csv');

% Write the table to CSV
writetable(results, output_csv);

fprintf('Results have been saved to %s\n', output_csv);