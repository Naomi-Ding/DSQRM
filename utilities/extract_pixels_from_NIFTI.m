%% process the input: NIFTI images

function [all_tumor, all_tumor_coord, sub_tumor_ratios, ...
    all_images, sub_tumor_coord, rm_idx] = ...
    extract_pixels_from_NIFTI(sample_name, image_folder)
% sample_name: a cell array (n,1), containing the sample name

n = length(sample_name);
% n = length(imagefiles); % sample size
all_images = cell(1);
all_tumor = cell(1);
all_tumor_coord = cell(1);
sub_tumor_ratios = zeros(1); % ratio of tumor subtypes
sub_tumor_coord = cell(1); % coordinates of tumor subtypes
idx = 0;
rm_idx = []; 
for i = 1:n
    image_name = dir(fullfile(image_folder, [sample_name{i}, '*']));
    if isempty(image_name)
        rm_idx = [rm_idx, i];
    else
        if length(image_name) == 1
            warning('The image & the mask are not match for the sample %s', sample_name{i});
            rm_idx = [rm_idx, i];
        elseif length(image_name) == 2
            idx = idx + 1;
            image = niftiread(fullfile(image_folder, image_name(1).name));
            mask = niftiread(fullfile(image_folder, image_name(2).name));
            
            if length(size(image)) == 3
                T = mask ~= 0;
                [~, loc] = max(sum(sum(T, 2), 1));
                mask_i = mask(:,:,loc);
                image_i = image(:,:,loc);
            elseif length(size(image)) == 2
                image_i = image;
                mask_i = mask;
            end
            all_images{idx} = image_i;
%             all_masks{idx} = mask_i;
            [z1, z2] = ind2sub(size(image_i), 1:numel(image_i));
            coord_i = [z1; z2]';
            %% extract the pixel intensities for the tumor region
            T_i = mask_i ~= 0;
            T_i = T_i(:);
            coord_tumor = coord_i(T_i, :);
            tumor = image_i(:);
            tumor = tumor(T_i);
            all_tumor{idx} = tumor;
            all_tumor_coord{idx} = coord_tumor;
            subtypes = sort(unique(mask_i));
            ntypes = length(subtypes) - 1;
            sub_tumor_coord{idx} = cell(ntypes,1);
            ct = 0;
            for j = subtypes'
                if j ~= 0
                    ct  = ct + 1;
                    sub_tumor_ratios(idx, ct) = sum(mask_i == j, 'all') / sum(T_i);
                    % coordinates for each subtype of tumors
                    T_j = mask_i == j;
                    T_j = T_j(:);
                    sub_tumor_coord{idx}{ct} = coord_i(T_j, :);
                end
            end
        end
    end
end

end


