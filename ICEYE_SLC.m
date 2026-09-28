clear;
clc;

%% Load ICEYE SLC data and metadata

filename = 'ICEYE_STT1R7_20251109T094946Z_6983691_X49_SLED_SLC.tif';
json_filename = 'ICEYE_STT1R7_20251109T094946Z_6983691_X49_SLED_SLC.json';

metadata = jsondecode(fileread(json_filename));

range_spacing   = metadata.properties.("sar_pixel_spacing_range");
azimuth_spacing = metadata.properties.("sar_pixel_spacing_azimuth");

PRF = metadata.properties.("iceye_acquisition_prf");

incidence_near = metadata.properties.("iceye_incidence_angle_near");
incidence_far = metadata.properties.("iceye_incidence_angle_far");

calibration_factor = metadata.properties.("iceye_calibration_factor");

fprintf('ICEYE SLC metadata loaded\n');
fprintf('Range spacing: %.9f m\n', range_spacing);
fprintf('Azimuth spacing: %.9f m\n', azimuth_spacing);
fprintf('Acquisition PRF: %.3f Hz\n', PRF);
fprintf('Near incidence angle: %.6f deg\n', incidence_near);
fprintf('Far incidence angle: %.6f deg\n', incidence_far);
fprintf('Calibration factor: %.15g\n', calibration_factor);

%% SLC TIFF data

t = Tiff(filename,'r');
t.setDirectory(1);

image_height = t.getTag('ImageLength');
image_width = t.getTag('ImageWidth');

tile_size = t.getTag('TileLength');

num_tiles_x = ceil(image_width / tile_size);
num_tiles_y = ceil(image_height / tile_size);

total_tiles = num_tiles_x * num_tiles_y;

fprintf('\nSLC TIFF information:\n');
fprintf('Image dimensions: %d x %d pixels\n', image_height, image_width);
fprintf('Tile size: %d x %d pixels\n', tile_size, tile_size);
fprintf('Tiles in range direction: %d\n', num_tiles_x);
fprintf('Tiles in azimuth direction: %d\n', num_tiles_y);
fprintf('Total tiles: %d\n', total_tiles);

%% Allocate for full complex SLC

SLC_full = complex(zeros(image_height, image_width, 'single'));
fprintf('\nAllocated SLC_full: %.2f GB\n',numel(SLC_full)*8/1e9);

%% Load complete SLC tile-by-tile

fprintf('\nLoading complete SLC tile-by-tile..\n');

for tile_row = 1:num_tiles_y

    for tile_col = 1:num_tiles_x

        tile_number = (tile_row - 1) * num_tiles_x + tile_col;
        tile = t.readEncodedTile(tile_number);
        row_start = (tile_row - 1) * tile_size + 1;
        col_start = (tile_col - 1) * tile_size + 1;
        row_end = min(row_start + tile_size - 1, image_height);
        col_end = min(col_start + tile_size - 1, image_width);

        nr = row_end - row_start + 1;
        nc = col_end - col_start + 1;

        complex_tile = tile(:,:,1) .* exp(1i * tile(:,:,2));

        SLC_full(row_start:row_end, col_start:col_end) = complex_tile(1:nr,1:nc);

    end

    fprintf('Loaded tile row %d / %d (%.1f %%)\n', tile_row, num_tiles_y, 100*tile_row / num_tiles_y);

end

fprintf('\nFull SLC loaded\n');
whos SLC_full

fprintf('SLC size: %d x %d\n', size(SLC_full,1), size(SLC_full, 2));
fprintf('Amplitude range: %.6f to %.6f\n', min(abs(SLC_full),[],'all'), max(abs(SLC_full),[],'all'));
fprintf('Phase range: %.6f to %.6f rad\n', min(angle(SLC_full),[],'all'), max(angle(SLC_full),[],'all'));

t.close();

%% Creating physical axes and incidence-angle profile

range_axis_full = (0:image_width-1) * range_spacing;
azimuth_axis_full = (0:image_width-1) * azimuth_spacing;
incidence_axis_full = linspace(incidence_near, incidence_far, image_height);

fprintf('\nFull-scene geometry:\n');
fprintf('Slant-range extent: %.2f m\n', range_axis_full(end));
fprintf('Azimuth extent: %.2f m\n', azimuth_axis_full(end));
fprintf('Near incidence angle: %.6f deg\n', incidence_axis_full(1));
fprintf('Far incidence angle: %.6f deg\n', incidence_axis_full(end));

%% Set up for image processing

block_rows = 1024;
num_blocks = ceil(image_height / block_rows);

fprintf('\nProcessing setup:\n');
fprintf('Block size: %d x %d pixels\n', block_rows, image_width);
fprintf('Number of blocks: %d\n', num_blocks);

%% Set up Five-look doppler processing

num_looks = 5;
num_azimuth = image_width;

doppler_axis = (-floor(num_azimuth/2):ceil(num_azimuth/2)-1) * (PRF / num_azimuth);
look_edges = linspace(doppler_axis(1), doppler_axis(end), num_looks + 1);

fprintf('\nFive-look processing setup:\n');
fprintf('Number of look: %d\n', num_looks);
fprintf('Doppler range: %.2f to %.2f Hz\n', doppler_axis(1), doppler_axis(end));

%% Allocate for full ground-range GRD image

ground_spacing_full = range_spacing ./ sind(incidence_axis_full);
ground_range_axis_full = zeros(size(ground_spacing_full));

for k = 2:length(ground_spacing_full)

    ground_range_axis_full(k) = ground_range_axis_full(k-1) + ground_spacing_full(k);

end

target_ground_spacing = 0.25;
full_ground_extent = ground_range_axis_full(end);
num_ground_pixels = floor(full_ground_extent / target_ground_spacing) + 1;
num_azimuth_pixels = image_width;

sigma0_ground_full = zeros(num_ground_pixels, num_azimuth_pixels, 'single');

fprintf('\nGRD setup:\n');
fprintf('Ground-range extent: %.2f m\n', full_ground_extent);
fprintf('Ground spacing: %.6f to %.6f m\n', ground_spacing_full(1), ground_spacing_full(end));
fprintf('GRD size: %d x %d pixels\n', size(sigma0_ground_full,1), size(sigma0_ground_full,2));
fprintf('GRD memory: %.2f GB\n', numel(sigma0_ground_full) * 4 / 1e9);

%% Process every range block

for block = 1:num_blocks

    row_start = (block - 1) * block_rows + 1;
    row_end = min(row_start + block_rows - 1, image_height);

    fprintf('\nProcessing block %d/%d: rows %d-%d\n', block, num_blocks, row_start, row_end);

    SLC_block = SLC_full(row_start:row_end, :);
    SLC_fft_az = fftshift(fft(SLC_block, [], 2), 2);
    
    fprintf('Azimuth FFT completed\n');

    look_power = zeros(size(SLC_block,1), size(SLC_block,2),'single');

    for look = 1:num_looks

        f_low = look_edges(look);
        f_high = look_edges(look + 1);

        mask = doppler_axis >= f_low & doppler_axis < f_high;

        SLC_look_fft = zeros(size(SLC_fft_az), 'like', SLC_fft_az);
        SLC_look_fft(:,mask) = SLC_fft_az(:,mask);
        SLC_look = ifft(ifftshift(SLC_look_fft,2),[],2);
        look_power = look_power + abs(SLC_look).^2;

        clear SLC_look_fft SLC_look

    end

    multilook_power = look_power / num_looks;

    fprintf('Five-look processing completed\n');

    block_incidence = incidence_axis_full(row_start:row_end);
    incidence_factor = sind(block_incidence(:));
    sigma0_block = calibration_factor .* multilook_power .* incidence_factor;

    fprintf('Radiometric calibration completed\n');

    block_ground_axis = ground_range_axis_full(row_start:row_end);
    block_ground_start = block_ground_axis(1);
    block_ground_end = block_ground_axis(end);

    ground_range_out = block_ground_start:target_ground_spacing:block_ground_end;

    sigma0_ground_block = zeros(length(ground_range_out), size(sigma0_block,2), 'single');

    for col = 1:size(sigma0_block,2)

        sigma0_ground_block(:,col) = interp1(block_ground_axis, sigma0_block(:,col), ground_range_out, 'linear', 0);

    end

    output_row_start = floor(block_ground_start / target_ground_spacing) + 1;
    output_row_end = min(output_row_start + size(sigma0_ground_block,1) -1, num_ground_pixels);
    num_output_rows = output_row_end - output_row_start + 1;

    sigma0_ground_full(output_row_start:output_row_end, :) = sigma0_ground_block(1:num_output_rows, :);

    fprintf('Ground-range resampling completed\n');
    fprintf('Stored GRD rows %d-%d\n', output_row_start, output_row_end);

    clear SLC_block SLC_fft_az look_power multilook_power sigma0_block sigma0_ground_block ground_range_out

end

clim_low_linear = prctile(sigma0_ground_full(:), 1);
clim_high_linear = prctile(sigma0_ground_full(:), 99);
clim_low = 10 * log10(max(clim_low_linear, eps));
clim_high = 10 * log10(max(clim_high_linear, eps));

fprintf('\n All range blocks processed\n');

%% Saving processed GRD

save('ICEYE_GRD_final.mat','clim_low', 'clim_high', 'sigma0_ground_full', 'azimuth_spacing', 'target_ground_spacing', '-v7.3');
fprintf('\nFinal GRD variables saved\n');