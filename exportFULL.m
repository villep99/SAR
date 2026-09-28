fprintf('\nPreparing export..\n');

range_step = 1;
azimuth_step = 5;

sigma0_export = sigma0_ground_full(1:range_step:end, 1:azimuth_step:end);
sigma0_export_dB = 10*log10(max(sigma0_export, eps));

fprintf('Export dimensions: %d x %d pixels\n', size(sigma0_export,1), size(sigma0_export,2));
fprintf('Range pixel spacing: %.6f m\n', target_ground_spacing * range_step);
fprintf('Azimuth pixel spacing: %.6f m\n', azimuth_spacing * azimuth_step);

%% Square pixels

fprintf('\nResampling to square pixels..\n');

range_spacing_current = target_ground_spacing * range_step;
azimuth_spacing_current = azimuth_spacing * azimuth_step;

target_spacing = azimuth_spacing_current;

range_current = (0:size(sigma0_export,1)-1) * range_spacing_current;
range_target = 0:target_spacing:range_current(end);

sigma0_square = zeros(length(range_target), size(sigma0_export,2), 'single');

for col = 1:size(sigma0_export,2)

    sigma0_square(:,col) = interp1(range_current, sigma0_export(:,col), range_target, 'linear', 0);

end

fprintf('Original export size: %d x %d pixels\n', size(sigma0_export,1), size(sigma0_export,2));
fprintf('Square-pixel size: %d x %d pixels\n', size(sigma0_square,1), size(sigma0_square,2));
fprintf('Square pixel spacing: %.6f m\n', target_spacing);

clear azimuth_spacing azimuth_spacing_current azimuth_step col range_current range_spacing_current range_step range_target sigma0_export sigma0_export_dB sigma0_ground_full target_ground_spacing target_spacing

%% Export

fprintf('\nPreparing PNG export..\n');

sigma0_square_dB = 10*log10(max(sigma0_square, eps));

display_min = double(clim_low);
display_max = double(clim_high);

png_image = uint8(255 * mat2gray(sigma0_square_dB, [display_min display_max]));
imwrite(png_image, 'ICEYE_GRD_image.png');

fprintf('PNG export completed\n');
fprintf('Image size: %d x %d pixels\n', size(png_image,1), size(png_image,2));
fprintf('Dynamic range: %.1f to %.1f dB\n', display_min, display_max);