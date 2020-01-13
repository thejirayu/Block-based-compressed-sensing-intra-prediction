clear;
close all;
clc;
profile on
UpFolder = fileparts(pwd);
addpath(fullfile('l1magic/Optimization'));
addpath(fullfile('l1magic/Measurements'));
addpath(fullfile('l1magic/Data'));
addpath(genpath('Ride'));
opts = 0;
imageOriginalPath = 'Ride';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles);
%___SIMULATION SETUPS___
sub_pixels       = 16;
n                = sub_pixels*sub_pixels; % NOTE: small error still present after increasing m;
bpp_buffer       = 0;

measurement_matrix_lists        = [128];
measurement_matrix_construction = 'binary_walsh_hadamard';
image_reconstruction_algorithm  = 'l1_eq_pd';
image_transformation_algorithm  = 'ifwht';
color_mode                      = 'gray';
quartus_interface               = 'off';

stack_size = 3;

for matrix_loop = 1:length(measurement_matrix_lists)
    reset          = 1;
    switch measurement_matrix_lists(matrix_loop)
        case {256}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 1;
        case {192}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.75;
        case {128}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.50;
        case {64}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.25;
        case {32}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.125;
        case {16}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.0625;
        case {8}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.03125;
        case {4}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.015625;
        case {2}
            m = measurement_matrix_lists(matrix_loop);
            sampling_rate = 0.0078125;
    end

    switch measurement_matrix_construction
        case 'binary_walsh_hadamard'
            hadamard_matrix          = hadamard(n);
            HadIdx                   = 0:n-1;                         % Hadamard index
            M                        = log2(n)+1;                     % Number of bits to represent the index
            binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0'; % Bit reversing of the binary index
            binSeqIdx                = zeros(n,M-1);                  % Pre-allocate memory
            for k = M:-1:2
                % Binary sequency index
                binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
            end
            SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');    % Binary to integer sequency index
            walshMatrix              = hadamard_matrix(SeqIdx+1,:); % 1-based indexing
            phi                      = max(walshMatrix(1:m,1:n), 0);
    end

%___THETA___
%___NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues___
theta = zeros(m,n);
for theta_loop = 1:n
    ek = zeros(1,n);
    ek(theta_loop) = 1;
    switch image_transformation_algorithm
        case 'ifwht'
            psi = ifwht(ek)';
    end
    theta(:,theta_loop) = phi*psi;
end
    frame_index = 1;
    for frame_number = 1:60
        %___LOAD IMAGE___
        load_frame = imread(imageFiles(frame_number).name);
        %load_frame = snapshot(cam);
        if(strcmp(color_mode,'rgb') || strcmp(color_mode,'RGB'))
            frame = load_frame(:,:,:);
            frame = imresize(frame,[720, 1280]);
            plane = 3;
        else
            frame = rgb2gray(load_frame);
            frame = imresize(frame,[720, 1280]);
            plane = 1;
        end

        %___RESET STATE___
        N_1 = zeros(1, size(frame,1)/sub_pixels) + sub_pixels;
        N_2 = zeros(1, size(frame,2)/sub_pixels) + sub_pixels;

        for i = 1:plane
            C(:,:,i) = mat2cell(double(frame(:,:,i)), N_1, N_2);
        end

        %___RESET STATE___
        if(reset == 1)
            reset = 0;
            bits_shift                  = zeros(size(N_1,2), size(N_2,2), plane)+4;
            threshold                   = zeros(size(N_1,2), size(N_2,2), plane);
            y                           = zeros(m,1);
            y_buffer_up_encoder         = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_left_encoder       = zeros(m, 1);
            y_buffer_dc_encoder         = zeros(m, 1);
            y_buffer_cp_encoder         = (zeros(m, 1));
            y_buffer_up_decoder         = zeros((m), size(frame,2)/sub_pixels);
            y_buffer_left_decoder       = zeros(m, 1);
            y_buffer_dc_decoder         = zeros(m, 1);
            y_buffer_cp_decoder         = (zeros(m, 1));
            y_predicted_encoder         = zeros(m, 1);
            plane_and                   = zeros(1,3);
            image_bpp                   = zeros(1,3);
            image_psnr                  = zeros(1,3);
            image_ssim                  = zeros(1,3);
            y_buffer_cp_encoder(1)      = 65280/2;
            y_buffer_cp_encoder(2:m)    = 65280/4;
            y_buffer_cp_decoder(1)      = 65280/2;
            y_buffer_cp_decoder(2:m)    = 65280/4;
            if(strcmp(color_mode,'rgb') || strcmp(color_mode,'RGB'))
                bpp_buffer              = zeros(1,3);
            else
                bpp_buffer              = zeros(1,1);
            end
            temp_padding                = zeros(size(frame,1)+2, size(frame,2)+2, plane);
            res_temp_padding            = zeros(size(frame,1)+2, size(frame,2)+2, plane);
            buffer                      = cell(size(N_1,2), size(N_2,2), plane);
            buffer_1                    = cell(size(N_1,2), size(N_2,2), plane);
            buffer_2                    = cell(size(N_1,2), size(N_2,2), plane);
            estimate_background_encoder = cell(size(N_1,2), size(N_2,2), plane);
            estimate_background_decoder = cell(size(N_1,2), size(N_2,2), plane);
            y_residual                  = cell(size(N_1,2), size(N_2,2), plane);
            combination_frame           = cell(size(N_1,2), size(N_2,2), plane);
            y_quantized                 = cell(size(N_1,2), size(N_2,2), plane);
            y_dequantized               = cell(size(N_1,2), size(N_2,2), plane);
            reconstructed_image         = cell(size(N_1,2), size(N_2,2), plane);
            res_reconstructed_image     = cell(size(N_1,2), size(N_2,2), plane);
            modes                       = zeros(size(frame,1)/sub_pixels, size(frame,2)/sub_pixels, plane);
            res_video_buffer            = zeros(size(frame,1), size(frame,2), plane, 100);
            video_buffer                = zeros(size(frame,1), size(frame,2), plane, 100);
            quartus_output              = zeros(1200, 1, 128);
            for i = 1:size(N_1, 2)
                for j = 1:size(N_2, 2)
                    for k = 1:plane
                        buffer{i,j,k}                      = zeros(m, 1);
                        buffer_1{i,j,k}                    = zeros(m, 1);
                        buffer_2{i,j,k}                    = zeros(m, 1);
                        estimate_background_encoder{i,j,k} = zeros(m, 1);
                        estimate_background_decoder{i,j,k} = zeros(m, 1);
                        combination_frame{i,j,k}           = zeros(m, 1);
                        y_residual{i,j,k}                  = zeros(m, 1);
                        y_quantized{i,j,k}                 = zeros(m, 1);
                        y_dequantized{i,j,k}               = zeros(m, 1);
                        reconstructed_image{i,j,k}         = zeros(m, 1);
                        res_reconstructed_image{i,j,k}     = zeros(m, 1);
                    end
                end
            end
        end

        %___THE RANDOM PROJECTION___
        for k = 1:plane
            for i = 1:1:size(frame,1)/sub_pixels
                for j = 1:1:size(frame,2)/sub_pixels
                   one_block_image(:,:,k) = reshape(C{i,j,k}.',1,[])';
                   y = BCS_encoder(one_block_image(:,:,k), phi);
                   [y_buffer_left_encoder, y_buffer_up_encoder(:,j), y_buffer_dc_encoder, y_buffer_cp_encoder, y_predicted_encoder, modes(i,j,k)] = coding_method(y, ...
                                                                                                                                                                  phi, ...
                                                                                                                                                                  i, ...
                                                                                                                                                                  j, ...
                                                                                                                                                                  sub_pixels, ...
                                                                                                                                                                  m, ...
                                                                                                                                                                  n, ...
                                                                                                                                                                  y_buffer_up_encoder(:,j), ...
                                                                                                                                                                  y_buffer_left_encoder, ...
                                                                                                                                                                  y_buffer_dc_encoder, ...
                                                                                                                                                                  y_buffer_cp_encoder, ...
                                                                                                                                                                  'intra_prediction');
                    y_residual{i,j,k} = y-y_predicted_encoder;
                    if(frame_number == 1)
                        M = 16;
                    else
                        M = 16;
                    end
                    bits_shift(i,j,k)  = 2^log2(M);
                    y_quantized{i,j,k} = floor(y_residual{i,j,k}./bits_shift(i,j,k));
                    bpp_buffer(k)      = bpp_buffer(k) + measurement_entropy(y_quantized{i,j,k}, 1280*720);
                end
            end
        end
        
                                            %%%%%%%%%%%%%%%%%%%%
                                            %%%%% CHANNALS %%%%%
                                            %%%%%%%%%%%%%%%%%%%%
                                            
        for k = 1:plane
            internal_slot=1;
            for i = 1:size(frame,1)/sub_pixels
                for j = 1:size(frame,2)/sub_pixels
                   y_dequantized{i,j,k} = (y_quantized{i,j,k}*bits_shift(i,j,k));
                   if(modes(i,j,k) == 0)
                      combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_left_decoder;
                   elseif(modes(i,j,k) == 1)
                      combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_up_decoder(:,j);
                   elseif(modes(i,j,k) == 2)
                      combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_dc_decoder;
                   elseif(modes(i,j,k) == 3)
                      combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                   else
                      combination_frame{i,j,k} = y_dequantized{i,j,k} + y_buffer_cp_decoder;
                   end
                  [y_buffer_left_decoder, y_buffer_up_decoder(:,j), y_buffer_dc_decoder, y_buffer_cp_decoder] = coding_method(combination_frame{i,j,k}, ...
                                                                                                                              phi, ...
                                                                                                                              i, ...
                                                                                                                              j, ...
                                                                                                                              sub_pixels, ...
                                                                                                                              m, ...
                                                                                                                              n, ...
                                                                                                                              y_buffer_up_decoder(:,j), ...
                                                                                                                              y_buffer_left_decoder, ...
                                                                                                                              y_buffer_dc_decoder, ...
                                                                                                                              y_buffer_cp_decoder, ...
                                                                                                                              'intra_prediction');

                  reconstructed_image{i,j,k}         = BCS_reconstruction(combination_frame{i,j,k}, ...
                                                                          theta, ...
                                                                          phi, ...
                                                                          image_transformation_algorithm, ...
                                                                          image_reconstruction_algorithm, ...
                                                                          sub_pixels, ...
                                                                          opts);
                end
            end
        end
        %___FINAL PROCESS ZONE___%
        for k = 1:plane
            if(plane == 3)
               temp_padding(:,:,k) = padarray(floor(cell2mat(reconstructed_image(:,:,k))),[1 1],'symmetric','both');
            else
               temp_padding = padarray(floor(cell2mat(reconstructed_image(:,:))),[1 1],'symmetric','both'); %padding
            end

            %___Overlapped Filtering___ Y X
            for i = 2:size(temp_padding,1)-1
                for j = 2:size(temp_padding,2)-1
                    video_buffer(i-1,j-1,k,frame_number) = floor((temp_padding(i,j,k)+temp_padding(i,j-1,k)+temp_padding(i,j+1,k))/3);
                end
            end
        end

        %___RESET FOR NEW FRAME___
        y_buffer_up_encoder      = zeros((m), size(frame,2)/sub_pixels);
        y_buffer_left_encoder    = zeros(m, 1);
        y_buffer_dc_encoder      = zeros(m, 1);
        y_buffer_cp_encoder      = (zeros(m, 1));
        y_buffer_up_decoder      = zeros((m), size(frame,2)/sub_pixels);
        y_buffer_left_decoder    = zeros(m, 1);
        y_buffer_dc_decoder      = zeros(m, 1);
        y_buffer_cp_decoder      = (zeros(m, 1));
        y_buffer_cp_encoder(1)   = 65280/2;
        y_buffer_cp_encoder(2:m) = 65280/4;
        y_buffer_cp_decoder(1)   = 65280/2;
        y_buffer_cp_decoder(2:m) = 65280/4;
        modes                    = zeros(size(frame,1)/sub_pixels, size(frame,2)/sub_pixels, plane);

        %___QUATITATIVE MATRICES___
        %block_per_frame(frame_number) = block_counting;
        image_bpp(frame_number)  = sum(bpp_buffer);
        image_psnr(frame_number) = PSNR(uint8(video_buffer(:,:,:,frame_number)), frame);
        image_ssim(frame_number) = ssim(uint8(video_buffer(:,:,:,frame_number)), frame);
        bpp_buffer               = zeros(k,1);
    end

    save(fullfile(strcat('PSNR_Ride_', num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'image_psnr');
    save(fullfile(strcat('SSIM_Ride_', num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'image_ssim');
    save(fullfile(strcat('BPP_Ride_' , num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'image_bpp');

    video_out = VideoWriter(fullfile(strcat('Ride', ...
                                            '_Frame_Skip_', num2str(sampling_rate), ...
                                            '_', measurement_matrix_construction, ...
                                            '_', image_reconstruction_algorithm, ...
                                            '_', image_transformation_algorithm, ...
                                            '_', color_mode, ...
                                            '_BITS_SHIFT_', num2str(bits_shift(bits_shift_loop)), ...
                                            '_Linear_Filter', ...
                                            '.avi')), 'Uncompressed AVI'); %create the video object
    video_out.FrameRate = 30;                                    
    open(video_out); %open the file for writing
    for loop = 1:frame_number-1
       writeVideo(video_out, uint8(video_buffer(:,:,:,loop)));
    end
    close(video_out);
end
profile report
profile off