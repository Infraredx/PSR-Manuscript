function BkgndCorrect_histology(image_dir)
%%
% Function will search for TIF files of blank slide (filename has "_blank"). 
% Use this image to correct histology images in the same folder. 
% Bright field images are corrected according to the following method:
%       Normalized Image = Original TIFF / Blank TIFF
% Polarized images are corrected according to the following method:
%       Normalized Image = Original TIFF - Blank TIFF
%
% Input 
% image_dir(optional): Provide path to folder containing images. User will be prompted to select folder if input is not provided
% 
% Output
% Background corrected histology images saved in image_dir with suffix "_C".  Original (non-background
% corrected) images are moved and saved to a subfolder called 'originals'
%
% Syntax
%   Correct_PONI_histology_images(image_dir)
%   Correct_PONI_histology_images
%
%% %%%%%%%%%%%%

if nargin == 0
    image_dir = uigetdir(pwd, 'Select folder containing images');
end

files = dir([image_dir '\*.TIF']); 

blank_files = strfind({files.name},'blank.TIF');
blank = ~cellfun(@isempty,blank_files); clear blank_files;
blank = find(blank == 1);

fn_blank = [image_dir '\' files(blank).name];

files = {files.name};
files(blank)=[];

for i = 1:size(files,2) 
    fn = [image_dir '\' files{i}];
    image_correction(fn_blank,fn)
end

movefile(fn_blank,[image_dir '\originals\']);
disp('Finished with background correction.')

end

%% %%%%%%%%%%%%%%%%%%%%%
 function image_correction(fn_blank,fn)
     BL = imread(fn_blank);
     try
        IM = imread(fn);
    catch
        disp([fn ': not available']);
    end
        
    if exist('IM','var')
        IM_norm = cast(IM,'double');
        BL_norm = cast(BL,'double');

        light_or_dark = mean(mean(mean(BL_norm(2:end-1,2:end-1,:)))); % Figure out whether blank slide is bright field or polarized

        if light_or_dark > (2^16)/3             % threshhold is set to 1/3 of full pixel range
            new_IM = cast((IM_norm ./ BL_norm) .* 2^16,'uint16');
        elseif light_or_dark < (2^16)/3
            new_IM = cast((IM_norm - BL_norm),'uint16');
        end

        [pth,nm] = fileparts(fn);
        if ~exist([pth '\originals'],'dir')
            mkdir([pth '\originals']);
        end
        movefile(fn,[pth '\originals\']);
         
        if ~exist([pth '\BkCorr'],'dir')
            mkdir([pth '\BkCorr']);
        end       
        imwrite(new_IM,[pth '\BkCorr\' nm '_C.TIF']);
        
    end
 end
