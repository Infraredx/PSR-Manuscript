function Coreg_BkgndCorrect_histology(image_dir)
%%
% Co-registers images acquired at different angles to the reference image, 0deg. 
% Performs file search for background corrected polarized images
% (filenames with *PR*_C*)
%
% Input 
% image_dir (optional): directory of PSR polarized background corrected images (TIF) 
% for same slide.  User will be prompted to select folder if input is not
% provided.
%
% Output
% Creates subfolder  "Rotated" "Composite_XX" where resulting co-registered image
% file is saved (suffix '_rotated'). Makes copy of 0deg iamge with suffix
% '_rotated'
%   
% Example Syntax 
%      overlay_histoimages_2deg('Sample1') 
%      overlay_histoimages_2deg
%
%% %%%%%%%%%

if nargin == 0
    image_dir = uigetdir(pwd, 'Select folder with Background corrected images');
end

flist = dir([image_dir '\*PR*deg_C*']); % find all TIF images;
flist = {flist.name};
if ~iscell(flist); flist = {flist}; end
flist = sort(flist);

save_dir = [fileparts(image_dir) '\Rotated']; %create subdirectory in root folder
if ~exist(save_dir,'dir'); mkdir(save_dir); end
        
for i = 1:length (flist)
    
    fn = flist{i};
    fn_comp = strsplit(fn,'_');
    rotangle = str2double((strrep(fn_comp{7},'deg', '')));
    disp(['rotation angle: ' num2str(rotangle)])

    if rotangle == 0
        IMfix = imread([image_dir '\' flist{i}]);
        imwrite (IMfix,[save_dir '\' strrep(flist{i}, '.TIF','_rotated.TIF')])
        IMfix = rgb2gray(IMfix);
    else
        IMmove = imread([image_dir '\' flist{i}]);
        
        %Pre-rotate to speed transformation search
        IMmove_rotated = imrotate(rgb2gray(IMmove),-1*rotangle); 
        
        %Determine transform matrix
        [optimizer, metric] = imregconfig('multimodal');
        optimizer.MaximumIterations = 1000;
        optimizer.Epsilon = 0.01;
        optimizer.InitialRadius = 0.001;
        optimizer.GrowthFactor = 1.05;
        metric.NumberOfSpatialSamples = 750;
        metric.NumberOfHistogramBins = 248;
        metric.UseAllPixels = 1;
        
        tform = imregtform(IMmove_rotated, IMfix, 'rigid', optimizer, metric);
        IMmove_registered = imwarp(IMmove_rotated, tform, 'OutputView',imref2d(size(IMfix)));
        
        %Perform transform on RGB image
        IMmove_registered = imwarp(imrotate(IMmove, -1*rotangle), tform,'OutputView',imref2d(size(IMfix)));
        
        fn_rotated = strrep(flist{i}, '.TIF','_rotated.TIF');
        imwrite(IMmove_registered,[save_dir '\' fn_rotated]);
        
    end
end

disp('Finished with image co-registration with 0deg.')


