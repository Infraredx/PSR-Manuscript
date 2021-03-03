function Coreg_composite(image_dir,refangle,fn_ref)
%%
% Co-registers composite images in image_dir.
% 
% 
% Input 
% image_dir (optional): directory of PSR polarized background corrected images (TIF) 
% for same slide.  User will be prompted to select folder if input is not
% provided.
% refangle(optional) = reference angle of  identify reference image that all other images will
% be co-registered to. Default = 60;  
% fn_ref(optional) = option to provide reference image used for co-registration of composite
% image
%
% Output
%   
% Example Syntax 
%      Coreg_composite(image_dir,60)
%      Coreg_composite(image_dir)
%
%% %%%%%%%%%

if ~exist('image_dir','var')
    image_dir = uigetdir(pwd, 'Select folder with Background corrected images');
end

if ~exist('refangle','var')
    refangle = 60;
end

flist = dir([image_dir '\*PR*deg_C*']); % find all TIF images;
mov_dir = [image_dir '\Nocoreg']; 
if ~exist(mov_dir,'dir'); mkdir(mov_dir); end

if ~exist('fn_ref','var')
    % find reference angle from image_dir
    strangle = ['_' num2str(refangle) 'deg_'];
    refidx = strfind({flist.name},strangle);
    refidx = ~cellfun(@isempty,refidx); 
    refidx = find(refidx == 1);

    if length(refidx)==1
        fn_ref = [image_dir '\' flist(refidx).name];
    else
        disp('Unable to find correct reference angle image'); return
    end
    
    % copy original reference image to subdirectory
    flist = {flist.name};
    IMfix = imread(fn_ref);
    movefile(fn_ref, [mov_dir '\' strrep(flist{refidx},'max.TIF', 'max_orig.TIF')]);
    imwrite(IMfix,fn_ref); 
    IMfix = rgb2gray(IMfix);
    
    flist(refidx)=[];
else
    % skip copy of original reference image to subdirectory
    flist = {flist.name};
    IMfix = rgb2gray(imread(fn_ref));
  
end
 

for i = 1:length (flist)
        fn = [image_dir '\' flist{i}];
        IMmoveRGB = imread(fn);
        IMmove = rgb2gray(IMmoveRGB);
       
        %Determine transform matrix
        [optimizer, metric] = imregconfig('multimodal');
        optimizer.MaximumIterations = 1500;
        optimizer.Epsilon = 0.005;
        optimizer.InitialRadius = 0.001;
        optimizer.GrowthFactor = 1.05;
        metric.NumberOfSpatialSamples = 750;
        metric.NumberOfHistogramBins = 248;
        metric.UseAllPixels = 1;
        
        tform = imregtform(IMmove, IMfix, 'rigid', optimizer, metric);
        IMmove_registered = imwarp(IMmove, tform, 'OutputView',imref2d(size(IMfix)));
    
        %Perform transform on RGB image
        IMmove_registered = imwarp(IMmoveRGB, tform,'OutputView',imref2d(size(IMfix)));
        
        movefile(fn,[mov_dir '\' strrep(flist{i},'max.TIF', 'max_orig.TIF')]);
        imwrite(IMmove_registered,fn);
        
end

disp('Finished with image co-registration')


