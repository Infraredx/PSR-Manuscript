function Call_ImageAnalysis(image_dir,contour_fn,pxthresh)
% Performs pixel intensity and pixel counting analysis on composite image
% files.
%
% Input
% image_dir (optional): directory of composite PSR images. user will be prompt if not
% provided
% contour_fn (optinoal): corresponding contour mat file. User will be prompt if not
% provided
% pxthresh (optinal) = pixel intneisty threshold for pixel coutning analysis. Default = 600.
%
% Output
% mat file with pixel indices of cap region, pixel intensity, pixel count
% and pixel intensity threshold for pixel counting (pxthresh). mat file
% will have suffix corresponding to pxthresh value and saved in image
% subdirectory called 'Analysis'
%
% Syntax
%  Call_ImageAnalysis
%  Call_ImageAnalysis(image_dir)
%  Call_ImageAnalysis(image_dir,contour_fn)
%  Call_ImageAnalysis(image_dir,contour_fn,pxthresh)
%% %%%%%

if ~exist('image_dir','var')
    image_dir = uigetdir('Select image directory for pixel intensity and pixel coutning analysis');
else
    if ~isfolder(image_dir)
        disp('Directory with composite images not found.'); return
    end
end

if ~exist('contour_fn','var')
    contour_fn= uigetfile ('*.mat','Select correpsonding contour mat file');
else
    if ~isfile(contour_fn)
        disp('Could not find contour mat file.'); return
    end
end

if ~exist('pxthresh','var')
    pxthresh = 600;
end
   
flist = dir([image_dir '\*composite_max.TIF']);
flist = {flist.name};

save_dir = [image_dir '\Analysis']; 
if ~exist(save_dir,'dir'); mkdir(save_dir); end

for i = 1:length(flist)
    fn = [image_dir '\' flist{i}];
    if isfile(fn)
        Output = Comp_ImageAnalysis(fn,contour_fn,pxthresh);
        
        Pxlidx = Output.PixelIndices;
        PxlInt = Output.PixelValues;
        PxlCnt = Output.numofpix;
        
        save_fn = strrep(flist{i},'composite_max.TIF', [num2str(pxthresh) '.mat']);
        save([save_dir '\' save_fn], 'Pxlidx', 'PxlInt', 'PxlCnt', 'pxthresh'); 
        clear Pxl*        
    else
        disp(['Composite image file not found: ' flist{i}]); return
    end
   
end

disp('Image analysis complete.')    
