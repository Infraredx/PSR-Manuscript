function Make_composite(image_dir,delta_angle,angle_range)
%% 
% Combines individual images based on delta angle to create composite image.
% Assumes images have been rotated and co-registered to 0 deg 
% (associated filename suffix:  _C_rotated.TIF) 
% 
% Input
% image_dir: directory of PSR polarized co-registered images (TIF) for same slide 
% delta_angle (optional): rotation angle increment.  Default=60  if
% values is not provided
% delta_angle = rotation angle step for composite formation. Default =60
% angle_range = range of rotation. Acceptable values
% 1= [0 360] (default)
% 2=[0 178]
% 3=[180 358]
% [start end]
%
% Output
% Creates subfolder  "Composite" where resulting composite image is saved
%
% Example syntax 
%      compositeoptions('C:\Test Data\TestOriginals\Rotated',2) 
%      compositeoptions('C:\Test Data\TestOriginals\Rotated')
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    image_dir = uigetdir(pwd, 'Select folder containing Rotated Co-registered images');
    delta_angle = 60;
    angle_range = [0 358];
elseif nargin==1
    delta_angle = 60; 
    angle_range = [0 358];
elseif nargin == 2
    angle_range = [0 358];
end

if  isequal(angle_range, [0 358]) || isequal(angle_range,1)
    angle = (0:delta_angle:358); 
    amin=0; amax=358;
elseif isequal(angle_range, [0 178]) ||isequal(angle_range, 2)
    angle = (0:delta_angle:178); 
    amin=0; amax=178;
elseif isequal(angle_range, [180 358]) ||isequal(angle_range,3)
    angle = (180:delta_angle:358); amin=180; amax=358;
else
    if length(angle_range)==2
        angle = (angle_range(1):delta_angle:angle_range(2));
        amin=angle_range(1); amax=angle_range(2);
    else
        disp ('Range of angles invalid.'); return
    end
end

disp(['Delta Angle = ' num2str(delta_angle)])
disp(['Range of rotation = [' num2str(amin) '-' num2str(amax) ']']);

files = dir([image_dir '\*PR*_C_rotated*']); 
files = {files.name}; 

% make list of files for composite and check that all files exists
fn0 = strfind(files,'_0deg_C_rotated');
fn0 = ~cellfun(@isempty,fn0); 
fn0 = find(fn0 == 1);
flist{1}=files{fn0};

for i = 2:length(angle)
        fn = strrep(flist{1},'_0deg_',['_' num2str(angle(i)) 'deg_']) ;

        if exist([image_dir '\' fn], 'file')    
            flist{i} = fn;
        else
            disp('File missing.Unable to make composite. '); return
        end     
end
      
% grab image files to make composite image
for i = 1:length(flist)     
         IM = imread([image_dir '\' flist{i}]);
         IM_max(:,:,:,i) = uint16(IM); 
end

IM_max = max(IM_max,[],4);

save_dir = [fileparts(image_dir) '\Composites_' num2str(amin) '_' num2str(amax) ]; %create subdirectory in root folder
if ~exist(save_dir,'dir'); mkdir(save_dir); end

IM_max_fn = strrep(flist{1},'0deg_C_rotated',[num2str(delta_angle) 'deg_composite_max']);
imwrite(uint16(IM_max),[save_dir '\' IM_max_fn]);

disp('Finished making composite image.' )
disp(['Composite saved in ' save_dir])
close all




