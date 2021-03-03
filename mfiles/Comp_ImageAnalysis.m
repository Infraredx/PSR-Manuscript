function OUT = Comp_ImageAnalysis(IMfile,contour_fn,pxthresh)
%
% Uses contours from pathologist (saved as mat) to compute collect pixel
% intensity values and pixel counts above thereshold in cap region under
% necrotic core. Analysis performed at 1 degree increment from 12 o'clock
% position in contour file.
% Uses dependent function "intersections" (Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz)
%
% Inputs
% Imfile= Composite PSR image(TIF). 
% contour_fn= corresponding contour mat file. 
% pxthresh = pixel intneisty threshold for pixel coutning analysis 
%
% Output
% OUT - structure with  cap pixel indices, intensity and pixel count at
% 1 degree increment.

%% %%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocate fields
contourfields = {'Lipid','Necro','Calcified'};
CTOR.PixelIndices = cell(1,3);
CTOR.PixelIndices{1} = cell(1,360);
CTOR.PixelIndices{2} = cell(1,360);
CTOR.PixelIndices{3} = cell(1,360);
CTOR.PixelContour = cell(1,3);
CTOR.PixelContour{1} = cell(1,360);
CTOR.PixelContour{2} = cell(1,360);
CTOR.PixelContour{3} = cell(1,360);
numofpix=zeros(360,1);
PixelValues = zeros(360,1);
L1x = cell(1,3);
L1x{1} = cell(1,360);
L1x{2} = cell(1,360);
L1x{3} = cell(1,360);
L1y = cell(1,3);
L1y{1} = cell(1,360);
L1y{2} = cell(1,360);
L1y{3} = cell(1,360);
Fx{1} = cell(1,360);
Fy{1} = cell(1,360);

CT = load(contour_fn);
CT.Lumen = [CT.Lumen; CT.Lumen(1,:)];
CTfields = fieldnames(CT);

if ~~isfield(CT,'Catheter')
    if ~isequal(CT.Catheter,[0 0])
        CTcent0temp = CT.cent0;
        CT.cent0 = CT.Catheter;
    end
end

new_top12 = CT.top12 - CT.cent0;
[theta,~] = cart2pol(new_top12(1),new_top12(2)); % Find 12 o'clock angle
theta = theta * 180/pi; %degrees

Comp_IM = imread(IMfile);
rho = max(size(Comp_IM));
next_theta = (0:359) + theta;
[x,y] = pol2cart(next_theta .* (pi/180),rho); % back to radians
x = x + CT.cent0(1); y = y + CT.cent0(2); % These are the coordinates of the 1 degree ray
[X2,Y2] = meshgrid(1:size(Comp_IM,2),1:size(Comp_IM,1));

tempNaN = cell(1,length(contourfields));

% Limit analysis to cap under necrotic core
for gg = 2 
    contour = strfind(CTfields,contourfields{gg});
    contour = ~cellfun(@isempty,contour);
    cnames = CTfields(contour);
    for hh = 1:length(cnames) % Loop through all contours of the same type
        if ~~isfield(CT,cnames{hh}) && ~isequal(CT.(cnames{hh}),[0 0])
            temp = [CT.(cnames{hh}); CT.(cnames{hh})(1,:)]; % close contour line back to beginning of contour.
            tempNaN{gg} = [tempNaN{gg}; NaN, NaN; temp];
            clear temp;
        end
    end
    clear contour cnames;
end

% Send out rays at 1 degree intervals
for i = 1:361 
    ff = mod(i-1,360) + 1; %current index "i";
    
    % Find where the 1 degree ray intersects the lumen contour
    [Fx{ff},Fy{ff}] = intersections([CT.cent0(1) x(ff)],[CT.cent0(2) y(ff)],CT.Lumen(:,1),CT.Lumen(:,2));
    
    % Find the closest lumen contour out from the lumen centroid. 
    [~,idx] = min(pdist2([CT.cent0(1),CT.cent0(2)],[Fx{ff},Fy{ff}],'euclidean'));
    clear Fx_d
    
    Fx_d = [Fx{ff}(idx), Fy{ff}(idx)]; %Lumen coordinate 
    Fx{ff} = Fx{ff}(idx);
    Fy{ff} = Fy{ff}(idx);
    clear idx;
    
    % Find lumen boundary in 1 degree segment
    if i > 1 
        onedeg_boundary = [CT.cent0(1), CT.cent0(2);x(ff) y(ff); x(i-1) y(i-1); CT.cent0(1), CT.cent0(2)];
        
        %find closest point to the 1 degree line boundaries intersection.
        [~,F1idx] = min(pdist2([Fx{i-1}, Fy{i-1}],[CT.Lumen(:,1),CT.Lumen(:,2)],'euclidean'));
        [~,F2idx] = min(pdist2([Fx{ff}, Fy{ff}],[CT.Lumen(:,1),CT.Lumen(:,2)],'euclidean'));
                
        if F2idx > F1idx && diff([F2idx,F1idx]) < length(CT.Lumen)/4
            [F2idx,F1idx] = deal(F1idx,F2idx);
        end
    end
    
    % Limit analysis to cap under necrotic core
    for gg = 2 
        edge_flag = 0;
        if ~isempty(tempNaN{gg})
            
            %Find intersection of the ray and the lipid contour
            [tempL1x,tempL1y] = intersections([CT.cent0(1) x(ff)],[CT.cent0(2) y(ff)],tempNaN{gg}(:,1),tempNaN{gg}(:,2)); 
            if i > 1
                [prevL1x,prevL1y] = intersections([CT.cent0(1) x(i-1)],[CT.cent0(2) y(i-1)],tempNaN{gg}(:,1),tempNaN{gg}(:,2)); 
            else
                prevL1x = []; prevL1y = [];
            end
            
            if (~isempty(tempL1x) && ~isempty(tempL1y)) || (~isempty(prevL1x) && ~isempty(prevL1y))
                
                all_L1d = pdist2(Fx_d,[tempL1x,tempL1y],'euclidean');
                if ~isempty(all_L1d)
                    [~, idx] = min(all_L1d(1:length(tempL1x)));
                    L1x{gg}{ff} = tempL1x(idx);
                    L1y{gg}{ff} = tempL1y(idx);
                    clear tempL1x tempL1y;
                end
                
                if i > 1 && (~isempty(L1x{gg}{ff}) || ~isempty(L1x{gg}{i-1}))
                    % Find lipid boundary of the 1 degree segment
                    
                    if ~isempty(L1x{gg}{ff})
                        [~,L1idx] = min(pdist2([L1x{gg}{ff},L1y{gg}{ff}],[tempNaN{gg}(:,1),tempNaN{gg}(:,2)],'euclidean'));
                    else
                        % situation where the 1 degree segment is on the edge of a contour (and the partial segment would otherwise be skipped)
                        p = [CT.cent0(1), CT.cent0(2); x(ff), y(ff)];
                        [~,L1idx] = orthoDistance(tempNaN{gg},p);
                        L1x{gg}{ff} = tempNaN{gg}(L1idx,1);
                        L1y{gg}{ff} = tempNaN{gg}(L1idx,2);
                        clear p;
                        edge_flag = 1;
                    end
                    
                    if ~isempty(L1x{gg}{i-1})
                        [jnk,L2idx] = min(pdist2([L1x{gg}{i-1},L1y{gg}{i-1}],[tempNaN{gg}(:,1),tempNaN{gg}(:,2)],'euclidean'));
                        clear jnk;
                    else
                        % situation where the 1 degree segment is on the edge of a contour (and the partial segment would otherwise be skipped)
                        p = [CT.cent0(1), CT.cent0(2); x(i-1), y(i-1)];
                        [~,L2idx] = orthoDistance(tempNaN{gg},p);
                        L1x{gg}{i-1} = tempNaN{gg}(L2idx,1);
                        L1y{gg}{i-1} = tempNaN{gg}(L2idx,2);
                        clear p;
                        edge_flag = 2;
                    end
                    
                    % for edge effect issues
                    clear Fx2 Fy2;
                    if ~~edge_flag
                        if edge_flag == 1
                            [Fx2,Fy2] = intersections([CT.cent0(1) L1x{gg}{ff}],[CT.cent0(2) L1y{gg}{ff}],CT.Lumen(:,1),CT.Lumen(:,2));
                        elseif edge_flag == 2
                            [Fx2,Fy2] = intersections([CT.cent0(1) L1x{gg}{i-1}],[CT.cent0(2) L1y{gg}{i-1}],CT.Lumen(:,1),CT.Lumen(:,2));
                        end
                        [~,idx] = min(pdist2([CT.cent0(1),CT.cent0(2)],[Fx2,Fy2],'euclidean'));
                        Fx2 = Fx2(idx);
                        Fy2 = Fy2(idx);
                        clear idx;
                    end
                    
                    %case if indices need to be reversed
                    if L2idx > L1idx 
                        [L2idx,L1idx] = deal(L1idx,L2idx);
                    end
                    
                    if numel(L2idx:L1idx) > length(tempNaN{gg})*0.5
                        fullidx = 1:length(tempNaN{gg});
                        notidx = L2idx:L1idx;
                        lip_idx = setdiff(fullidx,notidx);
                    else
                        lip_idx = L2idx:L1idx;
                    end
                    
                    if F2idx > F1idx && diff([F2idx,F1idx]) < length(CT.lumen)/4
                        [F2idx,F1idx] = deal(F1idx,F2idx);
                    end
                    
                    if numel(F2idx:F1idx) > length(CT.Lumen)*0.5
                        fullidx = 1:length(CT.Lumen);
                        notidx = F2idx:F1idx;
                        lum_idx = setdiff(fullidx,notidx);
                    else
                        lum_idx = F2idx:F1idx;
                    end
                    
                    %check to make sure contour lines are within the 1 degree boundary points:
                    Bidxy = contour_line_check(L1x{gg}{i-1},L1x{gg}{ff},L1y{gg}{i-1},L1y{gg}{ff},(tempNaN{gg}(lip_idx,:)));
                    
                    % for edge effect issues
                    if edge_flag == 1
                        Cidxy = contour_line_check(Fx{i-1},Fx2,Fy{i-1},Fy2,(CT.Lumen(lum_idx,:)));
                    elseif edge_flag == 2
                        Cidxy = contour_line_check(Fx2,Fx{ff},Fy2,Fy{ff},(CT.Lumen(lum_idx,:)));
                    else
                        Cidxy = contour_line_check(Fx{i-1},Fx{ff},Fy{i-1},Fy{ff},(CT.Lumen(lum_idx,:)));
                    end

                    temp2 = tempNaN{gg}(lip_idx(Bidxy),:);
                    
                    inCheck = inpolygon(temp2(:,1),temp2(:,2),onedeg_boundary(:,1),onedeg_boundary(:,2));
                    temp2 = temp2(inCheck,:);
                    
                    if length(tempNaN{gg}(lip_idx,:)) < 15
                        sorted_data = data_sorting_distance(temp2,L1x{gg}{i-1},L1y{gg}{i-1});
                        
                        if edge_flag == 1
                            xv = [Fx{i-1}, L1x{gg}{i-1}, sorted_data(:,1)', L1x{gg}{ff}, Fx2];
                            yv = [Fy{i-1}, L1y{gg}{i-1}, sorted_data(:,2)', L1y{gg}{ff}, Fy2];
                        elseif edge_flag == 2
                            xv = [Fx2, L1x{gg}{i-1}, sorted_data(:,1)', L1x{gg}{ff}, Fx{ff}];
                            yv = [Fy2, L1y{gg}{i-1}, sorted_data(:,2)', L1y{gg}{ff}, Fy{ff}];
                        else
                            xv = [Fx{i-1}, L1x{gg}{i-1}, sorted_data(:,1)', L1x{gg}{ff}, Fx{ff}];
                            yv = [Fy{i-1}, L1y{gg}{i-1}, sorted_data(:,2)', L1y{gg}{ff}, Fy{ff}];
                        end
                        
                        clear sorted_data Liptemp sortval sortorder;
                    else
                        if edge_flag == 1
                            xv = [Fx{i-1}, L1x{gg}{i-1}, L1x{gg}{ff}, Fx2];
                            yv = [Fy{i-1}, L1y{gg}{i-1}, L1y{gg}{ff}, Fy2];
                        elseif edge_flag == 2
                            xv = [Fx2, L1x{gg}{i-1}, L1x{gg}{ff}, Fx{ff}];
                            yv = [Fy2, L1y{gg}{i-1}, L1y{gg}{ff}, Fy{ff}];
                        else
                            xv = [Fx{i-1}, L1x{gg}{i-1}, L1x{gg}{ff}, Fx{ff}];
                            yv = [Fy{i-1}, L1y{gg}{i-1}, L1y{gg}{ff}, Fy{ff}];
                        end
                        
                    end
                    
                    if length(CT.Lumen(lum_idx,:)) < 15
                        sorted_data = data_sorting_distance(CT.Lumen(lum_idx(Cidxy),:),Fx{ff},Fy{ff});    
                        xv = [xv, sorted_data(:,1)'];
                        yv = [yv, sorted_data(:,2)'];
                        clear sorted_data Lumtemp sortval sortorder;
                    end
                    
                    if ~isempty(L1x{gg}{i-1}) && ~isempty(L1x{gg}{ff}) && ~isempty(L1y{gg}{i-1}) && ~isempty(L1y{gg}{ff})
                        in = inpolygon(X2,Y2,xv,yv);
                                               
                        CTOR.PixelContour{gg}{ff} = [xv',yv'];
                        CTOR.PixelIndices{gg}{ff} = [X2(in),Y2(in)];
                    end
                end
            end
        end
    end
end
clear Fx_d;

% Plot and perform pixel intensity and counting in cap region
show_contours(IMfile,contour_fn);
[~,tmp] = fileparts(IMfile); title(strrep(tmp,'_', ' ')); 
Comp_IM = rgb2gray(Comp_IM);
gg = 2;

for bb = 1:360    
    if ~isempty(CTOR.PixelIndices{gg}{bb})
        
        plot([CTOR.PixelContour{gg}{bb}(:,1); CTOR.PixelContour{gg}{bb}(1,1)],[CTOR.PixelContour{gg}{bb}(:,2); CTOR.PixelContour{gg}{bb}(1,2)],'--');
      
        linearInd = sub2ind(size(Comp_IM),CTOR.PixelIndices{gg}{bb}(:,2),CTOR.PixelIndices{gg}{bb}(:,1));       
        PixelValues(bb,1) = sum(Comp_IM(linearInd));
        
        if any((Comp_IM(linearInd)) > pxthresh)
            newarray = Comp_IM(linearInd)>pxthresh;
            numofpix(bb,1) = sum(newarray(:)); 
        end
    end
end

OUT.PixelValues=PixelValues;
OUT.numofpix = numofpix;
OUT.PixelIndices = CTOR.PixelIndices{1,2}';
end


function Fxy = contour_line_check(PPx1,PPx2,PPy1,PPy2,comp)
% Check to see if points on a line lie between two points PPx,y
%% %%%%%%%%%%%%%%

if PPx1 > PPx2
    Bx1 = PPx2;
    Bx2 = PPx1;
else
    Bx2 = PPx2;
    Bx1 = PPx1;
end
if PPy1 > PPy2
    By1 = PPy2;
    By2 = PPy1;
else
    By2 = PPy2;
    By1 = PPy1;
end
Fxy = ((comp(:,1)' <= Bx2 & comp(:,1)' >= Bx1) | (comp(:,2)' <= By2 & comp(:,2)' >= By1));
end

function sortd = data_sorting_distance(L_points,Px,Py)
% Sort L_points with respect to distance from Px,Py.
%% %%%%%%%%%%%%%%

sortval = (L_points(:,1) - Px).^2 + (L_points(:,2) - Py).^2;
[~,sortorder] = sort(sortval);
sortd = L_points(sortorder,:);
end

function [d, idx] = orthoDistance(curve,p)
% Finds the point on curve closest to the line p. 
% First checks that point on the curve is in the correct quadrant and not 180 degrees away.
% p(1,:) is assumed to be lumen centroid
%% %%%%%%%%%%%%%%%

centroid = p(1,:);
line = p(2,:) - centroid; 
init_quad = (atan2(line(:,2),line(:,1))); 
avoid_quad = init_quad - pi;
curve2 = curve - repmat(centroid,size(curve,1),1);
curve_quad = (atan2(curve2(:,2),curve2(:,1)));
curve_quad(curve_quad < 0) = curve_quad(curve_quad < 0) + 2*pi; 
if avoid_quad < 0; avoid_quad = avoid_quad + 2*pi; end 
avoid_quad = floor(avoid_quad ./ (pi/2)) + 1; 
CC = floor(curve_quad ./ (pi/2)) + 1; 
curve(CC == avoid_quad) = NaN; 

% Now process as normal.
nPoints = size(curve, 1);
dist = zeros(nPoints, 1);
for i = 1:nPoints
    if ~isnan(curve(i,:))
        dist(i) = abs(det([p(2,:)-p(1,:);curve(i,:)-p(1,:)]))/norm(p(2,:)-p(1,:));
    else
        dist(i) = NaN;
    end
end
[d,idx] = min(dist);
end

function show_contours(IMfile,CTfile)
% Displays histology image (IMfile) and any/all of its
% associated contours from CTfile. 
    
if ~exist(IMfile,'file') || ~exist(CTfile,'file')
    disp('Files do not exist');
    return
end

CT = load(CTfile);
Img = imread(IMfile);
CTfields = fieldnames(CT);

figure, imshow(Img); hold all;
set(gca,'Position',[0.05 0.05 .9 .9]);
hold all;

for ff = 1:length(CTfields)
    try
        if ~~isfield(CT,CTfields{ff}) && isempty(find(CT.(CTfields{ff}) == 0, 1))
            temp = CT.(CTfields{ff});
            if size(temp,1) > 1
                temp = [temp; temp(1,:)];
            end
            h = plot(temp(:,1),temp(:,2));
            if length(temp) > 2
                set(h,'LineWidth',3);
            else
                if ~~strfind(CTfields{ff},'Catheter')
                    set(h,'Marker','o','Color','r','Linewidth',2,'MarkerSize',17,'MarkerFaceColor','y');
                else
                    set(h,'Marker','+','MarkerSize',12,'LineWidth',5);
                end
            end
            if ~~strfind(CTfields{ff},'Lipid')
                set(h,'Color','y')
            elseif ~~strfind(CTfields{ff},'Necro')
                set(h,'Color','b')
            elseif ~~strfind(CTfields{ff},'Calcified')
                set(h,'Color','m')
            end
            
            clear h temp;
        end
    catch
        disp(CTfields{ff});
    end
end

end