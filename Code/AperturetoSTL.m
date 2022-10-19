close all;
clear all;

% NATURAL UPPER SURFACE (Select 1/uncomment if using natural preloaded surfaces)
%load('S2_upper_nogouge_3_Height_remout_strong_denoised_level.mat');
%s_upper = Sq;

% NATURAL LOWER SURFACE (Select 1/uncomment if using natural preloaded surfaces)
%load('S2_lower_nogouge_3_Height_remout_strong_denoised_level.mat');
%s_lower = Sq;

% SYNTHETIC UPPER SURFACES (Select 1/uncomment if using synthetic preloaded surfaces)
 load('Same_SyntheticAperture_H=0.9_V1_M=0.5mm_Disp=0mm_rms=1.1mm_fullsize_s_upper.mat')
% load('Same_SyntheticAperture_H=0.75_V1_M=0.5mm_Disp=0mm_rms=1.1mm_fullsize_s_upper.mat')

% SYNTHETIC LOWER SURFACES (Select 1/uncomment if using synthetic preloaded surfaces)
 load('Same_SyntheticAperture_H=0.9_V1_M=0.5mm_Disp=0mm_rms=1.1mm_fullsize_s_lower.mat')
% load('Same_SyntheticAperture_H=0.75_V1_M=0.5mm_Disp=0mm_rms=1.1mm_fullsize_s_lower.mat')

%% SURFACE DISPLACEMENT PARAMETERS %%

% ApertureSeperation:    Aperture seperation in 'Z' direction (normal to surface).
% ApertureDisplacementX: Aperture displacement along surface in 'X' direction (long dimension).
%                        Positive value = synestral motion viewed from top surface.
%                        Negative value = dextral motion viewed from top surface.
% ApertureDisplacementY: Aperture displacement along surface in 'Y' direction (short dimension).

ApertureSeperation    = 0.50000;   % [mm]: Initial aperture width surface seperation
ApertureDisplacementX = 1.50000;   % [mm]: X-direction (along rows) surface displacement magnitude
ApertureDisplacementY = 0.00000;   % [mm]: Y-direction (along columns) surface displacement magnitude


%% SYNTHETIC SURFACE GENERATION PARAMETERS (if gen_syn_fracture selected)

RMS_roughness = 0.0011;    % [m]: RMS roughness of generated surface
Hurst_exp     = 0.80000;   % [-]: Hurst exponent of generated surface, must be decimal/integer in range 0-1
x_length      = 0.05000;   % [m]: X-direction (along rows) surface length (DEFAULT > 0.03 m if displacement used)
side_length   = 5000;      % [pixels]: side length of surface in pixels (DEFAULT > 3000 if displacement used)


%% USER FUNCTION SELECTION %%

surf_to_STL         = 0;    % [ON/OFF]: Write .Stl file upon code execution (RUNTIME INTENSIVE)
    filename        = 'ApertureFileName.stl'; % [string]: Name of output .Stl file if using surf_to_STL 

hurst_exp_calc      = 0;    % [ON/OFF]: Calculate Hurst exponent of surfaces (upper, lower, directional, mean) 
RMS_calc            = 0;    % [ON/OFF]: Calculate RMS roughness of surfaces (upper, lower, directional, mean) 

figures             = 0;    % [ON/OFF]: Print figures
        
% FRACTURE SURFACE SELECTION
surface_type        = 0;    % [ON/OFF]: 0: preloaded synthetic surface 
                            %           1: preloaded natural surface 
                            %           0: generated synthetic surface
gen_syn_fracture    = 0;    % [ON/OFF]: Create aperture with generated synthetic surfaces
use_upper_surface   = 0;    % [ON/OFF]: Create aperture with identical surfaces using the upper surface (s_upper)
use_lower_surface   = 1;    % [ON/OFF]: Create aperture with identical surfaces using the lower surface (s_lower)
detrend             = 1;    % [ON/OFF]: Detrend surfaces prior to aperture generation (reccommended)

% SELF-INTERSECTING FACE CORRECTIONS FOR STL WRITING 
asperity_correction = 1;    % [ON/OFF]: Asperity correction to remove self-interecting faces
    min_aperture    = 0.2;  % [mm]: minimum aperture allowed (for meshing purposes)
edge_corrections    = 0;    % [ON/OFF]: Edge correction to improve COMSOL boundary selection
    bevel_height    = 0.05; % [mm]: additional bevel on each surface edge
    bevel_mag       = 30;   % [rows]: number of rows of bevel

    
%% PREPROCESSING/DEBUGGING

if (ApertureDisplacementX == 0) && (ApertureDisplacementY == 0)
   asperity_correction = 0; 
end

%% SYNTHETIC FRACTURE GENERATION

if gen_syn_fracture == 1
    % Generate synthetic upper surface using specified parameters
    [s_upper,PixelWidth_upper,PSD_upper] = artificial_surf(RMS_roughness,Hurst_exp,x_length,side_length,side_length);
    % Generate synthetic lower surface using specified parameters
    [s_lower,PixelWidth_lower,PSD_lower] = artificial_surf(RMS_roughness,Hurst_exp,x_length,side_length,side_length);  
end

%% UNIT CONVERSION AND FRACTURE COMPOSITION

ApertureSeperation = abs(ApertureSeperation*10^-3);  % units converted to 'm'
ApertureDisplacementX = ApertureDisplacementX*10^-3;     % units converted to 'm'
ApertureDisplacementY = ApertureDisplacementY*10^-3;     % units converted to 'm'

% % Convert Z data from 'mm' to 'm' 
if surface_type == 1 
s_lower = s_lower(:,501:4000).*10^-3;
s_upper = s_upper(:,501:4000).*10^-3;

% Orient dimensions of upper and lower fracture surfaces
s_upper = flipdim(s_upper,1); 
s_lower = double(-1.00000*s_lower);
end

%% TREND REMOVAL OPERATION

if detrend == 1
s_upper = detrend_2d(s_upper);
s_lower = detrend_2d(s_lower); 
end

%% APERTURE SEPERATION

% This section seperates the upper fracture surface from the lower surface
% in the positive Z direction by 'ApertureSeperation' as defined by the user.

if use_lower_surface == 1
   s_upper = s_lower;
end

if use_upper_surface == 1
   s_lower = s_upper;
end

% Center both surfaces on 0 in X-Y plane
s_upper = s_upper - mean(mean(s_upper));
s_lower = s_lower - mean(mean(s_lower));

s_upper_saved = s_upper;

 C(1:size(s_upper,1),1:size(s_upper,2)) = double(ApertureSeperation);
 s_upper = double(s_upper + C);
 
 clear C;
 
%% X-Y COMPOSITE FRACTURE SURFACE DISPLACEMENT

% This section displaces the upper and lower fracture surfaces by
% magnitude 'ApertureDisplacementX' [mm] in the slip-parallel (X) direction, and
% magnitude 'ApertureDisplacementY' [mm] in the slip-parallel (Y) direction.

fprintf(['Aperture slip X-direction: ', num2str(ApertureDisplacementX*10^3),' mm.\n'])
fprintf(['Aperture slip Y-direction: ', num2str(ApertureDisplacementY*10^3),' mm.\n'])
fprintf(['Total slip magnitude: ', num2str(sqrt((ApertureDisplacementY*10^3)^2+(ApertureDisplacementX*10^3)^2)),' mm.\n'])

if ApertureDisplacementX || ApertureDisplacementY ~= 0

    % converts displacement from meters to columns
    ApertureDisplacementX = ApertureDisplacementX*1*10^5;
    ApertureDisplacementY = ApertureDisplacementY*1*10^5;
    
    if ApertureDisplacementX > 0
        s_upper = s_upper(:,1:size(s_upper,2)-ApertureDisplacementX);
        s_lower = s_lower(:,(ApertureDisplacementX+1):size(s_lower,2));
    end
    
    if ApertureDisplacementX < 0
        ApertureDisplacementX = abs(ApertureDisplacementX);
        s_upper = s_upper(:,(ApertureDisplacementX+1):size(s_upper,2));
        s_lower = s_lower(:,1:size(s_lower,2)-ApertureDisplacementX);
    end
    
    if ApertureDisplacementY > 0
        s_upper = s_upper((ApertureDisplacementY+1):size(s_upper,1),:);
        s_lower = s_lower(1:(size(s_lower,1)-ApertureDisplacementY),:);
    end
    
     if ApertureDisplacementY < 0
        ApertureDisplacementY = abs(ApertureDisplacementY);
        s_upper = s_upper(1:(size(s_upper,1)-ApertureDisplacementY),:);
        s_lower = s_lower((ApertureDisplacementY+1):size(s_lower,1),:);
    end
    
end

%% APERTURE WIDTH CALCULATION

ApertureWidth(1:size(s_upper,1),1:size(s_upper,2)) = 0;
T(1:size(s_upper,1),1:size(s_upper,2)) = 0;
K(1:size(s_upper,1),1:size(s_upper,2)) = 0;

ApertureWidth = s_upper-s_lower; 
ApertureWidth(ApertureWidth<0) = 0; % [m]: Aperture region with min aperture = 0 
AW_intersect = s_upper-s_lower; % [m]: Aperture region with negative aperture where intersections exist 

%% TARGET APERTURE SELECTION

if surface_type == 1 % for natural surfaces, select region of interest
startpoint = 1;    
endpoint = 3000;
end

if surface_type == 0 % for synthetic surfaces, select region of interest
startpoint = 1001;  
endpoint = 4000;
end

aperture_width_select = ApertureWidth(startpoint:endpoint,startpoint:endpoint).*10^3; % [mm]: Targeted aperture region with min aperture = 0
aperture_width_select_intersect = AW_intersect(startpoint:endpoint,startpoint:endpoint).*10^3; % [mm]: Targeted aperture region with aperture intersections
target_upper = s_upper(startpoint:endpoint,startpoint:endpoint).*10^3; % [mm]: targeted upper fracture surface
target_lower = s_lower(startpoint:endpoint,startpoint:endpoint).*10^3; % [mm]: targeted lower fracture surface 

%% MEAN APERTURE CALCULATION

mean_aperture = mean(mean(aperture_width_select)); % units: [mm] mean aperture pre-asperity correction

%% ASPERITY CORRECTIONS
% This section removes asperities and replaces them with smoothed planar regions as
% to avoid self-intersecting faces in .stl file

if asperity_correction == 1
    
    mean_upper = mean(mean(target_upper));
    mean_lower = mean(mean(target_lower));
    
for i = 1:1:size(aperture_width_select,1)
    for j = 1:1:size(aperture_width_select,2)
        if  (aperture_width_select_intersect(i,j) <= 0)
            
            target_upper(i,j) = target_upper(i,j) + ((abs(aperture_width_select_intersect(i,j)))/2);         
            target_lower(i,j) = target_lower(i,j) - ((abs(aperture_width_select_intersect(i,j)))/2);
 
        end
    end 
end

target_upper = target_upper + (min_aperture/2);
target_lower = target_lower - (min_aperture/2);

mean_aperture_corrected =  mean(mean(target_upper-target_lower)); % units: [mm] mean aperture post-asperity correction

end

%% EDGE CORRECTIONS

bevel_scalar = (1:1:bevel_mag)./bevel_mag;

for i = 1:(3000-(2*bevel_mag))
bevel_scalar_upper(1:bevel_mag,i) = bevel_height.*bevel_scalar';
bevel_scalar_lower(1:bevel_mag,i) = bevel_height.*bevel_scalar';
end

bevel_scalar_upper = rot90(bevel_scalar_upper,2);
bevel_scalar_left = rot90(bevel_scalar_upper,1);
bevel_scalar_right = rot90(bevel_scalar_upper,3);


if edge_corrections == 1
target_upper(1:bevel_mag,(bevel_mag+1):(3000-bevel_mag)) = target_upper(1:bevel_mag,(bevel_mag+1):(3000-bevel_mag)) + bevel_scalar_upper;
target_upper(((3000-bevel_mag)+1):3000,(bevel_mag+1):(3000-bevel_mag)) = target_upper(((3000-bevel_mag)+1):3000,(bevel_mag+1):(3000-bevel_mag)) + bevel_scalar_lower;
target_upper((bevel_mag+1):(3000-bevel_mag),1:bevel_mag) = target_upper((bevel_mag+1):(3000-bevel_mag),1:bevel_mag) + bevel_scalar_left;
target_upper((bevel_mag+1):(3000-bevel_mag),((3000-bevel_mag)+1):3000) = target_upper((bevel_mag+1):(3000-bevel_mag),((3000-bevel_mag)+1):3000) + bevel_scalar_right;

bevel_scalar = [[0] bevel_scalar];
corner_scalar(1:bevel_mag,1:bevel_mag) = 0;
bevel_scalar = flip(bevel_scalar);

% CORNER CORRECTIONS
for i = 2:bevel_mag
corner_scalar(1:i,1:i) = corner_scalar(1:i,1:i) + (bevel_scalar(i).*bevel_height) - (bevel_scalar(i-1).*bevel_height);
end
corner_scalar(1,1) = -bevel_height;

corner_scalar_up_right = corner_scalar+(bevel_height*(1/0.99));
corner_scalar_up_left = fliplr(corner_scalar_up_right);
corner_scalar_down_right = flipud(corner_scalar_up_right);
corner_scalar_down_left = flipud(fliplr(corner_scalar_up_right));

% Corner corrections added to upper
target_upper(1:bevel_mag,1:bevel_mag) = target_upper(1:bevel_mag,1:bevel_mag) + corner_scalar_down_left;
target_upper(((3000-bevel_mag)+1):3000,1:bevel_mag) = target_upper(((3000-bevel_mag)+1):3000,1:bevel_mag) + corner_scalar_up_left;
target_upper(1:bevel_mag,((3000-bevel_mag)+1):3000) = target_upper(1:bevel_mag,((3000-bevel_mag)+1):3000) + corner_scalar_down_right;
target_upper(((3000-bevel_mag)+1):3000,((3000-bevel_mag)+1):3000) = target_upper(((3000-bevel_mag)+1):3000,((3000-bevel_mag)+1):3000) + corner_scalar_up_right;

% Corner corrections added to lower
target_lower(1:bevel_mag,1:bevel_mag) = target_lower(1:bevel_mag,1:bevel_mag) - corner_scalar_down_left;
target_lower(((3000-bevel_mag)+1):3000,1:bevel_mag) = target_lower(((3000-bevel_mag)+1):3000,1:bevel_mag) - corner_scalar_up_left;
target_lower(1:bevel_mag,((3000-bevel_mag)+1):3000) = target_lower(1:bevel_mag,((3000-bevel_mag)+1):3000) - corner_scalar_down_right;
target_lower(((3000-bevel_mag)+1):3000,((3000-bevel_mag)+1):3000) = target_lower(((3000-bevel_mag)+1):3000,((3000-bevel_mag)+1):3000) - corner_scalar_up_right;

target_lower(1:bevel_mag,(bevel_mag+1):(3000-bevel_mag)) = target_lower(1:bevel_mag,(bevel_mag+1):(3000-bevel_mag)) - bevel_scalar_upper;
target_lower(((3000-bevel_mag)+1):3000,(bevel_mag+1):(3000-bevel_mag)) = target_lower(((3000-bevel_mag)+1):3000,(bevel_mag+1):(3000-bevel_mag)) - bevel_scalar_lower;
target_lower((bevel_mag+1):(3000-bevel_mag),1:bevel_mag) = target_lower((bevel_mag+1):(3000-bevel_mag),1:bevel_mag) - bevel_scalar_left;
target_lower((bevel_mag+1):(3000-bevel_mag),((3000-bevel_mag)+1):3000) = target_lower((bevel_mag+1):(3000-bevel_mag),((3000-bevel_mag)+1):3000) - bevel_scalar_right;

end

%% Surface conversion from .MAT to .STL for COMSOL geometry input (RUNTIME/MEMORY INTENSIVE)

if surf_to_STL == 1
   
    surface_length = 30; % side length of export aperture, units: [mm]
    x = 0.01:0.01:surface_length;
    y = 0.01:0.01:surface_length;
    
    [F,V] = surf2solid(x,y,target_upper,'elevation',target_lower);
    stlwrite(filename,F,V);   
        
end

%% Hurst Exponent calculation (RUNTIME INTENSIVE)

if hurst_exp_calc == 1

    H_vec_upper_X(1:size(target_upper,1)) = 0;
    H_vec_lower_X(1:size(target_lower,1)) = 0;
    H_vec_upper_Y(1:size(target_upper,1)) = 0;
    H_vec_lower_Y(1:size(target_lower,1)) = 0;
    
    for i = 1:size(target_upper,1)
       H_vec_upper_X(i) = genhurst(target_upper(i,:)); 
       H_vec_lower_X(i) = genhurst(target_lower(i,:)); 
       
       H_vec_upper_Y(i) = genhurst(target_upper(:,i)); 
       H_vec_lower_Y(i) = genhurst(target_lower(:,i)); 
    end

    H_upper_X = mean(H_vec_upper_X); H_upper_Y = mean(H_vec_upper_Y);
    H_lower_X = mean(H_vec_lower_X); H_lower_Y = mean(H_vec_lower_Y);
    
    H_mean_X = mean([H_upper_X H_lower_X]); % Mean Hurst exponent in X direction (along rows)
    H_mean_Y = mean([H_upper_Y H_lower_Y]); % Mean Hurst exponent in Y direction (along columns)

end

%% RMS Roughness calculation (RUNTIME INTENSIVE)

if RMS_calc == 1

    RMS_vec_upper_X(1:size(target_upper,1)) = 0;
    RMS_vec_lower_X(1:size(target_lower,1)) = 0;
    RMS_vec_upper_Y(1:size(target_upper,1)) = 0;
    RMS_vec_lower_Y(1:size(target_lower,1)) = 0;
    
    for i = 1:size(target_upper,1)
       RMS_vec_upper_X(i) = rms(target_upper(i,:)); 
       RMS_vec_lower_X(i) = rms(target_lower(i,:)); 
       
       RMS_vec_upper_Y(i) = rms(target_upper(:,i)); 
       RMS_vec_lower_Y(i) = rms(target_lower(:,i)); 
    end

    RMS_upper_X = mean(RMS_vec_upper_X);  RMS_upper_Y = mean(RMS_vec_upper_Y);
    RMS_lower_X = mean(RMS_vec_lower_X);  RMS_lower_Y = mean(RMS_vec_lower_Y);
    
    RMS_mean_X = mean([RMS_upper_X RMS_lower_X]); % Mean RMS roughness in X direction (along rows)
    RMS_mean_Y = mean([RMS_upper_Y RMS_lower_Y]); % Mean RMS roughness in Y direction (along columns)
    
end

%% Figure production

if figures == 1
    
figure
hold on
surf(target_upper,'edgecolor','none','FaceColor','interp','FaceLighting','gouraud');
surf(target_lower,'edgecolor','none','FaceColor','interp','FaceLighting','gouraud');
title(['[Aperture] (H = ',num2str(mean([H_mean_X H_mean_Y])),' RMS = ',num2str(mean([RMS_mean_X RMS_mean_Y])),')']);
xticks([0,500,1000,1500,2000,2500,3000]);
xticklabels({'0','5','10','15','20','25','30'});
yticks([0,500,1000,1500,2000,2500,3000]);
yticklabels({'0','5','10','15','20','25','30'});
pbaspect([1 1 1/4]);
%axis equal
hold off

figure
hold on
contourf(target_upper-target_lower);
title('Aperture Width [mm]');
xticks([0,500,1000,1500,2000,2500,3000]);
xticklabels({'0','5','10','15','20','25','30'});
yticks([0,500,1000,1500,2000,2500,3000]);
yticklabels({'0','5','10','15','20','25','30'});
xlabel('X [mm]')
ylabel('Y [mm]')
axis equal
hold off

end