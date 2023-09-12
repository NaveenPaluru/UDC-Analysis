
% Sample run file for computing porosity and thickness of the 
% under deposit corrosion from optical microscopy (OM) images.

clc; clear all; close all;

% load the input OM image
I = imread('image.tif');

% assign the resolution
pixel = 0.1107;  % (μm)

% segmented image
[S,centroids] = imsegkmeans(I,3);       

% getting tube deposit
indices = [1;2;3];
template = [round(mean(centroids,2)) indices];
sortedtemplate = sortrows(template);
idx = sortedtemplate(2,2);
T = zeros([size(I,1), size(I,2)]); T(S==idx)=1;

% holes filling
H1 = logical(imfill(T));
H2 = bwareafilt(H1,[0.1 15000]);
H  = H1 & (~H2);            

% extracting pores globally
C = H & logical(T);
P = xor(C,H);
vv = sum(P>0);
vt = sum(C>0);
porosity = vv/vt;
tortuosity = 1 + 0.8 * (1-porosity);

% extracting pores locally
mask = ones(55,55);
nr = conv2(double(P), mask, 'same');
C = or(H,P);
dr = conv2(double(C), mask, 'same');         
locpor = nr./dr;
locpor(isnan(locpor)) = 0;
locpor(isinf(locpor)) = 0;
locpor(P==0) = 0;
avg_locpor = mean(locpor(P>0));   
avg_loctor = 1 + 0.8 * (1-avg_locpor);

% thickness estimation        
temp1 = sum(H,1);
loc1  = find(temp1==max(temp1));      
for j=1:length(loc1)        
    temp2 = H(:,loc1(j)); [x,y] = find(temp2==1);
    dist(j) = (x(end)-x(1));
end        
loc2 = find(dist==max(dist));
loc2 = loc2(1);
temp3 = H(:,loc1(loc2)); [x,y] = find(temp3==1);        
thickness = (x(end)-x(1))*pixel;   

% layerwise porosity
num_pores_row = sum(P,2);
row_porosity  = sum(locpor,2);        
layer_porosity = row_porosity./num_pores_row;     
lp = [layer_porosity (1:size(layer_porosity,1))'];
pltdata1 = lp(~isnan(lp(:,1)),:);
pltdata2  = pltdata1(~isinf(pltdata1(:,1)),:);
pltdata2(:,2) = pltdata2(:,2)-min(pltdata2(:,2));
            
% saving the results 
mkdir('Results');
cd('./Results/');

figure('Visible', 'off');      
imshow(I);
imwrite(getframe(gcf).cdata, 'original.png');
close all;

figure('Visible', 'off'); 
imshow(mat2gray(H));  
imwrite(getframe(gcf).cdata, 'deposit.png');
close all;            

figure('Visible', 'off');
imshow(mat2gray(P));
imwrite(getframe(gcf).cdata, 'pores.png');
close all;

figure('Visible', 'off');
imagesc(locpor); colorbar; axis off;
imwrite(getframe(gcf).cdata, 'locpores.png');
close all;

figure('Visible', 'off');
coord = [loc1(loc2) x(1) loc1(loc2) x(end)];
RGB = insertShape(I, 'Line', coord ,'LineWidth', 10 , 'Color', { 'blue'});
thickness = round(thickness,2);
RGB = insertText(RGB , [2000, 1700] , strcat(num2str(thickness), ' µm'),  'FontSize', 100, 'TextColor', 'blue', 'BoxOpacity', 0);
imshow(RGB); 
imwrite(getframe(gcf).cdata, 'thickness.png');
close all;

figure('Visible', 'off');
plot(pltdata2(:,2)*pixel, flipud(smooth(pltdata2(:,1)*100)), 'm-', 'LineWidth', 1.5);
xlabel('Thickness (µm)', 'FontSize', 20);
ylabel('Porosity   (%)', 'FontSize', 20);
ax = gca;
ax.FontSize = 16; 
grid on;
title('Porosity VS Thickness');
imwrite(getframe(gcf).cdata, 'porosity_VS_thickness.png');
close all;            

logdata{3,1} = strcat('Porosity : ', num2str(porosity));
logdata{5,1} = strcat('Local Porosity : ', num2str(avg_locpor));
logdata{7,1} = strcat('Tortuosity : ', num2str(tortuosity));
logdata{9,1} = strcat('Local Tortuosity : ', num2str(avg_loctor));
logdata{11,1} = strcat('Thickness : ', num2str(thickness));
writecell(logdata,'result.csv');
cd('../');
   
        
      














