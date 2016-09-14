function qt1_vfa_lin_fit(img, alpha, TR, outT1, outM0, B1img, maskimg)
% Calculates T1 map using DESPOT1
% see Deoni, Peters and Rutt, Magn Reson Med 53: 237-241 (2005)
%
%   img1            :   minc file: SPGR (FLASH) image with alpha1
%   img2            :   minc file: SPGR image with alpha2
%   alpha1          :   excitation flip angle for first image in degrees
%   alpha2          :   excitation flip angle for second image in degrees
%   B1img           :   B1 map mnc file
%   TR              :   repetition time, constant for both images
%   output          :   name for output file


%  October 16th 2013: Fixed NAN T1 output
%  June 18th 2014: Removed hack for single slice images at the eof,
%  replaced niak_write_minc with niak_write_minc_ss which is a modified
%  file to handle 2D images. See in mboudrea-phd-code/MATLAB/general/ or
%  the mril3 share folder.
%

% test if output file already exists
fid = fopen(outT1,'r');
if(fid ~= -1)
  fprintf('\nError T1map: cannot overwrite T1map output file %s\n', ...
      outT1);
  fclose(fid);
  return
end

fid = fopen(outM0,'r');
if(fid ~= -1)
  fprintf('\nError T1map: cannot overwrite M0map output file %s\n', ...
      outM0);
  fclose(fid);
  return
end

if (length(img) ~= length(alpha))
  fprintf('\nError T1map: img and alpha are not the same length \n');
  fclose(fid);
  return
end    

% open images
h = cell(length(img));
h_head = cell(length(img));

for i = 1:length(img)
    [h_head{i},h{i}] = niak_read_minc(img{i});
end


W = h_head{1}.info.dimensions(1);
L = h_head{1}.info.dimensions(2);
N = h_head{1}.info.dimensions(3);
fprintf('Image Width %d\n', W)
fprintf('Image Height %d\n', L) 
fprintf('Image Slices %d\n', N) 

if (~isempty(B1img))
  [h_b1_hdr,h_b1] = niak_read_minc(B1img);
  
  if ( (h_b1_hdr.info.dimensions(1) ~= W) || (h_b1_hdr.info.dimensions(2) ~= L) || (h_b1_hdr.info.dimensions(3) ~= N))
    fprintf('\nError T1map: img and B1map must have the same dimensions \n');
    %fclose(fid);
    return
  end 
end;
   
if (~isempty(maskimg))
  [h_mask_hdr,h_mask] = niak_read_minc(maskimg);
  
  if ( (h_mask_hdr.info.dimensions(1) ~= W) || (h_mask_hdr.info.dimensions(2) ~= L) || (h_mask_hdr.info.dimensions(3) ~= N))
    fprintf('\nError T1map: img and mask must have the same dimensions \n');
    %fclose(fid);
    return
  end 
end;



%------------------
% calculate T1 map
%------------------

T1 = zeros(W*L,N);
M0 = zeros(W*L,N);


for n = 1:N
    
    disp(n)
    
    for i = 1:length(img)
        temp_h = h{i}(:,:,n);
        image(i,:) = temp_h(:);
    end
    
    if (~isempty(B1img))
        temp_h_b1 = h_b1(:,:,n); 
        b1 = temp_h_b1(:);    
    else
        b1 = ones(W*L,1);
    end;
    
    if (~isempty(maskimg))
        temp_h_mask = h_mask(:,:,n);
        mask = temp_h_mask(:);    
    else
        mask = ones(W*L,1);
    end;
    
    
    for i = 1:length(img)
        alpha_corr(i,:) = alpha(i).*b1;
    end  
    
    y = image./sind(alpha_corr);
    x = image./tand(alpha_corr);
        
    for w = 1:1:W*L
            
        if (mask(w))
            p = polyfit(x(:,w), y(:,w), 1);

            T1(w,n) = -TR/log(p(1));
            
            M0(w,n) = p(2)/(1-p(1));
        else
            T1(w,n) = 0;
            M0(w,n) = 0;
        end	
        
    end

end

bounds = [0.0 5.0];
T1 = T1 + (T1 > bounds(2)).*(bounds(2)-T1) + (T1 < bounds(1)).*(bounds(1)-T1);

bounds = [0.0 20000.0];
M0 = M0 + (M0 > bounds(2)).*(bounds(2)-M0) + (M0 < bounds(1)).*(bounds(1)-M0);


%----------------------------------
% save T1 and M0 map to minc format
%----------------------------------

T1image = zeros(W,L,N);
M0image = zeros(W,L,N);
for ii=1:N
    T1image(:,:,ii) = reshape(T1(:,ii),W,L);
    M0image(:,:,ii) = reshape(M0(:,ii),W,L);
end

T1image(isnan(T1image))=0;

t1hdr = niak_read_hdr_minc(img{1});
t1hdr.file_name = outT1;
niak_write_minc_ss(t1hdr,T1image);

t1hdr.file_name = outM0;
niak_write_minc_ss(t1hdr,M0image);

return;


