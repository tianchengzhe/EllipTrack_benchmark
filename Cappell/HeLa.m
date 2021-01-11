function Timelapse_universal_tif( movie_id, num_frames )
%% define basic stuffs
maskwrite = 1;
num_channels = 1;
channel_names = {'H2B'};
ringcalc = 0;

addpath('Functions/');
image_dir = ['Z:/projects/tracking_code/submission2/ctc_hela/training_raw/', sprintf('%02d', movie_id), '/'];

%% define mask directory, load cmos offset, define dir of bias
base_output_dir = 'Z:/projects/tracking_code/submission2/ctc_hela/mingyu/';
cmosoffset_dir = '';
bias_dir = '';

% mask dir
mask_dir = [base_output_dir, 'mask/'];
if (~exist(mask_dir, 'dir') && maskwrite)
    mkdir(mask_dir);
end

% cmosoffset
try
    h = load([cmosoffset_dir, 'cmosoffset.mat']); cmos = h.cmosoffset;
catch
    disp('No cmosoffset is found. Use 0 instead.'); cmos = 0;
end

% bias
bias = cell(num_channels, 1);
for i=1:num_channels
    try
        h = load([bias_dir, channel_names{i}, '_1.mat']); bias{i} = h.bias;
    catch
        disp(['Bias for ', channel_names{i}, ' is not found. Use 1 instead.']); bias{i} = 1;
    end
end

%% general settings that usually do not need changes%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namenucedge = 'nucedge_';

%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%
nucr = 24;  %nucr=6 for 10X and 2x2 bin
debrisarea = 100;
boulderarea = 1500;
blobthreshold = -0.03;

%%% tracking parameters %%%%%%%%%%%%%%%%%%%%%%%%%
masschangethreshold = 0.60;
areachangethreshold = 0.60;
daughtervariance = 0.10;
jitterframes = [];
maxjump = nucr*4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames = 1:num_frames;
blocksize = 10000;
maxcellnum = blocksize;
parameternum = 4 + 4*sum(ringcalc==0) + 9*sum(ringcalc==1); 
entry_names = {'nuc_center_x', 'nuc_center_y', 'nuc_area', 'nuc_mass'};
for i=1:num_channels
    entry_names = cat(2, entry_names, ...
        {[channel_names{i}, '_nuc_mean'], ...
        [channel_names{i}, '_nuc_median'], ...
        [channel_names{i}, '_nuc_75th'], ...
        [channel_names{i}, '_nuc_sum']});
    if (ringcalc(i) == 1)
        entry_names = cat(2, entry_names, ...
            {[channel_names{i}, '_cytoring_mean'], ...
            [channel_names{i}, '_cytoring_median'], ...
            [channel_names{i}, '_cytoring_fgmedian'], ...
            [channel_names{i}, '_cytoring_75th'], ...
            [channel_names{i}, '_cytoring_sum']});
    end
end
timetotal = tic;

%% analyze the image sequence
% well and site position
shot = num2str(movie_id);

% if tracedata already exist, do not redo
if exist([base_output_dir, 'tracedata_', shot, '.mat'], 'file')
    return;
end

% set empty matrixes for data storage
tracedata = nan(maxcellnum,num_frames,parameternum);
tracking = nan(maxcellnum,5);
badframes = nan(num_frames,1);
jitters = zeros(num_frames,2);

% read the files
raw = imread([image_dir, 't000.tif']);
raw = cat(3, raw, imread([image_dir, 't001.tif']));
raw = cat(3, raw, imread([image_dir, 't002.tif']));
for i=1:3
    raw(:, :, i) = max(raw(:, :, i) - min(min(raw(:, :, i))), 1);
end

try
    [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes]=timelapsesetup_nd2(raw,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
catch
    return;
end

% loop through frames for this well_site
for i=1:num_frames
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;        

    %% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % actully reading images
    raw = cell(num_channels, 1);
    for j=1:num_channels
        raw{j} = max((double(imread([image_dir, sprintf('t%03d.tif', f-1)])) - cmos)./bias{j}, 1);
        raw{j} = max(raw{j} - min(min(raw{j})), 1);
    end

    %% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == firstgoodindex
        firstsegmethod = 'single';
        switch firstsegmethod
            case 'log'
                nuc_mask = blobdetector_4(log(raw{1}),nucr,blobthreshold,debrisarea);
            case 'single'
                blurradius = 3;
                nuc_mask = threshmask(raw{1},blurradius);
                nuc_mask = markershed(nuc_mask,round(nucr*2/3));
            case 'double'
                blurradius = 3;
                nuc_mask = threshmask(raw{1},blurradius);
                nuc_mask = markershed(nuc_mask,round(nucr*2/3));
                nuc_mask = secondthresh(raw{1},blurradius,nuc_mask,boulderarea*2);
        end
        foreground = nuc_mask;
        nuc_mask = bwareaopen(nuc_mask,debrisarea);
        
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        eccentricitythresh = 0.85;
        nuc_mask = excludelargeandwarped_3(nuc_mask,boulderarea,eccentricitythresh);
    else
        nuc_mask = threshmask(raw{1},1);
        
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe = find(badframes==0,1,'last');
        if ismember(f,jitterframes)
            [reljitx,reljity]=registerimages(imfill(extractmask,'holes'),nuc_mask);
        else
            reljitx = 0;
            reljity = 0;
        end
        jitters(f,:) = jitters(lastgoodframe,:)+[reljitx,reljity];
        secondsegmethod = 'apriori';
        switch secondsegmethod
            case 'log'
                nuc_mask = blobdetector_4(log(raw{1}),nucr,blobthreshold,debrisarea);
            case 'double'
                nuc_mask = threshmask(raw{1},blurradius);
                nuc_mask = markershed(nuc_mask,round(nucr*2/3));
                nuc_mask = secondthresh(raw{1},blurradius,nuc_mask,boulderarea*2);
            case 'apriori'
                nuc_mask = apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
        end
        foreground = nuc_mask;
        nuc_mask = bwareaopen(nuc_mask,debrisarea);
    end
    
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask = imclearborder(nuc_mask);

    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression = 4;
    nanmask = imdilate(foreground,strel('disk',nucr/2));
    nanmaskcyto = imdilate(foreground,strel('disk',nucr*2));
    for j=1:num_channels
        try
            if (ringcalc(j) == 0)
                raw{j} = max(bgsubmasked_global_2(imfilter(raw{j},fspecial('disk',3),'symmetric'), nanmask, 11, compression, 30), 1);
            elseif (ringcalc(j) == 1)
                raw{j} = max(bgsubmasked_global_2(imfilter(raw{j},fspecial('disk',3),'symmetric'), nanmaskcyto, 11, compression, 30), 1);
            end
        catch
        end
    end

    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells] = bwlabel(nuc_mask);
    nuc_info = struct2cell(regionprops(nuc_mask,raw{1},'Area','Centroid','MeanIntensity')');
    nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
    
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc = median(nuc_area);
    if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f) = 1;
        badextractmask = bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[mask_dir,shot,'_',namenucedge,num2str(f),'bad.tif']);
        end
        continue;
    end
    blurthreshhigh = 1.2*mednuc;
    blurthreshlow = 0.8*mednuc;
    numthresh = 0.5*numcells;
    nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density = squeeze(cell2mat(nuc_info(3,1,:)));
    
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass = nuc_density.*nuc_area;
    curdata = [nuc_center(:,1:2),nuc_area,nuc_mass];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        nuc_center(:,1) = nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2) = nuc_center(:,2)+jitters(f,2);
        
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        debugpackage = {extractmask,jitters(lastgoodframe,:), [reljitx,reljity]};
        
        %%% track & correct merges (update centers, masses and labels) %%%%
        trackparams = {nucr,maxjump,debrisarea,masschangethreshold,areachangethreshold,daughtervariance};
        [tracedata,curdata,tracking,nuc_label] = adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,raw{1},nuc_label,jitters(f,:),trackparams,debugpackage);
        badframes(f)=0;
    end
    
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask = bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(nuc_label),[mask_dir,shot,'_',namenucedge,num2str(f),'.tif']);
    end
    
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid = find(~isnan(curdata(:,1)));
    nuc_center = curdata(cellid, 1:2);
    num_live_cells = length(cellid);
    nuc_info = regionprops(nuc_label, 'PixelIdxList');
    
    % compute ring info
    ring_label = cell(num_channels, 1); ring_info = cell(num_channels, 1);
    for j=1:num_channels
        if (ringcalc(j) == 1)
            ring_label{j} = getcytoring_3(nuc_label, 4, raw{j});
            ring_info{j} = regionprops(ring_label{j}, 'PixelIdxList');
        end
    end
    
    % for following stuffs, always arranged in: 
    % nuc_mean, nuc_median, nuc_75th, nuc_sum
    % cytoring_mean, cytomedian, cytoring_fgmedian, cytoring_75th, cytoring_sum, if ringcalc(j)=1
    signal_per_channel = cell(num_channels, 1);
    for n=1:num_live_cells
        curr_cell_id = cellid(n);
        for j=1:num_channels
            curr_signal = [ nanmean(raw{j}(nuc_info(curr_cell_id).PixelIdxList)), ...
                nanmedian(raw{j}(nuc_info(curr_cell_id).PixelIdxList)), ...
                quantile(raw{j}(nuc_info(curr_cell_id).PixelIdxList), 0.75), ...
                nansum(raw{j}(nuc_info(curr_cell_id).PixelIdxList)) ];
            
            if (ringcalc(j) == 1)
                if curr_cell_id > numel(ring_info{j})
                    curr_signal = cat(2, curr_signal, nan(1, 5));
                else
                    ring_all = raw{j}(ring_info{j}(curr_cell_id).PixelIdxList);
                    ring_all = ring_all(ring_all <= quantile(ring_all, 0.98));
                    if (sum(ring_all > 100) >= 30) % 100 and 30 are thresholds
                        ring_foreground = ring_all(ring_all > 100);
                    elseif (length(ring_all) > 30)
                        ring_foreground = ring_all;
                    else
                        ring_foreground = [];
                    end

                    curr_signal = cat(2, curr_signal, [ ...
                        nanmean(raw{j}(ring_info{j}(curr_cell_id).PixelIdxList)), ...
                        nanmedian(raw{j}(ring_info{j}(curr_cell_id).PixelIdxList)), ...
                        nanmedian(ring_foreground), ...
                        quantile(ring_all, 0.75), ...
                        sum(raw{j}(ring_info{j}(curr_cell_id).PixelIdxList))]);
                end
            end
            signal_per_channel{j} = cat(1, signal_per_channel{j}, curr_signal);
        end
    end
    
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata(cellid,f,:) = [ curdata(cellid, 1:4), cell2mat(signal_per_channel') ];
    if maxcellnum-max(cellid) < blocksize
        tempdata = nan(blocksize,num_frames,parameternum);
        temptrack = nan(blocksize,5);
        tracedata = [tracedata;tempdata];
        tracking = [tracking;temptrack];
        maxcellnum = maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
[tracedata,genealogy,jitters] = postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);

%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([base_output_dir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters','entry_names');
toc(timetotal);

end
