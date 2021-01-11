function Timelapse_universal_nd2( row_id, col_id, site_id )
%% define basic stuffs
maskwrite = 1;
num_channels = 2;
channel_names = {'CFP', 'mCherry'};
ringcalc = [0, 1];
num_frames = 288; % start from frame 1, until frame num_frames

addpath('Functions/');

base_folders = {'Z:/projects/tracking_code/20181115/Raw/'};

% get the readers for nd2 files
num_folders = length(base_folders);
nd2_frame_range = nan(num_folders, 2);
all_bfReader = cell(num_folders, 1);
leading_filename = ['Well', char(row_id-1+'A'), sprintf('%02d', col_id)];
for i=1:num_folders
    % get bfReader
    all_filenames = dir(base_folders{i});
    recorded_fileid = [];
    for k=1:length(all_filenames)
        if (length(all_filenames(k).name) > 8 && strcmp(all_filenames(k).name(1:length(leading_filename)), leading_filename))
            recorded_fileid = k;
            break;
        end
    end
    if (isempty(recorded_fileid))
        error('ND2 file not found.');
    end
    all_bfReader{i} = BioformatsImage([base_folders{i}, all_filenames(recorded_fileid).name]);
    
    % get frame ranges
    if (i==1)
        nd2_frame_range(i,1) = 1;
    else
        nd2_frame_range(i,1) = nd2_frame_range(i-1,2)+1;
    end
    nd2_frame_range(i,2) = nd2_frame_range(i,1) - 1 + all_bfReader{i}.sizeT;
end

%% define mask directory, load cmos offset, define dir of bias
base_output_dir = 'Z:/projects/tracking_code/submission2/my_MCF10A/mingyu/';
cmosoffset_dir = 'Z:/microscope_mat_files/nikon1_matfiles/';
bias_dir = 'Z:/microscope_mat_files/nikon1_matfiles/';

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
        h = load([bias_dir, channel_names{i}, '.mat']); bias{i} = h.bias;
    catch
        disp(['Bias for ', channel_names{i}, ' is not found. Use 1 instead.']); bias{i} = 1;
    end
end

%% general settings that usually do not need changes%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namenucedge = 'nucedge_';

%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%
nucr = 12;  %nucr=6 for 10X and 2x2 bin
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
shot = [num2str(row_id), '_', num2str(col_id), '_', num2str(site_id)];

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
raw = all_bfReader{1}.getXYplane(1, site_id, 1);
raw = cat(3, raw, all_bfReader{1}.getXYplane(1, site_id, 2));
raw = cat(3, raw, all_bfReader{1}.getXYplane(1, site_id, 3));

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
    % select the right reader and frame id
    reader_id = find(nd2_frame_range(:,1) <= f & nd2_frame_range(:,2) >= f);
    frame_within_reader = f - nd2_frame_range(reader_id, 1) + 1;

    % actully reading images
    raw = cell(num_channels, 1);
    for j=1:num_channels
        raw{j} = max((double(all_bfReader{reader_id}.getXYplane(j,site_id,frame_within_reader)) - cmos)./bias{j}, 1);
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
            imwrite(uint16(badextractmask),[mask_dir,'\',shot,'_',namenucedge,num2str(f),'bad.tif']);
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
        imwrite(uint16(nuc_label),[mask_dir,'\',shot,'_',namenucedge,num2str(f),'.tif']);
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
