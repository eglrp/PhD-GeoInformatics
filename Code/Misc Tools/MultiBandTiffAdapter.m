% This undocumented class may be removed in a future release.

% TIFFADAPTER ImageAdapter for TIFF images.
%   ADAPTER = TIFFADAPTER(FILENAME,MODE) creates a TiffAdapter object,
%   ADAPTER, associated with the TIFF file, FILENAME.  TiffAdapter inherits
%   from the ImageAdapter class and provides an interface for reading and
%   writing regions to the TIFF file on disk.  The MODE argument can be:
%
%      'r'  : read-only
%      'r+' : read & write
%      'w'  : read & write to a new file (see syntax below)
%
%   ADAPTER = TIFFADAPTER(FILENAME,'w',IMAGE_SIZE,FILL_VALUE) creates a
%   TiffAdapter object associated with a new file, which is created during
%   the object construction.  The new file, FILENAME, has size IMAGE_SIZE
%   and data type as determined by the FILL_VALUE.  The TiffAdapter is
%   ready for writing after construction, but is not pre-filled with the
%   FILL_VALUE.  When the close() method is invoked any unwritten pixels
%   will be filled with the FILL_VALUE.
%
%   Tag Restrictions for Writing
%   ----------------------------
%   The TiffAdapter class cannot write to certain TIFF files.  Writable
%   TIFF files must satisfy the following TIFF Tag restrictions:
%
%      Compression         : Must be None.
%      PlanarConfiguration : Cannot be Separate.
%      Photometric         : Cannot be YCbCr.
%
%   Files that violate any of these conditions will be read-only when
%   accessed via the TiffAdapter class regardless of the value of the MODE
%   argument.
%
%   Data Types Supported for Creating New Images
%   --------------------------------------------
%   grayscale : any non-64-bit numeric type
%   rgb       : uint8 only
%
%   Data Sizes Supported for Creating New Images
%   --------------------------------------------
%   M-by-N or M-by-N-by-3
%
%   Examples
%   --------
%   This example uses the BLOCKPROC function to resize a matrix using a
%   TiffAdapter as a read-only input.
%
%       fun = @(block_struct) imresize(block_struct.data,0.15);
%       adpt = TiffAdapter('pout.tif','r');
%       I2 = blockproc(adpt,[100 100],fun);
%       adpt.close();
%       figure;
%       imshow('pout.tif');
%       figure;
%       imshow(I2);
%
%   See also BLOCKPROC, IMAGEADAPTER.

%   Copyright 2009-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.10.2.1 $  $Date: 2011/07/18 00:35:20 $

classdef MultiBandTiffAdapter < ImageAdapter
    
    properties (Constant = true)
        
        DefaultTileHeight = 512;   % tile height for writing
        DefaultTileWidth  = 512;   % tile width for writing
        
    end % constant public properties
    
    
    properties (Access = private)
        
        % File related properties
        Filename         % Filename
        Mode             % Mode
        
        % TIFF related properties
        TiffObj          % Tiff objects
        Photom           % Photometric Interpretation
        Datatype         % MATLAB data type
        
        UseImread        % true if IMREAD used for reading
        Info             % IMFINFO (for use with IMREAD)
        TileCache        % Tile cache object
        ReadOnly         % writable status for file
        
        BlockWidth       % Block/Tile Width
        BlockHeight      % Block/Tile Height
        
        % Properties for writing
        FillVal          % Fill value
        IsTileWritten    % Array of logical of number of tiles/blocks
        
    end % private properties
    
    
    methods
        
        %---------------------------------------------------
        function obj = MultiBandTiffAdapter(filename, mode, varargin)
            
%             validatestring(mode,{'r','r+','w'},'MultiBandTiffAdapter','Mode',2);
            validatestring(mode,{'w'},'MultiBandTiffAdapter','Mode',2); %make it write only
            
            obj.Filename = filename;
            obj.Mode = mode;
            obj.TileCache = internal.images.TileCache();
            
            if strcmpi(mode, 'w')
                % create a new file, discard if existing
                
                if nargin ~= 4
                    error(message('images:MultiBandTiffAdapter:invalidInputParamCount'));
                end
                
                % verify that our file is either MxN or MxNx3
                obj.ImageSize = varargin{1};
                if length(obj.ImageSize) < 2 || length(obj.ImageSize) > 3 %||...
                    error(message('images:MultiBandTiffAdapter:invalidDims'));
                end
                
                if (length(obj.ImageSize) == 2)
                    obj.ImageSize(3) = 1;
                end
%                 if (length(obj.ImageSize) == 3 && obj.ImageSize(3) == 1)
%                     obj.ImageSize = obj.ImageSize(1:2);
%                 end                
                
                % validate data type
                obj.FillVal = varargin{2};
                obj.Datatype = class(obj.FillVal);
                
                non_writable_types = {'int64','uint64'};
                if any(strcmpi(obj.Datatype,non_writable_types))
                    error(message('images:MultiBandTiffAdapter:invalidType', obj.Datatype));
                end
                
                %is this really necessary - can't we have single/double multi band images??
                is_rgb = (length(obj.ImageSize) == 3 && obj.ImageSize(3) == 3);
                if is_rgb
                    if isa(obj.FillVal,'single') || ...
                        isa(obj.FillVal,'double') || ...
                         isa(obj.FillVal,'logical') 
%                        error(message('images:MultiBandTiffAdapter:invalidType', obj.Datatype));
                    end                    
                end
                
                obj.create();
                
            elseif strcmpi(mode, 'r')

                % Read-only existing file.  To reduce the number of calls
                % to the mex routines, we use the 'PixelRegion' syntax of
                % IMREAD for read-only TIFF files instead of the Tiff
                % class.
                obj.TiffObj = [];
                obj.Info = imfinfo(obj.Filename);
                
                % use IMREAD
                obj.ReadOnly = true;
                obj.UseImread = true;
                
                % Get image size
                imageHeight = obj.Info.Height;
                imageWidth  = obj.Info.Width;
                imageChannels = obj.Info.SamplesPerPixel;
                if (imageChannels == 1)
                    obj.ImageSize = [imageHeight imageWidth];
                else
                    obj.ImageSize = [imageHeight imageWidth imageChannels];
                end
                
            else % mode == 'r+'
                
                % existing file
                obj.TiffObj = Tiff(obj.Filename, obj.Mode);
                obj.Info = imfinfo(obj.Filename);
                
                % get some properties of the TIFF file
                is_compressed = ~isequal(obj.TiffObj.getTag('Compression'),...
                    Tiff.Compression.None);
                is_planar_config_separate = isequal(obj.TiffObj.getTag(Tiff.TagID.PlanarConfiguration),...
                    Tiff.PlanarConfiguration.Separate);
                obj.Photom = obj.TiffObj.getTag(Tiff.TagID.Photometric);
                is_YCbCr = isequal(obj.Photom, Tiff.Photometric.YCbCr);
                
                % cache data type
                obj.Datatype = tags2Datatype(obj.TiffObj.getTag(Tiff.TagID.BitsPerSample),...
                    obj.TiffObj.getTag(Tiff.TagID.SampleFormat));
                
                % use object to read by default
                obj.UseImread = false;
                obj.ReadOnly = false;
                
                % fall back to read-only mode if we cannot write
                if is_compressed
                    warning(message('images:MultiBandTiffAdapter:CannotWriteCompressed'));
                    obj.ReadOnly = true;
                    % for performance, we also use IMREAD here...
                    obj.UseImread = true;
                end
                if is_planar_config_separate
                    warning(message('images:MultiBandTiffAdapter:planarConfigurationError'));
                    obj.ReadOnly = true;
                    obj.UseImread = true;
                end
                if is_YCbCr
                    warning(message('images:MultiBandTiffAdapter:photometricError'));
                    obj.ReadOnly = true;
                    obj.UseImread = true;
                end
                
                % Get image size
                imageHeight = obj.TiffObj.getTag(Tiff.TagID.ImageLength);
                imageWidth  = obj.TiffObj.getTag(Tiff.TagID.ImageWidth);
                imageChannels = obj.TiffObj.getTag(Tiff.TagID.SamplesPerPixel);
                if (imageChannels == 1)
                    obj.ImageSize = [imageHeight imageWidth];
                else
                    obj.ImageSize = [imageHeight imageWidth  imageChannels];
                end
                
                % get TIFF block size
                [obj.BlockWidth obj.BlockHeight] = obj.getBlockSize();
                
                % If using IMREAD, close the Tiff object.  We no longer
                % need it and it interferes with calls to IMREAD.
                if obj.UseImread
                    obj.TiffObj.close();
                end
                
                % Assume that all blocks are written
                [Nx Ny] = obj.getSizeInBlocks();
                obj.IsTileWritten = true(1, Nx * Ny);
                
            end
            
        end % Constructor
        
        
        %----------------------------------------------
        function result = readRegion(obj, start, count)
            
            % Get selected region
            if (start(1) + count(1) - 1 > obj.ImageSize(1) || ...
                    start(2) + count(2) - 1 > obj.ImageSize(2))
                error(message('images:MultiBandTiffAdapter:countOutOfRangeRead'));
            end
            
            if obj.UseImread
                % If imread is to be used
                result = imread(obj.Filename, 'Info', obj.Info, 'PixelRegion', ...
                    {[start(1), start(1) + count(1) - 1], ...
                    [start(2), start(2) + count(2) - 1]});
                return;
            end
            
            % Calculate block indexes at top-left corner of region
            [bStartX bStartY] = obj.getBlockIndex(start(2), start(1));
            % Calculate block indexes at bottom-right corner of region
            [bEndX bEndY] = obj.getBlockIndex(start(2) + count(2) - 1, ...
                start(1) + count(1) - 1);
            
            % iterate over these block indexes
            result = [];
            [Nx ~] = obj.getSizeInBlocks();
            for y = bStartY:bEndY
                % Top pixel coordinate
                py = (y - 1) * obj.BlockHeight + 1;
                
                for x = bStartX:bEndX
                    % Left pixel coordinate
                    px = (x - 1) * obj.BlockWidth + 1;
                    
                    % get Tile Handle object for blockId
                    blockId = (y - 1) * Nx + x;
                    th = obj.getTileHandle(blockId, true);
                    
                    if isempty(result)
                        % compute the region size (get num bands from tile)
                        rdims = size(th.tile);
                        rdims(1:2) = count(1:2);
                        
                        % preallocate region result
                        if islogical(th.tile)
                            result = false(rdims);
                        else
                            result = zeros(rdims, class(th.tile));
                        end
                    end
                    
                    % calculate valid pixels along Y
                    ry = max(py, start(1)) : min(py + size(th.tile, 1) - 1, ...
                        start(1) + count(1) - 1);
                    % calculate valid pixels along X
                    rx = max(px, start(2)) : min(px + size(th.tile, 2) - 1, ...
                        start(2) + count(2) - 1);
                    % copy valid pixels
                    result(ry - start(1) + 1, rx - start(2) + 1, :) = ...
                        th.tile(ry - py + 1, rx - px + 1, :); %#ok<AGROW>
                end
            end
            
            % If minimum is interpreted as white for logicals, complement
            % the output result
            if obj.Photom == Tiff.Photometric.MinIsWhite && ...
                    islogical(result)
                result = ~result;
            end
            
        end % readRegion
        
        
        %------------------------------------------
        function [] = writeRegion(obj, start, data)
            
            if obj.ReadOnly
                error(message('images:MultiBandTiffAdapter:readOnly'));
            end
            
%                 (length(obj.ImageSize) == 3 && length(size(data)) == 3 && obj.ImageSize(3) == size(data,3)),...
%             assert((length(obj.ImageSize) == 2 && size(data,3) == 1) || ...
%                 (length(obj.ImageSize) == 3 && obj.ImageSize(3) == size(data,3)),...
%                 'images:MultiBandTiffAdapter:invalidDataSize', ...
%                 'MultiBandTiffAdapter: Data size should match the file dimensions.');
            assert(obj.ImageSize(3) == size(data, 3),...
                'images:MultiBandTiffAdapter:invalidDataSize', ...
                'MultiBandTiffAdapter: Data size should match the file dimensions.');
            
            count = [size(data,1) size(data,2)];
            if any(count==0)
                return;
            end
            
            if (start(1) + count(1) - 1 > obj.ImageSize(1) || ...
                    start(2) + count(2) - 1 > obj.ImageSize(2))
                error(message('images:MultiBandTiffAdapter:countOutOfRangeWrite'));
            end
            
            [Nx] = obj.getSizeInBlocks();
            % Calculate block indexes at top-left corner of region
            [bStartX bStartY] = obj.getBlockIndex(start(2), start(1));
            % Calculate block indexes at bottom-right corner of region
            [bEndX bEndY] = obj.getBlockIndex(start(2) + count(2) - 1, ...
                start(1) + count(1) - 1);
            
            % If minimum is interpreted as white for logicals, complement
            % the input data
            if obj.Photom == Tiff.Photometric.MinIsWhite && ...
                    islogical(data)
                data = ~data;
            end
            
            % iterate over these block indexes
            for y = bStartY:bEndY
                py = (y - 1) * obj.BlockHeight + 1;
                blockId = (y - 1) * Nx + bStartX;
                
                for x = bStartX:bEndX
                    px = (x - 1) * obj.BlockWidth + 1;
                    
                    % check if reading block is necessary
                    if (py < start(1) || px < start(2) || ...
                            py + obj.BlockHeight > start(1) + count(1) || ...
                            px + obj.BlockWidth > start(2) + count(2))
                        
                        % read block/tile
                        th = obj.getTileHandle(blockId, obj.ReadOnly);

                        % calculate valid pixels along Y
                        ry = max(py, start(1)) : min(py + size(th.tile, 1) - 1, ...
                            start(1) + count(1) - 1);
                        % calculate valid pixels along X
                        rx = max(px, start(2)) : min(px + size(th.tile, 2) - 1, ...
                            start(2) + count(2) - 1);
                        % copy valid pixels
                        th.tile(ry - py + 1, rx - px + 1, :) = ...
                            data(ry - start(1) + 1, rx - start(2) + 1, :);
                        
                        % block/tile is automatically written py when
                        % th (tileHandle) goes out of scope from tileCache
                        
                    else
                        % calculate valid pixels along Y
                        ry = py : py + obj.BlockHeight - 1;
                        % calculate valid pixels along X
                        rx = px : px + obj.BlockWidth - 1;
                        
                        % since whole block/tile is written no caching will
                        % be done here for efficiency
                        obj.writeBlock(blockId, data(ry - start(1) + 1, ...
                            rx - start(2) + 1, :));
                        
                        % mark tile written
                        obj.IsTileWritten(blockId) = true;
                        
                    end
                    
                    blockId = blockId + 1;
                end
            end
            fprintf('w.');
        end % writeRegion
        
        
        %------------------
        function close(obj)
            
            % This makes sure that tiles are flushed before the TiffObj
            % closes
            obj.TileCache.close();
            
            % The TiffObj can be empty in some cases so we use the
            % functional calling syntax
            for i = 1:obj.ImageSize(3)
                close(obj.TiffObj{i});
            end
            
        end % close
        
        
        %---------------------------------------
        function block = readBlock(obj, blockId)
            if obj.TiffObj{1}.isTiled()
                block = obj.TiffObj{1}.readEncodedTile(blockId);
            else
                block = obj.TiffObj{1}.readEncodedStrip(blockId);
            end
            
        end % readBlock
        
        
        %--------------------------------------------
        function [] = writeBlock(obj, blockId, block)
            for i = 1:obj.ImageSize(3)
                if obj.TiffObj{i}.isTiled()
                    obj.TiffObj{i}.writeEncodedTile(blockId, block(:,:,i));
                else
                    obj.TiffObj{i}.writeEncodedStrip(blockId, block(:,:,i));
                end
            end
            
        end % writeBlock
        
    end % public methods
    
    
    methods (Access = private)
        
        %------------------------
        function [] = create(obj)
            obj.ReadOnly = false;
            
            %make one tiff for each band            
            for i = 1:obj.ImageSize(3)
                obj.TiffObj{i} = Tiff(sprintf('%s_Band%d.tif', obj.Filename(1:end-4), i), 'w');
            end

            % Select appropriate tile size
            tagStruct.TileWidth  = obj.DefaultTileWidth;
            tagStruct.TileLength = obj.DefaultTileHeight;
            obj.BlockWidth       = obj.DefaultTileWidth;
            obj.BlockHeight      = obj.DefaultTileHeight;
            
            % Set known image attributes
            tagStruct.ImageWidth = obj.ImageSize(2);
            tagStruct.ImageLength = obj.ImageSize(1);
            
            % Set Photometric Interpretation
            if obj.ImageSize(3) == 3
                % Unsigned datatypes require RGB Photometric Interp
                obj.Photom = Tiff.Photometric.RGB;
                if isa(obj.FillVal, 'int8') || ...
                    isa(obj.FillVal, 'int16') || ...
                     isa(obj.FillVal, 'int32')
                    % Signed datatypes require MinIsBlack Photometric Interp
                    obj.Photom = Tiff.Photometric.MinIsBlack;
                end
                tagStruct.Photometric = obj.Photom;
            elseif obj.ImageSize(3) < 3
                if islogical(obj.FillVal)
                    obj.Photom = Tiff.Photometric.MinIsWhite;
                    obj.FillVal = ~obj.FillVal;
                else
                    obj.Photom = Tiff.Photometric.MinIsBlack;
                end
                tagStruct.Photometric = obj.Photom;
            elseif obj.ImageSize(3) > 3
%                  obj.Photom = Tiff.Photometric.RGB;
%                  
                obj.Photom = Tiff.Photometric.MinIsBlack;  
%                 obj.Photom = Tiff.Photometric.RGB;
%                 tagStruct.Photometric = obj.Photom;
%                 xs(1:obj.ImageSize(3) - 3) = Tiff.ExtraSamples.Unspecified;
%                 tagStruct.ExtraSamples = xs;
                tagStruct.Photometric = obj.Photom;
            else
                error(message('images:MultiBandTiffAdapter:unsupportedNumChannels'));
            end
            tagStruct.SamplesPerPixel = 1; %prod(obj.ImageSize) / prod(obj.ImageSize(1:2));
            
            % set datatype
            dtype = class(obj.FillVal);
            [bits_per_sample sample_format] = datatype2Tags(dtype);
            tagStruct.BitsPerSample = bits_per_sample;
            tagStruct.SampleFormat = sample_format;
            
            % Set default parameters
            tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagStruct.Compression = Tiff.Compression.None;
            
            % set all tags
            for i = 1:obj.ImageSize(3)
                obj.TiffObj{i}.setTag(tagStruct);
            end            
            
            % mark all tiles unwritten
            [Nx Ny] = obj.getSizeInBlocks();
            obj.IsTileWritten = false(1, Nx * Ny);
            % Create the Tiff object.
%             obj.TiffObject = Tiff(obj.Filename, 'w');
            
            % Setup the tiff file properties
            % See "Exporting Image Data and Metadata to TIFF files
            % http://www.mathworks.com/help/techdoc/import_export/f5-123068.html#br_c_iz-1
            %
%             obj.TiffObject.setTag('ImageLength',   obj.ImageSize(1));
%             obj.TiffObject.setTag('ImageWidth',    obj.ImageSize(2));
%             obj.TiffObject.setTag('TileLength',    obj.TileLength);
%             obj.TiffObject.setTag('TileWidth',     obj.TileWidth);
%             obj.TiffObject.setTag('Photometric',   Tiff.Photometric.RGB);
%             obj.TiffObject.setTag('BitsPerSample', 8);
%             obj.TiffObject.setTag('SampleFormat',  Tiff.SampleFormat.UInt);
%             obj.TiffObject.setTag('SamplesPerPixel', 3);
%             obj.TiffObject.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);            
        end % create
        
        
        %--------------------------------------------
        function [bWidth bHeight] = getBlockSize(obj)
            if obj.TiffObj.isTiled()
                bWidth = obj.TiffObj{1}.getTag(Tiff.TagID.TileWidth);
                bHeight = obj.TiffObj{1}.getTag(Tiff.TagID.TileLength);
            else
                bWidth = obj.TiffObj{1}.getTag(Tiff.TagID.ImageWidth);
                bHeight = obj.TiffObj{1}.getTag(Tiff.TagID.RowsPerStrip);
            end
            
        end % getBlockSize
        
        
        %--------------------------------------
        function [Nx Ny] = getSizeInBlocks(obj)
            % Taking care of rounding issues
            Nx = floor((obj.ImageSize(2) - 1) / obj.BlockWidth) + 1;
            Ny = floor((obj.ImageSize(1) - 1) / obj.BlockHeight) + 1;
        end
        
        
        %--------------------------------------------------------
        function [indX indY] = getBlockIndex(obj, pixelX, pixelY)
            % Taking care of rounding issues
            indX = floor((pixelX - 1) / obj.BlockWidth) + 1;
            indY = floor((pixelY - 1) / obj.BlockHeight) + 1;
        end
        
        
        %------------------------------------------------------
        function th = getTileHandle(obj, blockId, isReadRegion)
            
            if isReadRegion
                
                % The only case where getTileHandle is called from a
                % readRegion method is when the adapter was opened in 'r+'
                % mode.  Pure 'r'-mode readRegions are handled through
                % IMREAD.
                if obj.IsTileWritten(blockId)
                    % if tile is written, check the cache and read
                    % block from disk if necessary.
                    th = obj.TileCache.getTile(blockId);
                    if isempty(th)
                        th.tile = obj.readBlock(blockId);
                    end
                else
                    % if tile is not written, make a buffer of FillVal
                    % and associate with a tile handle so that it gets
                    % written eventually. Also add new TileHandle to
                    % the cache.
                    tiledims = obj.ImageSize;
                    tiledims(1:2) = [obj.BlockHeight obj.BlockWidth];
                    tileBuf = repmat(obj.FillVal, tiledims);
                    th = internal.images.TileHandle(obj, blockId, obj.ReadOnly, tileBuf);
                    obj.TileCache.setTile(th);
                    
                    % Set tile written to be true
                    obj.IsTileWritten(blockId) = true;
                end
                
            else
                
                % writeRegion function call
                if obj.IsTileWritten(blockId)
                    % if tile is written, check the cache and if
                    % necessary create a tile handle object. Also add
                    % the tileHandle to the cache.
                    th = obj.TileCache.getTile(blockId);
                    if isempty(th)
                        th = internal.images.TileHandle(obj, blockId, obj.ReadOnly);
                        obj.TileCache.setTile(th);
                    end
                else
                    % if tile is not written, make a buffer of FillVal
                    % and associate with a tile handle so that it gets
                    % written eventually. Also add new tileHandle to
                    % the cache.
                    tiledims = obj.ImageSize;
                    tiledims(1:2) = [obj.BlockHeight obj.BlockWidth];
                    tileBuf = repmat(obj.FillVal, tiledims);
                    th = internal.images.TileHandle(obj, blockId, obj.ReadOnly, tileBuf);
                    obj.TileCache.setTile(th);

                    % Set tile written to be true
                    obj.IsTileWritten(blockId) = true;
                end
            end
            
        end % getTileHandle
        
        
    end % private methods
    
end % MultiBandTiffAdapter


%-----------------------------------------------------------------
function [bits_per_sample sample_format] = datatype2Tags(datatype)

if strcmp(datatype, 'logical')
    bits_per_sample = 1;
    sample_format   = Tiff.SampleFormat.UInt;
elseif strcmp(datatype, 'int8')
    bits_per_sample = 8;
    sample_format   = Tiff.SampleFormat.Int;
elseif strcmp(datatype, 'uint8')
    bits_per_sample = 8;
    sample_format   = Tiff.SampleFormat.UInt;
elseif strcmp(datatype, 'int16')
    bits_per_sample = 16;
    sample_format   = Tiff.SampleFormat.Int;
elseif strcmp(datatype, 'uint16')
    bits_per_sample = 16;
    sample_format   = Tiff.SampleFormat.UInt;
elseif strcmp(datatype, 'int32')
    bits_per_sample = 32;
    sample_format   = Tiff.SampleFormat.Int;
elseif strcmp(datatype, 'uint32')
    bits_per_sample = 32;
    sample_format   = Tiff.SampleFormat.UInt;
elseif strcmp(datatype, 'single')
    bits_per_sample = 32;
    sample_format   = Tiff.SampleFormat.IEEEFP;
elseif strcmp(datatype, 'double')
    bits_per_sample = 64;
    sample_format   = Tiff.SampleFormat.IEEEFP;
else
    error(message('images:MultiBandTiffAdapter:unsupportedDatatype'));
end

end % datatype2Tags


%---------------------------------------------------------------
function datatype = tags2Datatype(bits_per_sample,sample_format)

if isequal(bits_per_sample,1)      && isequal(sample_format,Tiff.SampleFormat.UInt)
    datatype = 'logical';
elseif isequal(bits_per_sample,8)  && isequal(sample_format,Tiff.SampleFormat.Int)
    datatype = 'int8';
elseif isequal(bits_per_sample,8)  && isequal(sample_format,Tiff.SampleFormat.UInt)
    datatype = 'uint8';
elseif isequal(bits_per_sample,16) && isequal(sample_format,Tiff.SampleFormat.Int)
    datatype = 'int16';
elseif isequal(bits_per_sample,16) && isequal(sample_format,Tiff.SampleFormat.UInt)
    datatype = 'uint16';
elseif isequal(bits_per_sample,32) && isequal(sample_format,Tiff.SampleFormat.Int)
    datatype = 'int32';
elseif isequal(bits_per_sample,32) && isequal(sample_format,Tiff.SampleFormat.UInt)
    datatype = 'uint32';
elseif isequal(bits_per_sample,32) && isequal(sample_format,Tiff.SampleFormat.IEEEFP)
    datatype = 'single';
elseif isequal(bits_per_sample,64) && isequal(sample_format,Tiff.SampleFormat.IEEEFP)
    datatype = 'double';
else
    error(message('images:MultiBandTiffAdapter:unsupportedDatatype'));
end

end % tags2Datatype



