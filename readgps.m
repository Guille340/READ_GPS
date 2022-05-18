%  GpsData = READGPS(fpath,varargin)
%
%  DESCRIPTION
%  Reads the navigation information from one or more GPS files and saves
%  it in output structure GPSDATA. The fields in the output structure
%  vary depending on the file format and type of GPS sentence. Two file
%  formats are supported: SeicheSSV (*.gpstext) and PamGuard (*.csv). For
%  the first type, the GPS identifier can be set ($GPRMC by default). The
%  output only includes data from valid GPS sentences. The validity of
%  the sentences is assessed by the NMEA checksum (GPGGA,GPGSA) or by
%  simple inspection of the data "status" field (GPRMC).
%
%  INPUT VARIABLES
%  - fpath: character string of cell vector of character strings specifying
%    the full path of the GPS file(s).
%  - gpid (varargin{1}) [string]: type of GPS sentence. Three options:
%    ¬ 'GPRMC': recommended minimum specific GPS and transit data
%    ¬ 'GPGGA': essential fix data (3D location and accuracy)
%    ¬ 'GPGSA': dilution of precision and active satellites
%
%  OUTPUT VARIABLES
%  - GpsData [structure]: navigation data. The fields in the struct change
%    depending on the sentence type:
%
%      SeicheSSV       SeicheSSV      SeicheSSV       PamGuard
%       'GPRMC'         'GPGGA'        'GPGSA'
%       ---------------------------------------------------------------
%        pctick         pctick          pctick        pctick
%        utctick        utctick         fmode         utctick
%        lat            lat             pdop          lat
%        lon            lon             hdop          lon
%        cmg            fixq            vdop          sog
%        sog            nsat                          hea
%                       hdop
%                       alt
%                       hog
%
%    * GPSDATA FIELDS
%      ¬ pctick: PC clock timestamp (referred to '01-jan-0001 00:00:00')[s]
%      ¬ utctick: UTC GPS timestamp (referred to '01-jan-0001 00:00:00')[s]
%      ¬ lat: latitude [deg]
%      ¬ lon: longitude [deg]
%      ¬ cmg: course made good [deg]
%      ¬ sog: speed over ground [kts]
%      ¬ hea: heading [deg, re. North]
%      ¬ fixq: fix quality (0 = invalid, 1 = GPS fix, 2 = DGPS fix)
%      ¬ nsat: number of satellites in view
%      ¬ alt: altitude above mean sea level [m]
%      ¬ hog: height of geoid above WGS84 ellipsoid [m]
%      ¬ hdop: horizontal dilution of precision
%      ¬ vdop: vertical dilution of precision
%      ¬ pdop: position dilution of precision
%      ¬ fmode: fix mode (1 = fix not available, 2 = 2D, 3 = 3D)
%
%  INTERNALLY CALLED FUNCTIONS
%  - nmeaChecksum
%  - numLines
%
%  CONSIDERATIONS & LIMITATIONS
%  - A checksum is applied in GPGGA and GPGSA sentences, but not in GPRMC.
%    In GPRMC the "navigation receiver warning" is used to determine
%    whether a sentence is valid or not (A = OK, V = warning). This is a
%    much faster method than the NMEA checksum.
%  - All files must be of the same type (SeicheSSV, *.gpstext; PamGuard,
%    *.csv).
%
%  FUNCTION CALLS
%  1. GpsData = readgps (fpath)
%  2. GpsData = readgps (fpath,gpid)
%
%  See also NMEACHECKSUM, NUMLINES

%  VERSION 2.1
%  Date: 16 May 2019
%  Author: Guillermo Jimenez Arranz
%  - utctick assigned to "GpsDate" column for .csv files (PamGuard),
%    instead of "UTC" column.
%  - pctick assigned to "UTC" column for .csv files (PamGuard), instead
%    of "PCTime" column.
%
%  VERSION 2.0
%  Date: 14 Apr 2018
%  Author: Guillermo Jimenez Arranz
%  - Added support for PamGuard GPS format (*.csv)
%  - Multiple GPS files can be read.
%
%  VERSION 1.1
%  Date: 20 Apr 2016
%  Author: Guillermo Jimenez Arranz
%  - The output argument is a struct, instead a variable output argument.
%    This leads to simpler code and more intuitive use of the function and
%    the output data.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  12 Feb 2016

function GpsData = readgps(fpath,varargin)

% Number of Selected Files
if iscell(fpath)
    nFiles = numel(fpath);
end

if ischar(fpath)
    fpath = {fpath};
    nFiles = 1;
end

% Error Control: File Extension
fext = cell(1,nFiles);
for k = 1:nFiles
    [~,~,fext{k}] = fileparts(fpath{k});
end

iscsv = strcmp(fext,'.csv');
isgpstext = strcmp(fext,'.gpstext');
allcsv = all(iscsv);
allgpstext = all(isgpstext);

if any(~iscsv & ~isgpstext)
    error('One or more selected files have an unrecognised or unsupported format')
end

if ~allcsv && ~allgpstext
    error('All selected files must have the same format')
end

% File Format
if allcsv
    fext = '.csv';
else
    fext = '.gpstext';
end

% Extract Positioning Data
switch fext
    case '.gpstext' % Read GPS Data (SeicheSSV)
        % ERROR CONTROL: Number of Input Arguments
        switch nargin
            case 0
                error('Not enough input arguments')
            case 1
                gpid = 'GPRMC';
            case 2
                gpid =  varargin{1};
            otherwise
                error('Too many input arguments')
        end

        % READ FILE
        switch gpid
            case 'GPRMC'

                % Initialise GPRMC structure
                GpsData = struct('pctick',[],'utctick',[],'lat',[],'lon',[],'sog',[],'cmg',[]);

                % Read Positioning Data from Files
                i1 = 1;
                for k = 1:nFiles
                    % Open File
                    fid = fopen(fpath{k},'r'); % open GPS file in 'read' mode
                    nlin = numLines(fpath{k}); % number of sentences in the GPS file

                    % Initialise Variables
                    pctick = cell(nlin,1); % pc tick (ref. to '00-jan-0000 00:00:00') [s]
                    utctick = cell(nlin,1); % utc tick (ref. to '00-jan-0000 00:00:00') [s]
                    lat = NaN(nlin,1); % latitude [deg]
                    lon = NaN(nlin,1); % longitude [deg]
                    sog = NaN(nlin,1); % speed over ground [kts]
                    cmg = NaN(nlin,1); % course made good [deg]
                    m = 0; % counter for valid sentences

                    % Read File
                    for n = 1:nlin
                        lin = fgets(fid); % get line
                        data = textscan(lin,'%s','Delimiter',',');
                        data = data{1}; % cell array with the data from the GPS sentence
                        nElements = length(data);
                        
                        cs = false;
                        if nElements >= 4
                            cs = strcmp(data{4},'A'); % true for valid checksum (FASTER)
                        end
                        
                        gpid0 = '';
                        if nElements >= 2
                            gpid0 = data{2}; % GPS identifier ($GP...)
                            gpid0 = gpid0(2:end); % GPS identifier (GP...)                            
                        end

                        if strcmp(gpid0,gpid) && cs
                            m = m + 1; % update counter
                            
                            % Read Times
                            pctick{m} = data{1}; % PC timestamp (char)
                            utctick{m} = [data{11} data{3}]; % GPS timestamp (char)
                            
                            % Read Latitude
                            latstr = data{5}; % latitude (char)
                            lat_dm = sscanf(latstr,'%2f%7f'); % latitude [deg min]
                            
                            % Sign of Latitude (+1 North, -1 South)
                            if strcmp(data{6},'S'), ns = -1; else, ns = 1; end 
                            if ~isempty(lat_dm)
                                lat(m) = ns*(lat_dm(1) + lat_dm(2)/60); 
                            else
                                lat(m) = NaN; 
                            end 
                            
                            % Read Longitude
                            lonstr = data{7}; % longitude (char)
                            lon_dm = sscanf(lonstr,'%3f%7f'); % longitude [deg min]
                            
                            % Sign of Longitude (+1 North, -1 South)
                            if strcmp(data{8},'W'), ew = -1; else, ew = 1; 
                            end
                            if ~isempty(lon_dm)
                                lon(m) = ew*(lon_dm(1) + lon_dm(2)/60);
                            else 
                                lon(m) = NaN; 
                            end 
                            
                            % Read Speed Over Ground (SOG)
                            sog0 = sscanf(data{9},'%f'); % sog [kts]
                            if ~isempty(sog0)
                                sog(m) = sog0; 
                            end 
                            
                            % Read Course Made Good (CMG)
                            cmg0 = sscanf(data{10},'%f'); % cmg [deg]
                            if ~isempty(cmg0)
                                cmg(m) = cmg0; 
                            end 
                        end
                    end
                    fclose(fid); % close GPS file

                    if m > 1
                        pctick = datenum(pctick(1:m),...
                            'yyyymmdd_HHMMSS_FFF')*86400; % PC tick time vector (ref. '01-jan-0001 00:00:00') [s]
                        utctick = datenum(utctick(1:m),...
                            'ddmmyyHHMMSS.FFF')*86400; % GPS tick time vector (ref. '01-jan-0001 00:00:00') [s]
                    else
                        pctick = [];
                        utctick = [];
                    end
                    lat = lat(1:m); % latitude vector [deg]
                    lon = lon(1:m); % longitude vector [deg]
                    sog = sog(1:m); % speed over ground vector [kts]
                    cmg = cmg(1:m); % course made good vector [deg]

                    % Generate Output Structure
                    i2 = i1 + m - 1;
                    GpsData.pctick(i1:i2,1) = pctick;
                    GpsData.utctick(i1:i2,1) = utctick;
                    GpsData.lat(i1:i2,1) = lat;
                    GpsData.lon(i1:i2,1) = lon;
                    GpsData.sog(i1:i2,1) = sog;
                    GpsData.cmg(i1:i2,1) = cmg;

                    i1 = i2 + 1;
                end

                % Sort Data by Time
                [~,isort] = sort(GpsData.utctick);
                GpsData.utctick = GpsData.utctick(isort);
                GpsData.pctick = GpsData.pctick(isort);
                GpsData.lat = GpsData.lat(isort);
                GpsData.lon = GpsData.lon(isort);
                GpsData.sog = GpsData.sog(isort);
                GpsData.cmg = GpsData.cmg(isort);

            case 'GPGGA'

                % Initialise GPGGA structure
                GpsData = struct('pctick',[],'utctick',[],'lat',[],'lon',[],...
                    'fixq',[],'nsat',[],'hdop',[],'alt',[],'hog',[]) ;

                % Read Positioning Data from Files
                i1 = 1;
                for k = 1:nFiles
                    % Open File
                    fid = fopen(fpath{k},'r'); % open GPS file in 'read' mode
                    nlin = numLines(fpath{k}); % number of sentences in the GPS file

                    % Initialise Variables
                    pctick = cell(nlin,1); % pc tick (ref. to '01-jan-0001 00:00:00') [s]
                    utctick = cell(nlin,1); % utc tick (ref. to '00:00:00' of current day) [s]
                    lat = NaN(nlin,1); % latitude
                    lon = NaN(nlin,1); % longitude
                    fixq = NaN(nlin,1); % fix quality (0 = invalid, 1 = GPS fix, 2 = DGPS fix)
                    nsat = NaN(nlin,1); % number of satellites
                    hdop = NaN(nlin,1); % Horizontal Dilution Of Precision
                    alt = NaN(nlin,1); % altitude
                    hog = NaN(nlin,1); % height of geoid above WGS84 ellipsoid
                    m = 0; % counter for valid sentences

                    % Read File
                    for n = 1:nlin
                        lin = fgets(fid); % get line
                        data = textscan(lin,'%s','Delimiter',',');
                        data = data{1}; % cell array with the data from the GPS sentence
                        
                        gpid0 = '';
                        nElements = length(data);
                        if nElements >= 2
                            gpid0 = data{2}; % GPS identifier ($GP...)
                            gpid0 = gpid0(2:end); % GPS identifier (GP...)                            
                        end

                        if strcmp(gpid0,gpid)
                            if nmeaChecksum(lin) % true for valid checksum
                                m = m + 1; % update counter
                                
                                % Read Times
                                pctick{m} = data{1}; % PC timestamp (char)
                                utctick{m} = ['00010000' data{3}]; % time referred to 00:00 am of current day
                                
                                % Read Latitude
                                latstr = data{4}; % (char)
                                lat_dm = sscanf(latstr,'%2f%7f'); % [deg min]
                                
                                % Sign of Latitude (+1 North, -1 South)
                                if strcmp(data{5},'S')
                                    ns = -1; 
                                else
                                    ns = 1; 
                                end 
                                if ~isempty(lat_dm)
                                    lat(m) = ns*(lat_dm(1) + lat_dm(2)/60); 
                                end 
                                
                                % Read Longitude
                                lonstr = data{6}; % longitude (char)
                                lon_dm = sscanf(lonstr,'%3f%7f'); % [deg min]
                                
                                % Sign of Longitude (+1 East, -1 West)
                                if strcmp(data{7},'W')
                                    ew = -1; 
                                else
                                    ew = 1; 
                                end 
                                if ~isempty(lon_dm)
                                    lon(m) = ew*(lon_dm(1) + lon_dm(2)/60); % [deg]
                                end 
                                
                                % Read Fix Quality (0 = invalid, 1 = GPS fix, 
                                % 2 = DGPS fix)
                                fixq0 = sscanf(data{8},'%f'); % fix quality 
                                if ~isempty(fixq0), fixq(m) = fixq0; end
                                
                                % Read Number of Satellites
                                nsat0 = sscanf(data{9},'%f'); 
                                if ~isempty(nsat0), nsat(m) = nsat0; end 
                                
                                % Read Horizontal Dilution of Precision
                                hdop0 = sscanf(data{10},'%f'); 
                                if ~isempty(hdop0), hdop(m) = hdop0; end 
                                
                                % Read Altitude (above MSL)
                                alt0 = sscanf(data{11},'%f'); 
                                if ~isempty(alt0), alt(m) = alt0; end 
                                
                                % Read Height of Geoid (above WGS84 ellipsoid)
                                hog0 = sscanf(data{12},'%f'); 
                                if ~isempty(hog0), hog(m) = hog0; end
                            end
                        end
                    end
                    fclose(fid);

                    if m > 1
                        % PC tick time vector (ref. '01-jan-0001 00:00:00') [s]
                        pctick = datenum(pctick(1:m),...
                            'yyyymmdd_HHMMSS_FFF')*86400; 
                        % PC tick time vector (ref. 00:00:00 of current day) [s]
                        utctick = datenum(utctick(1:m),...
                            'ddmmyyyyHHMMSS.FFF')*86400; 
                    else
                        pctick = [];
                        utctick = [];
                    end
                    lat = lat(1:m); % latitude vector [deg]
                    lon = lon(1:m); % longitude vector [deg]
                    fixq = fixq(1:m); % fix quality vector
                    nsat = nsat(1:m); % number of satellites vector
                    hdop = hdop(1:m); % horizontal dilution of precision vector
                    alt = alt(1:m); % altitude vector [m]
                    hog = hog(1:m); % height over geoid vector [m]

                    % Generate Output Structure
                    i2 = i1 + m - 1;
                    GpsData.pctick(i1:i2,1) = pctick;
                    GpsData.utctick(i1:i2,1) = utctick;
                    GpsData.lat(i1:i2,1) = lat;
                    GpsData.lon(i1:i2,1) = lon;
                    GpsData.fixq(i1:i2,1) = fixq;
                    GpsData.nsat(i1:i2,1) = nsat;
                    GpsData.hdop(i1:i2,1) = hdop;
                    GpsData.alt(i1:i2,1) = alt;
                    GpsData.hog(i1:i2,1) = hog;

                    i1 = i2 + 1;
                end

                % Sort Data by Time
                [~,isort] = sort(GpsData.utctick);
                GpsData.utctick = GpsData.utctick(isort);
                GpsData.pctick = GpsData.pctick(isort);
                GpsData.lat = GpsData.lat(isort);
                GpsData.lon = GpsData.lon(isort);
                GpsData.fixq = GpsData.fixq(isort);
                GpsData.nsat = GpsData.nsat(isort);
                GpsData.hdop = GpsData.hdop(isort);
                GpsData.alt = GpsData.alt(isort);
                GpsData.hog = GpsData.hog(isort);

                case 'GPGSA'

                    % Initialise GPGSA structure
                    GpsData = struct('pctick',[],'fmode',[],'pdop',[],...
                        'hdop',[],'vdop',[]);

                    % Read Positioning Data from Files
                    i1 = 1;
                    for k = 1:nFiles
                        % Open File
                        fid = fopen(fpath{k},'r'); % open GPS file in 'read' mode
                        nlin = numLines(fpath{k}); % number of sentences in the GPS file

                        % Initialise Variables
                        pctick = cell(nlin,1); % pc tick (ref. to '01-jan-0001 00:00:00') [s]
                        fmode = NaN(nlin,1); % fix mode (0 = fix not available, 1 = 2D, 2 = 3D)
                        pdop = NaN(nlin,1); % Position Dilution Of Precision (combination of hdop and vdop)
                        hdop = NaN(nlin,1); % Horizontal Dilution Of Precision
                        vdop = NaN(nlin,1); % Vertical Dilution Of Precision
                        m = 0; % counter for valid sentences

                        % Read File
                        for n =1:nlin
                            lin=fgets(fid); % get line
                            data = textscan(lin,'%s','Delimiter',',');
                            data = data{1}; % cell array data from GPS sentence
                            
                            gpid0 = '';
                            nElements = length(data);
                            if nElements >= 2
                                gpid0 = data{2}; % GPS identifier ($GP...)
                                gpid0 = gpid0(2:end); % GPS identifier (GP...)                            
                            end

                            if strcmp(gpid0,gpid)
                                if nmeaChecksum(lin) % true for valid checksum
                                    m = m + 1; % update counter
                                    
                                    % Read Time
                                    pctick{m} = data{1}; % PC timestamp (char)
                                    
                                    % Read Fix Mode (0 = fix not available, 1 = 2D, 2 = 3D)
                                    fmode0 = sscanf(data{4},'%d'); 
                                    if ~isempty(fmode), fmode(m) = fmode0; end 
                                    
                                    % Read Position Dilution Of Precision
                                    pdop0 = sscanf(data{17},'%f');
                                    if ~isempty(pdop0), pdop(m) = pdop0; end 
                                    
                                    % Read Horizontal Dilution Of Precision
                                    hdop0 = sscanf(data{18},'%f'); 
                                    if ~isempty(hdop0), hdop(m) = hdop0; end 
                                    
                                    % Read Vertical Dilution Of Precision
                                    vdop0 = sscanf(data{19},'%f*%*s');
                                    if ~isempty(vdop0), vdop(m) = vdop0; end 
                                end

                            end
                        end
                        fclose(fid);

                        if m > 1
                            % PC tick time vector (ref. '01-jan-0001 00:00:00') [s]
                            pctick = datenum(pctick(1:m),...
                                'yyyymmdd_HHMMSS_FFF')*86400; 
                        else
                            pctick = [];
                        end
                        fmode = fmode(1:m); % fix mode vector
                        pdop = pdop(1:m); % Position Dilution of Precision vector
                        hdop = hdop(1:m); % Horizontal Dilution of Precision vector
                        vdop = vdop(1:m); % Vertical Dilution of Precision vector

                        % Generate Output Structure
                        i2 = i1 + m - 1;
                        GpsData.pctick(i1:i2,1) = pctick;
                        GpsData.fmode(i1:i2,1) = fmode;
                        GpsData.pdop(i1:i2,1) = pdop;
                        GpsData.hdop(i1:i2,1) = hdop;
                        GpsData.vdop(i1:i2,1) = vdop;

                        i1 = i2 + 1;
                    end

                    % Sort Data by Time
                    [~,isort] = sort(GpsData.pctick);
                    GpsData.pctick = GpsData.pctick(isort);
                    GpsData.fmode = GpsData.fmode(isort);
                    GpsData.pdop = GpsData.pdop(isort);
                    GpsData.hdop = GpsData.hdop(isort);
                    GpsData.vdop = GpsData.vdop(isort);

            otherwise
                error('Invalid GPID')
        end

    case '.csv' % Read GPS Data (PamGuard)

        % Error Control: Number of Input Arguments
        if nargin < 1, error('Not enough input arguments'); end
        if nargin > 1, error('Too many input arguments'); end

        % Initialise GPS structure
        GpsData = struct('pctick',[],'utctick',[],'lat',[],'lon',[],...
            'sog',[],'hea',[]); % initialise GPS structure

        % Read Positioning Data from Files
        i1 = 1;
        for k = 1:nFiles

            % Open .csv File
            fid = fopen(fpath{k});
            ncol = length(strfind(fgets(fid),',')) + 1; % no. columns in CSV
            frewind(fid)
            data = textscan(fid,'%s','delimiter',',');
            fclose(fid);

            % Convert Cell Vector to Char Array
            nrow = length(data{1})/ncol;
            datatxt = reshape(data{1},ncol,nrow)';

            % Find Valid Positioning Sentences
            gpserror = str2double(cell2mat(datatxt(2:nrow,17)));
            datastatus = deblank(cell2mat(datatxt(2:nrow,18)));
            ivalid = find(~gpserror & datastatus=='A')+1;

            % Assign Data to Variables
            pctick = datenum(datatxt(ivalid,2),'yyyy-mm-dd HH:MM:SS.FFF')*86400; % PC tick time vector (ref. '01-jan-0001 00:00:00') [s]
            utctick = datenum(datatxt(ivalid,6),'yyyy-mm-dd HH:MM:SS.FFF')*86400; % GPS tick time vector (ref. '01-jan-0001 00:00:00') [s]
            lat = str2double(datatxt(ivalid,8)); % latitude vector [deg]
            lon = str2double(datatxt(ivalid,9)); % longitude vector [deg]
            sog = str2double(datatxt(ivalid,10)); % speed over ground vector [kts]
            hea = str2double(datatxt(ivalid,12)); % true heading vector [deg]

            % Generate Output Structure
            i2 = i1 + length(ivalid) - 1;
            GpsData.pctick(i1:i2,1) = pctick;
            GpsData.utctick(i1:i2,1) = utctick;
            GpsData.lat(i1:i2,1) = lat;
            GpsData.lon(i1:i2,1) = lon;
            GpsData.sog(i1:i2,1) = sog;
            GpsData.hea(i1:i2,1) = hea;

            i1 = i2 + 1;
        end

        % Sort Data by Time
        [~,isort] = sort(GpsData.utctick);
        GpsData.utctick = GpsData.utctick(isort);
        GpsData.pctick = GpsData.pctick(isort);
        GpsData.lat = GpsData.lat(isort);
        GpsData.lon = GpsData.lon(isort);
        GpsData.sog = GpsData.sog(isort);
        GpsData.hea = GpsData.hea(isort);

        % Remove Duplicated Ticks
        [~,iuni,~] = unique(GpsData.utctick,'stable');
        GpsData.utctick = GpsData.utctick(iuni);
        GpsData.pctick = GpsData.pctick(iuni);
        GpsData.lat = GpsData.lat(iuni);
        GpsData.lon = GpsData.lon(iuni);
        GpsData.sog = GpsData.sog(iuni);
        GpsData.hea = GpsData.hea(iuni);
end

