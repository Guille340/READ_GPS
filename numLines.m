%  nlines = NUMLINES(filePath)
%
%  DESCRIPTION
%  Returns the number of lines in a file.
%
%  NUMLINES is useful to preallocate memory for a variable that is going
%  to store the information of each line found in the file. 
%  
%  Pre-allocating will reduce the time used to populate a variable that 
%  changes dynamically in each iteration of a loop, specially when the 
%  number of iterations is particularly large.
%
%  INPUT VARIABLES
%  - filePath: absolute path of the file
%
%  OUTPUT VARIABLES
%  - nlines: number of lines in the file
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  FUNCTION CALLS
%  nlines = numLines(filePath)

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  14 Nov 2014

function nlines = numLines(filePath)

% Open File
fid = fopen(filePath, 'rb');

% Get File Size
fseek(fid, 0, 'eof');
fileSize = ftell(fid);
frewind(fid);

% Read the File and Count Lines
data = fread(fid, fileSize, 'uint8'); % read the whole file.
nlines = sum(data == 10); % count number of line-feeds and increase by one.
fclose(fid);

% % ALTERNATIVE METHOD (a bit slower)
% fid = fopen(filePath);
% allText = textscan(fid,'%s','delimiter','\n');
% nlines = length(allText{1});
% fclose(fid)