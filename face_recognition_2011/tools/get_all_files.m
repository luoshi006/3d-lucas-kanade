function results = get_all_files(dirName, varargin)
%% Parse input
inp = inputParser;

inp.addRequired('dirName', @(x)ischar(x));
inp.addOptional('regex', '', @(x)ischar(x));

inp.parse(dirName, varargin{:});
arg = inp.Results;
regex = arg.regex;
clear('inp');

%% Get all file names in directory

dirData = dir(dirName);     
dirIndex = [dirData.isdir]; 
fileList = {dirData(~dirIndex).name}';

results = {};

for i=1:length(fileList)
    s = fileList(i);
    s = s{1};
    if (isempty(regex) || ~isempty(regexpi(s, regex)))
        results{end + 1} = s;
    end
end