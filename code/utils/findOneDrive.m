function pathOut = findOneDrive()

% Check if we're on Windows
if ispc
    % Try to get OneDrive path using PowerShell
    [status, result] = system(['powershell -command ' ...
        '"$env:OneDriveCommercial"']);
    if status == 0 && ~isempty(strtrim(result))
        pathOut = strtrim(result);
        pathOut = strrep(pathOut, '"', '');
        return;
    end

    % If that fails, try personal OneDrive
    [status, result] = system('powershell -command "$env:OneDrive"');
    if status == 0 && ~isempty(strtrim(result))
        pathOut = strtrim(result);
        pathOut = strrep(pathOut, '"', '');
        return;
    end
else
    % For macOS and Linux
    home = getenv('HOME');
    switch home
        case '/Users/bellapatel'
            pathOut = ['/Users/bellapatel/Library/CloudStorage/' ...
                'OneDrive-UniversityofPittsburgh/'];
        otherwise
            homeDirList = mydirsub(home);
            pathOut = [home filesep ...
                homeDirList{contains(homeDirList, 'OneDrive')}];
    end

end
end

function dirlist = mydirsub(string)

% dirlist = mydirsub(string)

% use dir to get the names
fc = dir(string);

% return a cell array of names that do not start with '.'
dirlist = {fc(cellfun(@(x)~any(x == 1), strfind({fc.name}, '.'))).name};
end