function rslt = funcGetFilePath(rootpath, contentFilter)

%% get path of imu data and video data


%% get paths of IMU and Video data
filelist = dir(rootpath);

filename = {filelist.name}.';


%% Content Filter

%contentFilter = "20211117_1413_43";
filename = filename(contains(filename, contentFilter));

%% check number of IMU files and Video files

if ((size(filename,1)) < 1)
    %if number of files is equal
    disp("No Files Found.");
    rslt = 0;
    return
end


rslt = filename;
end
