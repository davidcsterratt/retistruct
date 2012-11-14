function dirs = listdirs(tldir) 
files = dir(tldir);
dirs = {};
for n = 1:size(files, 1)
    f = files(n);
    if f.isdir
        if ~strcmp(f.name, '.') & ~strcmp(f.name, '..')
            dirs = {dirs{:}, [tldir '/' f.name]};
        end
    end
end
newdirs = {};
for n = 1:size(dirs, 2) 
    subdirs = listdirs(dirs{n});
    newdirs = {newdirs{:}, subdirs{:}};
end
dirs = {dirs{:}, newdirs{:}};

