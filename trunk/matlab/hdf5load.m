function s = hdf5load(filename, varargin)
info = hdf5info(filename, varargin{:});
g = info.GroupHierarchy;
s = get_struct(filename, g, [], varargin{:});

function s = get_struct(filename, g, s, varargin) 
ds = g.Datasets;
gs = g.Groups;
for d=ds
    n = d.Name(2:end);
    data = hdf5read(filename, n, varargin{:});
    if strcmp(d.Datatype.Class, 'H5T_STRING')
        data = data.Data;
    end
    sn = ['s.' strrep(n, '/', '.')];
    eval([sn '= data;']);
end
for g=gs
    s = get_struct(filename, g, s, varargin{:});
end

