function makefigures(tldir, outputdir)
datasets = listdirs(tldir);
mkdir(outputdir)
for n = 1:size(datasets, 2)
    d = datasets{n};
    f = [d, '/r.mat'];
    if exist(f)
        disp(['Loading ', f])
        clear Dss Sss phi0 Tss DssMean
        load(f)
        clf
        plot_datapoints_polar(Dss)
        hold on
        plot_landmarks_polar(Sss)
        plot_outline_polar(Tss)
        
        % Make a nice name for the file
        sname = d;
        sname = sname((length(tldir)+2):end);
        sname(sname=='/') = '_';
        sname = [outputdir, '/', sname, '_polar.pdf'];
        
        print('-dpdf', sname);
    end
end
