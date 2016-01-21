function [ M ] = SignalToMatrix( file )

    fid=fopen(file);
    inenum = 1;
    C = textscan(fid, '%d', 1, 'delimiter', '\n', 'headerlines', inenum-1);
    c = cell2mat(C);
    M = dlmread(file,'\t',c+1,0);
%     M = M(end,2:end);

end

