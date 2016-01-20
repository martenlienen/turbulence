function [ text ] = getline( file,line )

    fid = fopen(file,'rt');
    s   = textscan(fid,'%s','delimiter','\n');
    textt = s{1}(line);
    text = textt{1,1};
end

