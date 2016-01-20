function [ id ] = findline( file,text )

    fid = fopen(file,'rt');
    s   = textscan(fid,'%s','delimiter','\n');
    id = find(strcmp(s{1},text),1,'first');
    
end

