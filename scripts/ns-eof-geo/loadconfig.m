%% function for reading the configuaration file and saving the needed settings in a struct
function [ s ] = loadconfig(folder, file )

    s ={};

    xml    = XMLCapsule();
    xml.init(folder,file);
    s.lx = str2double(xml.getAttribute({'geometry'},'lengthX'));
    s.ly = str2double(xml.getAttribute({'geometry'},'lengthY'));
    s.sx = str2double(xml.getAttribute({'geometry'},'sizeX'));
    s.sy = str2double(xml.getAttribute({'geometry'},'sizeY'));
    
    s.mesh = xml.getAttribute({'geometry','mesh'},'t');
end

