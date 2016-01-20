%% Function for creating a mesh (only uniform possible)
function [ x,y ] = createmesh( conf )


    lengthx = conf.lx;
    lengthy = conf.ly;

    cellsx  = conf.sx;
    cellsy  = conf.sy;

    if strcmp(conf.mesh,'uniform')
        x = linspace(0,lengthx,cellsx); 
        y = linspace(0,lengthy,cellsy); 
    else
        
    end

end

