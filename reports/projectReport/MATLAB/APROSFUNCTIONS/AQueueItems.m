classdef AOPCItems < handle
    %AOPCITEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        AW
        TIT
        M
    end
    
    methods
        function self = AOPCItems(AW,file)
            self.AW = AW;

            fid = fopen(file); 
            TITL = strsplit(fgetl(fid),','); 
            TITL(length(TITL)) = [];
            fclose(fid);

            self.TIT = TITL;

            self.M = dlmread(file,',',1,0);

        end
        
        function setValue(self,v)
%             temp = self.M(:,1);
%             
%             if temp(1) >= v
%                 v = temp(1);
%             elseif temp(end) <= v
%                 v = temp(end);
%             else
%                 for i = 1:length(temp)-1
%                     if temp(i)==v
%                         v = temp(i);
%                     elseif temp(i)<=v && temp(i+1)>=v
%                         v = [v temp(i) temp(i+1)];
%                     end
%                 end
%             end
            
            
            for i=1:length(self.TIT)
                vq = interp1(self.M(:,1),self.M(:,i),v);
                self.AW.writei(sprintf('MODI %s %d',char(self.TIT(i)),vq));  
%                 self.AW.writei(sprintf('MODI %s %d',char(self.TIT(i)),self.M(v,i)));                
            end
        end

        function renameProperty(self,propo,propn)
            for i=1:length(self.TIT)
                self.TIT(i)
                if strcmp(self.TIT(i),propo)
                    self.TIT(i)={propn};
                    self.TIT(i)
                    break;
                end
            end
        end

        function deleteProperty(self,prop)

        end


    end
    
end

