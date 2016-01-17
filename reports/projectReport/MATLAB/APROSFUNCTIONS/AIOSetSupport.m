classdef AIOSetSupport  < handle
    %IOSUPPORT Object to manage IO of Apros
    %   Detailed explanation goes here
    
    properties
        IONAME
        AW
        AQS
        COUNTERI
        COUNTERO
    end
    
    methods
        function self = AIOSetSupport(AQS,AW)
            self.AQS = AQS;
            self.AW = AW;
        end
        
        function writei(self,text)
            self.AW.writei( text);
        end
        
        function init(self,IONAMET)
            self.COUNTERI = 100;
            self.COUNTERO = 0;
            
            self.IONAME = IONAMET;
            
            % Create IO_SET
            self.writei(sprintf('ADD IO_SET %s -',self.IONAME));
            self.writei('IO_STOP_TIME 2.0E9 -');
            self.writei('IO_TIME_STEP 0.2 -');
            self.writei('IO_DEVICE_PARAMS NULL -');

            f = @(x) self.writei(sprintf('IO_DB_NAMES(%d) NULL -',x));
            for i=1:20
               f(i); 
            end

            f = @(x) self.writei(sprintf('IO_CONVERSION(%d) NULL -',x));
            for i=1:20
               f(i); 
            end

            self.writei('IO_FILE_NAME  -');
            self.writei('IO_START_TIME 0.0 -');

            f = @(x) self.writei(sprintf('IO_EXT_NAMES(%d) NULL -',x));
            for i=1:20
               f(i); 
            end

            self.writei('IO_HEADER  -');
            self.writei('IO_TYPE output');
            self.writei(sprintf('INCLUDE COMP_ROOT %s',self.IONAME));          

            self.writei(sprintf('UPDATE IO_SET %s -',self.IONAME));
            self.writei('IO_ACL_DOUBLE_PREC 1');
        end
        
        function monitor(self,module, properties)
            
            for property=properties
                
                if(self.COUNTERI==100)
                   self.COUNTERI = 0;
                   self.COUNTERO = self.COUNTERO + 1;
                   
                    % Create DB_NAME
                    self.writei(sprintf('ADD DB_NAMES %s%d -',self.IONAME,self.COUNTERO));
                    f = @(x) self.writei(sprintf('DB_NAME(%d) NULL NULL 0 0 -',x));
                    for i=1:99
                        f(i); 
                    end
                    self.writei('DB_NAME(100) NULL NULL 0 0');
                    self.writei(sprintf('INCLUDE COMP_ROOT %s%d',self.IONAME,self.COUNTERO));
            

                    self.writei(sprintf('UPDATE IO_SET %s -',self.IONAME));
                    self.writei(sprintf('IO_DB_NAMES(%d) %s%d',self.COUNTERO,self.IONAME,self.COUNTERO));
                end
                
                self.COUNTERI = self.COUNTERI + 1;
                
                self.writei(sprintf('UPDATE DB_NAMES %s%d -',self.IONAME,self.COUNTERO));
                self.writei(sprintf('DB_NAME(%d) %s %s 1 1',self.COUNTERI,module,char(property)));
            end
            
        end
        
        function title(self,tit)
            self.writei(sprintf('UPDATE IO_SET %s -',self.IONAME));
            self.writei(sprintf('IO_FILE_NAME %s',tit));
        end
        
        function reset(self)
%             self.writei('MODI ECCO SIMULATION_CURRENT_TIME 0.0');
%             self.writei('SET SISECO 0');
        end
        
        function timestep(self,dt)
            self.writei(sprintf('MODI ECCO CURRENT_TIME_STEP %f',dt));
            self.writei(sprintf('MODI ECCO MAXIMUM_TIME_STEP %f',dt));
            self.writei(sprintf('UPDATE IO_SET %s -',self.IONAME));
            self.writei(sprintf('IO_TIME_STEP %f',dt));
        end
        
        function run(self,nr)
            self.AQS.do(nr);
%             self.writei(sprintf('do %d',nr));
        end
        
        function runAndWrite(self,nr)
            self.writei(sprintf('io open %s',self.IONAME));
            self.AQS.do(nr);
%             self.writei(sprintf('do %d',nr));
            self.writei(sprintf('io close %s',self.IONAME));
        end
        
        function finalize(self)
            for i=1:self.COUNTERO
                self.writei(sprintf('DELETE/ALL %s%d',self.IONAME,i));
            end
            
            self.writei(sprintf('DELETE/ALL %s',self.IONAME));
        end
        
        
    end
    
end

