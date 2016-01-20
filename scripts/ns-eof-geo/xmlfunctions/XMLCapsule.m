classdef XMLCapsule < handle
    
    properties
        folder
        file
        root
    end
    
    methods
        
        function init(self, folder, file)
            self.folder = folder;
            self.file   = file;
            % laod xml-file
            self.root   = XMLParser.getRoot(sprintf('%s/%s', folder,file));
            XMLParser.clean(self.root);
        end
        
        function res = getRoot(self)
           res = self.root; 
        end
        
        function res = getChildren(self,sts)
            % find node
            node = self.root;

            for st=sts
                node = XMLParser.getChildren(node,char(st));
                node = node{1};
            end
            
            % get children
           res = XMLParser.getValidChildren(node); 
        end
        
        function node = getNode(self,sts)
            % find node
            node = self.root;

            for st=sts
                node = XMLParser.getChildren(node,char(st));
                node = node{1};
            end
        end
        
        function value = getValue(self,sts)
            % find node
            node = self.root;

            for st=sts
                node = XMLParser.getChildren(node,char(st));
                node = node{1};
            end
            
            % get value
            value = XMLParser.getValue(node);
        end
        
        function value = getAttribute(self,sts,attr)
            % find node
            node = self.root;

            for st=sts
                node = XMLParser.getChildren(node,char(st));
                node = node{1};
            end
            
            % get value
            value = XMLParser.getAttribute(node,attr);
        end
        
        function setValue(self,sts,value)
            % find node
            node = self.root;

            for st=sts
                node = XMLParser.getChildren(node,char(st));
                node = node{1};
            end
            
            XMLParser.setAttribute(node,'value',value);
        end
        
        function f = createFile(self, format)
            % new file name
            f = strcat(self.folder,'\\', sprintf(format,strrep(self.file,'.xml','')));
        end
        
        function save(self)
             sprintf('%s/%s', self.folder,self.file)
             XMLParser.setRoot(sprintf('%s/%s', self.folder,self.file),self.root);
        end
        
    end
    
end

