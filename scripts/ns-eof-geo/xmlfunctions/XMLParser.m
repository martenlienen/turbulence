classdef XMLParser

    methods(Static)
        
        function root = getRoot(file)
            % return root of xml-file
            xDoc = xmlread(file);
            root = xDoc.getDocumentElement;
        end
        
        function setRoot(file,root)
            xmlwrite(file,root);
        end
        
        function nodes = getAllChildren(node)
            childNodes = node.getChildNodes();
            ls = {};

            for i=0:childNodes.getLength()-1 
                item = childNodes.item(i);
                ls(end+1) = item;             
            end
            nodes = ls;
        end
        
        function nodes = getChildren(node, name)
            childNodes = node.getChildNodes();
            childNodes = node.getElementsByTagName(name);
            ls = {};

            for i=0:childNodes.getLength()-1 
                item = childNodes.item(i);
                if item.getNodeName() == name
                    ls(end+1) = item;
                end
                
            end
            nodes = ls;
            
%             nodes = node.getElementsByTagName(name);
        end
        
        function clean(node)
            childs = XMLParser.getAllChildren(node);
            
            if isempty(childs)
               return 
            end
            
            for i = 1:length(childs)
               child = childs{i};
               if ~strcmp(child.getNodeName(),'#text')
                  XMLParser.clean(child);
               else
                   node.removeChild(child);
               end
               
            end
        end
        
        function deleteChildren(node)
            childs = XMLParser.getAllChildren(node);
            for i = 1:length(childs)
               child = childs{i};
               node.removeChild(child);
            end
        end
        
        function nn = appdendChild(node,str)
            nn = node.getOwnerDocument().createElement(str);
            node.appendChild(nn);
        end
        
        function nodes = getValidChildren(node)
            childNodes = node.getChildNodes();
            ls = {}; name = '#text';
            
            for i=0:childNodes.getLength()-1 
                item = childNodes.item(i);
                
                if ~strcmp(item.getNodeName(),name)
                    ls(end+1) = item;
                end
                
            end
            nodes = ls;
        end
        
        function attr = getName(node)
            attr = node.getNodeName();
        end
        
        function attr = getValue(node)
            attr = char(node.getNodeValue());
            attr
            char(node.getNodeName())
            char(node.toString)
            node.getTextContent
            char(node.getTextContent)
        end
        
        function attr = getAttribute(node, name)
            attr = char(node.getAttributes.getNamedItem(name).getNodeValue());
        end
        
        function appendAttribute(node, name,attr)
            node.setAttribute(name,attr);
        end
        
        function setAttribute(node, name,attr)
            node.getAttributes.getNamedItem(name).setNodeValue(attr);
        end
        
        
    end
    
end

