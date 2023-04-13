function [data] = readTrodesExtractedDataFile(filename, varargin)  
%fields = readTrodesExtractedDataFile(filename)
%filename-- a string containing the name of the binary file with extracted data
%
%data--the output structure, which contains information about the file and
%a data field called fields, a structure of length n, where n is the number of data types in the
%file. Each entry has the following fields:
%
%data.fields(n).name -- the name of the data
%data.fields(n).type -- the data type (i.e., uint32, uint16, int16, or double)
%data.fields(n).data -- the data, an n by 1 vector


%varargins
forcetype = ''; %forces datatype and ignores the header 'type' field value.. This is a temp hack instead of editing the header
if (~isempty(varargin))
    assign(varargin{:});
end


%fid = fopen(filename,'r','ieee-be');
fid = fopen(filename,'rb','ieee-le');
%fseek(fid,'eof');
%fileSize = ftell(fid);
%fseek(fid,'bof');
%headerReadSize = 10000;
%if (fileSize < 10000)
%    headerReadSize = fileSize-1;
%end


headerText = fread(fid,10000,'uint8');
headerText = char(headerText');
endHeaderLoc = strfind(headerText,'<End settings>');
                
startHeaderLoc = strfind(headerText,'<Start settings>')+17;
data = [];

%default fields in the binary data
fields(1).name = 'time';
fields(1).type = 'uint32';
fields(1).bytesPerItem = 4;
fields(1).columns = 1;
fields(1).data = [];
bytesPerPacket = 4;

clockRate = [];
voltageScaling = [];

if (~isempty(endHeaderLoc))
    headersize = endHeaderLoc+14;
    
    %Read in the header info to the output structure
    %headerLines = strread(char(headerText(startHeaderLoc:(endHeaderLoc-1))),'%s', 'delimiter', sprintf('\n'));
    %headerLines = textscan(char(headerText(startHeaderLoc:(endHeaderLoc-1))),'%s');
    headerLines = strsplit(char(headerText(startHeaderLoc:(endHeaderLoc-1))),sprintf('\n'));
    
    for l=1:length(headerLines)
       colonLoc = strfind(headerLines{l},':');
       if (~isempty(colonLoc))
           tmpFieldName = lower(regexprep(headerLines{l}(1:colonLoc(1)-1),'[^\w'']',''));
           tmpFieldVal = strtrim(headerLines{l}(colonLoc(1)+1:end));
           
           %make sure this is not the field structure description line
           if (~isequal(lower(tmpFieldName),'fields'))
               %if the value is convertable to a number, do it
               trynum = str2num(tmpFieldVal);
               if (~isempty(trynum))
                   data = setfield(data,tmpFieldName,trynum);
               else                 
                   data = setfield(data,tmpFieldName,tmpFieldVal);
               end
           end
       end
    end
    
    
%     voltageScaleLoc  = strfind(headerText,'Voltage_scaling:');
%     if (~isempty(voltageScaleLoc))
%         voltageScaling = str2num(char(strtok(headerText(voltageScaleLoc+17:end))));
%     end
    
%     clockRateLoc  = strfind(headerText,'Clock rate:');
%     if (~isempty(clockRateLoc))
%         clockRate = str2num(char(strtok(headerText(clockRateLoc+12:end))));
%     end
    
     %See if the file designates fields, and assign them
    
    fieldsLoc  = strfind(headerText,'Fields:');
    if (~isempty(fieldsLoc)) 
        fields = [];
        bytesPerPacket = 0;
        fseek(fid, fieldsLoc+7, -1);
        fieldString = fgetl(fid);
        remainder = fieldString;
        currentFieldNum = 1;
        while (~isempty(remainder))
            %each field is encapsulated in < >
            %inside, there is a fieldname, followed by a space, then the datatype
            %Example: <time uint32><voltage1 int16>
            %If a field is an array, the datatype has a multiplier in front like this:
            %<waveformCh1 40*int16>
            
            [token, remainder] = strtok(remainder,'<>');
            if (~isempty(token))
                [tmpField rem] = strtok(token);
                fields(currentFieldNum).name = tmpField;
                if ~isempty(forcetype)
                    fields(currentFieldNum).type = forcetype;
                else
                    fields(currentFieldNum).type = strtok(rem);
                end
                fields(currentFieldNum).columns = 1;
                fields(currentFieldNum).bytesPerItem = 0;
                multiplier = 1;
                multLoc = strfind(fields(currentFieldNum).type,'*');
                if ~isempty(multLoc)
                  tmpMultiplier = str2num(fields(currentFieldNum).type(1:(multLoc-1)));
                  if (~isempty(tmpMultiplier))
                      multiplier = tmpMultiplier;
                      fields(currentFieldNum).type = fields(currentFieldNum).type(multLoc+1:end);
                      fields(currentFieldNum).columns = multiplier;
                      
                  else
                     disp('Error: number of points in array not found');
                  end
                end 
                
                
                if isequal(fields(currentFieldNum).type,'uint32')
                     fields(currentFieldNum).bytesPerItem = 4;
                elseif isequal(fields(currentFieldNum).type,'int32')
                     fields(currentFieldNum).bytesPerItem = 4;
                elseif isequal(fields(currentFieldNum).type,'uint64')
                     fields(currentFieldNum).bytesPerItem = 8;
                elseif isequal(fields(currentFieldNum).type,'int64')
                     fields(currentFieldNum).bytesPerItem = 8;
                elseif isequal(fields(currentFieldNum).type,'double')
                     fields(currentFieldNum).bytesPerItem = 8;
                elseif isequal(fields(currentFieldNum).type,'uint16')
                     fields(currentFieldNum).bytesPerItem = 2;
                elseif isequal(fields(currentFieldNum).type,'int16')
                     fields(currentFieldNum).bytesPerItem = 2;
                elseif isequal(fields(currentFieldNum).type,'uint8')
                     fields(currentFieldNum).bytesPerItem = 1;
                elseif isequal(fields(currentFieldNum).type,'single')
                     fields(currentFieldNum).bytesPerItem = 4;
                else
                     error(['Field datatype not supported: ', fields(currentFieldNum).type]); 
                end
                
                bytesPerPacket = bytesPerPacket+(fields(currentFieldNum).bytesPerItem*multiplier);     
                               
                currentFieldNum = currentFieldNum+1; 
            end                           
        end
        data.fields = fields;
    else
        data.fields = fields;
        %error('Fields section of header not found.');
    end
else
    error('Error reading header section of file.');
end
      
    
      
%read in each field
byteOffset = 0;
for i = 1:length(fields)
    frewind(fid);
    
    junk = fread(fid,headersize+byteOffset,'uint8');
    skipBytes = bytesPerPacket-(fields(i).bytesPerItem*fields(i).columns);
    byteOffset = byteOffset+(fields(i).bytesPerItem*fields(i).columns);
    
    %skipBytes = 0;
    %if isequal(fields(i).type,'uint32')
    %    skipBytes = bytesPerPacket-4;
    %    byteOffset = byteOffset+4;
    %elseif isequal(fields(i).type,'uint16')
    %    skipBytes = bytesPerPacket-2;
    %    byteOffset = byteOffset+2;
    %elseif isequal(fields(i).type,'int16')
    %    skipBytes = bytesPerPacket-2;
    %    byteOffset = byteOffset+2;    
    %end
    
 
    tmpData = fread(fid,[fields(i).columns,inf],[num2str(fields(i).columns),'*', fields(i).type,'=>',fields(i).type],skipBytes);
    %if the field is time, convert to seconds
    
    if isequal(fields(i).name,'time') && (~isempty(clockRate))       
        %tmpData = double(tmpData)/clockRate;
        %data.fields(i).type = 'double';        
    end
%     if isequal(fields(i).name,'voltage') && (~isempty(voltageScaling))       
%         tmpData = double(tmpData)*voltageScaling;
%         data.fields(i).type = 'double';        
%     end
    
    %if (fields(i).length > 1)
    %  tmpData = reshape(tmpData, [fields(i).length length(tmpData)/fields(i).length])';
    %end
    data.fields(i).data = tmpData';
   
end
    

fclose(fid);


