function [image,xyz,filename,s,res] = OpenMaskXML(filename)
%% function [image,xyz,filename,s,res] = OpenMaskXML(filename)
% opens *Mask.xml files

if nargin == 0
    [filename,PathName] = uigetfile('*Mask.xml','Alte Kalibrierung wählen');
    filename = [PathName, filename];
end

%filename = '326103009_M12Slices.xml';
%filename = '326103009_bla_M12Slices.xml';
%filename = '326103009_M12_rotatedSlices.xml';
fid = fopen(filename,'r','l');
s = zeros(3,1);
xyz = zeros(3,2);
res = zeros(3,1);
% xml auslesen:
while true
    tline = fgetl(fid);
    % auf Tags überprüfen:
    [isSubstr,number] = GetSubstr(tline,'<SizeX>','</SizeX>');
    if isSubstr
        s(1) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<SizeY>','</SizeY>');
    if isSubstr
        s(2) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<SizeZ>','</SizeZ>');
    if isSubstr
        s(3) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<X1>','</X1>');
    if isSubstr
        xyz(1,1) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<X2>','</X2>');
    if isSubstr
        xyz(1,2) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<Y1>','</Y1>');
    if isSubstr
        xyz(2,1) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<Y2>','</Y2>');
    if isSubstr
        xyz(2,2) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<Z1>','</Z1>');
    if isSubstr
        xyz(3,1) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<Z2>','</Z2>');
    if isSubstr
        xyz(3,2) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<SpacingX>','</SpacingX>');
    if isSubstr
        res(1) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<SpacingY>','</SpacingY>');
    if isSubstr
        res(2) = number;
    end
    [isSubstr,number] = GetSubstr(tline,'<SpacingZ>','</SpacingZ>');
    if isSubstr
        res(3) = number;
    end
    % abbrechen
    if (strfind(tline,'</Header>')),   break,   end
end

% Bild auslesen:
A = fread(fid,prod(s),'*int8');
fclose(fid);

A = reshape(A,s');
image = zeros(s(2),s(1),s(3));
for i=1:s(3)
    image(:,:,i) = A(:,:,i)';
end
end

function [isSubstr,number] = GetSubstr(line,startTag,endTag)
isSubstr = false;
number = 0;
if (strfind(line,startTag))
    if (strfind(line,endTag))
    startPos = strfind(line,startTag);
    endPos = strfind(line,endTag)-1;
    len = length(startTag);
    substr = line(startPos+len:endPos);
    number = str2double(substr);
    isSubstr = true;
    end
end
end