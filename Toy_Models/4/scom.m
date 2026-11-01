function line = scom(fid)
% SCOM Legge la riga successiva ignorando commenti che iniziano con '!' e righe vuote
line = fgetl(fid);
if ~ischar(line)
    line = '';
    return;
end

tmp = sscanf(line, '%s', 1);
while isempty(tmp) || tmp(1) == '!'
    if feof(fid)
        line = '';
        return;
    end
    line = fgetl(fid);
    if ~ischar(line)
        line = '';
        return;
    end
    tmp = sscanf(line, '%s', 1);
end
