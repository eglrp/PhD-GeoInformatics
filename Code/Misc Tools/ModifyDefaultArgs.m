function ModifyDefaultArgs(arg)
if (length(arg) <= 1)
    return;
end;
for i = 1:2:length(arg)
    if (~ischar(arg{i}))
        error('arg should be formatted as: <''varname 1''>, <var value 1>,...<''varname n''>, <var value n>');
    end
    if (~evalin('caller', sprintf('exist(''%s'', ''var'')', arg{i})))
        error('Variable "%s:" does not exist', arg{i});
    end
    assignin('caller', arg{i}, arg{i+1});
end

end
