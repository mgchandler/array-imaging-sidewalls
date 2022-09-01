function memorytable = fn_struct_mem_breakdown(struct, structname)
% Returns a list of fields within a structure with comparable

fields = strings(0);
mb = [];

names = fieldnames(struct);
% structinfo = whos('struct');
% fields(end+1) = structname;
% mb     = [mb, structinfo.bytes / 1024^2];
for n = 1:length(names)
    name = names{n};
    field = struct.(name);
    
    if isstruct(field)
        newtable = fn_struct_mem_breakdown(field, strcat(structname, '.', name));
        fields = [fields, newtable.fields.'];
        mb     = [mb, newtable.mb.'];
    else
        fieldinfo = whos('field');
        fields(end+1) = strcat(structname, '.', name);
        mb     = [mb, fieldinfo.bytes / 1024^2];
    end
end

fields = fields.';
mb = mb.';
memorytable = table(fields, mb);

end 