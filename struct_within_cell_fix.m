function varargout = struct_within_cell_fix(varargin)

%% Moves the cell inside the fields of the struct for easier use
for k = 1:length(varargin)
    currentvar = varargin{k};
    sz = size(currentvar);
    structlength = length(currentvar{1});
    
    fields = fieldnames(currentvar{1});
    temp = [];
    for f = 1:length(fields)
        for i = 1:sz(1)
            for j = 1:sz(2)
                for n = 1:structlength
                    temp.(fields{f}){i,j}(n,:) = currentvar{i,j}(n).(fields{f});
                end
            end
        end
    end
    varargout{k} = temp;
end


