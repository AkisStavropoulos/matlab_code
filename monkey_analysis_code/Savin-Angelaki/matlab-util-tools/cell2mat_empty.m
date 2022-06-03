function result = cell2mat_empty(x)
    if length(x) == 1
        result = cell2mat(x);
    else
        
        result = [];
        for j = 1:length(x)
            if ~isempty(x{j})
                break
            end
        end
        
        if length(size(x{1})) ~= 2
            ME = MException('VerifyIInput:WrongDimNumber', ...
                 'cell array must contain matrices');
            throw(ME);
        end
        
        for i = j+1:length(x)
            if ~isempty(x{i})
                break
            end
        end
        
            
        if size(x{i},1) == size(x{j},1)
            for j =1:length(x)
                result = [result; x{j}];
            end
        else
            for j =1:length(x)
                result = [result, x{j}];
            end
        end
     end  
    
end
