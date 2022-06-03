function y = cell2mat_singleDouble(x)
    y = cell(size(x));
    for i =1:numel(x)
         y{i} = double(x{i});
    end
    y = cell2mat(y);
end