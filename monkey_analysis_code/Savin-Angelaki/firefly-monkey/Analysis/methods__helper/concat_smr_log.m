function concat_str = concat_smr_log(smr,log)
   % check if there is a prs.xfp field and substitute in smr (smr will
   % either be noise or be full of zeros if xfp exists in the log)
   if isfield(log(1).prs,'xfp')
       for tr = 1:length(log)
            smr(tr).continuous.xfp = log(tr).prs.xfp.*ones(length(smr(tr).continuous.xfp), 1);
            smr(tr).continuous.yfp = log(tr).prs.yfp.*ones(length(smr(tr).continuous.yfp), 1);
       end
   end
   % get fields of smr
   fields_smr = fieldnames(smr);
   % empty struct creation
   concat_str = struct();
   for ii = 1:length(fields_smr)
       % check if field is present in both
       if isfield(log,fields_smr{ii})
            % concat structure
            concat_str.(fields_smr{ii}) = catstruct(log.(fields_smr{ii}),smr.(fields_smr{ii}));
       % otherwise store field
       else
           concat_str.(fields_smr{ii}) = smr.(fields_smr{ii});
       end
   end
   % get fields of log
   fields_log = fieldnames(log);
   for ii = 1:length(fields_log)
       % check if field is not present in concatenated and add
       if ~isfield(concat_str,fields_log{ii})
           concat_str.(fields_log{ii}) = log.(fields_log{ii});
       end
   end
   
end