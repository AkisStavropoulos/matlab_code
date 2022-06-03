%% Submit jobs to NYU Greene cluster
clear c

% Get a handle to the cluster
c = parcluster;
% Specify the walltime (e.g., 5 hours)
c.AdditionalProperties.WallTime = '10:00:00';

% Specify e-mail address to receive notifications about your job

c.AdditionalProperties.EmailAddress = 'ges6@nyu.edu';
% Specify number of GPUs of a particular GPU card
c.AdditionalProperties.GpusPerNode = 1;
c.AdditionalProperties.GpuCard = 'gpu-card';
% Specify memory to use for MATLAB jobs, per core (default: 4gb)
c.AdditionalProperties.MemUsage = '35gb';
% Specify a queue to use for MATLAB jobs
c.AdditionalProperties.QueueName = 'queue-matlab';
% Request entire nodes (default: false)
c.AdditionalProperties.RequireExclusiveNode = true;
% Specify reservation to run under
c.AdditionalProperties.Reservation = 'reservation-data-extraction';

c.saveProfile


%% Send independent batch jobs

% Get a handle to the cluster
c = parcluster;
% Submit job to query where MATLAB is running on the cluster
job = c.batch(@pwd, 1, {}, 'CurrentFolder','.', 'AutoAddClientPath',false);
% Query job for state
job.State
% If state is finished, fetch the results
job.fetchOutputs{:}
% Delete the job after results are no longer needed
job.delete