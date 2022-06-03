%% add trials from Unity VR setup
function AddTrials_VR(this,prs)
    cd(prs.filepath_behv);
    this.trials = AddTrials2Behaviour_VR(prs);
end