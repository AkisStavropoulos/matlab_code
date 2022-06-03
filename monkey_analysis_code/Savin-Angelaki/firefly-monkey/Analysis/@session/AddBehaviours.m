%% add behaviour
function AddBehaviours(this,prs)
    cd(prs.filepath_behv);
    this.behaviours = behaviour(prs.comments);
    if prs.unityVR
        this.behaviours.AddTrials_VR(prs);
    else
        this.behaviours.AddTrials(prs);
    end
    this.behaviours.AnalyseBehaviour(prs);
    this.behaviours.UseDatatype('single');
end