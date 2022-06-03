function [data, disc_data] = check_sizes_monkey_VR(data, disc_data)

if length(data.stop_trial) == length(disc_data.trial_num)
    % do nothing
elseif length(data.stop_trial) > length(disc_data.trial_num)
    data.start_trial(end) = [];
    data.end_trial(end) = [];
    data.stop_trial(end) = [];
elseif length(data.stop_trial) < length(disc_data.trial_num)
    disc_data.trial_num(end) = [];
    disc_data.maxV(end) = [];
    disc_data.maxW(end) = [];
    disc_data.ffv(end) = [];
    disc_data.PosXo(end) = [];
    disc_data.PosYo(end) = [];
    disc_data.PosZo(end) = [];
    disc_data.RotXo(end) = [];
    disc_data.RotYo(end) = [];
    disc_data.RotZo(end) = [];
    disc_data.RotWo(end) = [];
    disc_data.FFx(end) = [];
    disc_data.FFy(end) = [];
    disc_data.FFz(end) = [];
    disc_data.pcheckX(end) = [];
    disc_data.pcheckY(end) = [];
    disc_data.pcheckZ(end) = [];
    disc_data.rcheckX(end) = [];
    disc_data.rcheckY(end) = [];
    disc_data.rcheckZ(end) = [];
    disc_data.rcheckW(end) = [];
    disc_data.distToFF(end) = [];
    disc_data.rewarded(end) = [];
    disc_data.timeout(end) = [];
    disc_data.beginTime(end) = [];
    disc_data.checkTime(end) = [];
    disc_data.rewardTime(end) = [];
    disc_data.endTime(end) = [];
    disc_data.waitTime(end) = [];
    disc_data.ITI(end) = [];
end

end

