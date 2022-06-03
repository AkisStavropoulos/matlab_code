function [d_m_nFB, disc_d_m_nFB] = check_sizes_monkey_VR(d_m_nFB, disc_d_m_nFB)

if length(d_m_nFB.stop_trial) == length(disc_d_m_nFB.trial_num)
    % do nothing
elseif length(d_m_nFB.stop_trial) > length(disc_d_m_nFB.trial_num)
    d_m_nFB.start_trial(end) = [];
    d_m_nFB.end_trial(end) = [];
    d_m_nFB.stop_trial(end) = [];
elseif length(d_m_nFB.stop_trial) < length(disc_d_m_nFB.trial_num)
    disc_d_m_nFB.trial_num(end) = [];
    disc_d_m_nFB.maxV(end) = [];
    disc_d_m_nFB.maxW(end) = [];
    disc_d_m_nFB.ffv(end) = [];
    disc_d_m_nFB.PosXo(end) = [];
    disc_d_m_nFB.PosYo(end) = [];
    disc_d_m_nFB.PosZo(end) = [];
    disc_d_m_nFB.RotXo(end) = [];
    disc_d_m_nFB.RotYo(end) = [];
    disc_d_m_nFB.RotZo(end) = [];
    disc_d_m_nFB.RotWo(end) = [];
    disc_d_m_nFB.FFx(end) = [];
    disc_d_m_nFB.FFy(end) = [];
    disc_d_m_nFB.FFz(end) = [];
    disc_d_m_nFB.pcheckX(end) = [];
    disc_d_m_nFB.pcheckY(end) = [];
    disc_d_m_nFB.pcheckZ(end) = [];
    disc_d_m_nFB.rcheckX(end) = [];
    disc_d_m_nFB.rcheckY(end) = [];
    disc_d_m_nFB.rcheckZ(end) = [];
    disc_d_m_nFB.rcheckW(end) = [];
    disc_d_m_nFB.distToFF(end) = [];
    disc_d_m_nFB.rewarded(end) = [];
    disc_d_m_nFB.timeout(end) = [];
    disc_d_m_nFB.beginTime(end) = [];
    disc_d_m_nFB.checkTime(end) = [];
    disc_d_m_nFB.duration(end) = [];
    disc_d_m_nFB.delays(end) = [];
    disc_d_m_nFB.ITI(end) = [];
end

end

