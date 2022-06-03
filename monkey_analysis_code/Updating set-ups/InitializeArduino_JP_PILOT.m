function [ mr, master_x1, master_x2, master_y1, master_y2 ] = InitializeArduino_JP_PILOT( RigParameters )
    delete(instrfindall)
    mr = MouseReader_2sensors(RigParameters.arduinoPort);
    master_x1 = [];
    master_x2 = [];
    master_y1  = [];
    master_y2 = [];

end

