function [ x1, y1, x2, y2 ] = ReadArduino_JP_PILOT( mr )
    mr.poll_mouse();
    [x1, y1, x2,y2] = mr.get_xy_change();

end

