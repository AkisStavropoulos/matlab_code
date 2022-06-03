function [s, flag] = INIT_Serial(comPort)
    flag = 1; 
    s = serial(comPort);
    fopen(s); 
    set(s, 'Baudrate', 115200); 

end

