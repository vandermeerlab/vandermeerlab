function [ leds ] = leds_position( position_data_x,position_data_y,led1_position_x,led1_position_y,led2_position_x,led2_position_y,threshold )

led1_position_x = and(position_data_x < (led1_position_x + threshold), position_data_x > (led1_position_x - threshold));
led1_position_y = and(position_data_y < (led1_position_y + threshold), position_data_y > (led1_position_y - threshold));

led2_position_x = and(position_data_x < (led2_position_x + threshold), position_data_x > (led2_position_x - threshold));
led2_position_y = and(position_data_y < (led2_position_y + threshold), position_data_y > (led2_position_y - threshold));

leds1 = led1_position_x | led2_position_y;
leds2 = led2_position_x | led1_position_y;

leds = leds1 & leds2;

end