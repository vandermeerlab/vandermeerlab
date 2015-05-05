function [ expday ] = expday(rat_id)

if strcmp(rat_id,'R066_EI')
    expday.one = 'R066-2014-11-27_recording';
    expday.two = 'R066-2014-11-28_recording';
    expday.three = 'R066-2014-11-29_recording';
    expday.four = 'R066-2014-12-01_recording';
    expday.five = 'R066-2014-12-02_recording';
    expday.six = 'R066-2014-12-03_recording';
    expday.seven = 'R066-2014-12-04_recording';
    expday.eight = 'R066-2014-12-05_recording';
end

if strcmp(rat_id,'R067_EI')
    expday.one = 'R067-2014-11-29_recording';
    expday.two = 'R067-2014-12-01_recording';
    expday.three = 'R067-2014-12-02_recording';
    expday.four = 'R067-2014-12-03_recording';
    expday.five = 'R067-2014-12-04_recording';
    expday.six = 'R067-2014-12-05_recording';
    expday.seven = 'R067-2014-12-06_recording';
    expday.eight = 'R067-2014-12-07_recording';
end

if strcmp(rat_id,'R068_EI')
    expday.one = 'R068-2014-12-01_recording';
    expday.two = 'R068-2014-12-04_recording';
    expday.three = 'R068-2014-12-05_recording';
    expday.four = 'R068-2014-12-06_recording';
    expday.five = 'R068-2014-12-07_recording';
    expday.six = 'R068-2014-12-08_recording';
    expday.seven = 'R068-2014-12-09_recording';
    expday.eight = 'R068-2014-12-10_recording';
end

if strcmp(rat_id,'R063_EI')
    expday.one = 'R063-2015-03-18_recording';
    expday.two = 'R063-2015-03-20_recording';
    expday.three = 'R063-2015-03-22_recording';
    expday.four = 'R063-2015-03-23_recording';
    expday.five = 'R063-2015-03-24_recording';
    expday.six = 'R063-2015-03-25_recording';
    expday.seven = 'R063-2015-03-26_recording';
    expday.eight = 'R063-2015-03-27_recording';
end

end