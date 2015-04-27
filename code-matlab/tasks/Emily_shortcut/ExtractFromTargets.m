%ExtractFromTargets Extract various values from array of neuralynx video tracker led/diode/bright spot bitfields ("targets" or "points").
%   [x,y,color,valid_targets] = ExtractFromTargets( nlx_targets ), for M-by-N matrix nlx_targets.
%   nlx_targets is composed of M number nlx records by N max targets per record.
%   The function returns 4 parameters:
%       x - an M by N matrix where M = the number of records and N = number of x coordinates per record.
%           this value represents the x coordinate of a given target.
%       y - an M by N matrix where M = the number of records and N = number of y coordinates per record.
%           this value represents the y coordinate of a given target.
%       color - an M by N by 7 matrix where M = the number of records and N = number of colors per record.
%           The 7 represents the 7 types of color information that is stored.
%           this value represents the color of a given target, or zero if that color is not present, or one if present.
%           color(:,:,1) = Pure Red
%           color(:,:,2) = Pure  Green
%           color(:,:,3) = Pure Blue
%           color(:,:,4) = Raw Red
%           color(:,:,5) = Raw Green
%           color(:,:,6) = Raw Blue
%           color(:,:,7) = Intensity
%       valid_targets - a matrix containing the number of values for each N in the above variables.
%           this value represents the number of targets per record (which varies from record to record).
%
%   Example:
%       Import the data from an NVT file:
%	    >>[ts, x, y, angle, targets, points] = Nlx2MatVT('VT1.nvt', [1,1,1,1,1,1], 0, 1);
%
%       Due to Matlab limitations, processing more than about 1000 targets will take an extremely long time.  Breaking
%		the imported targets into smaller matricies is recommended.
%       >> small_target_matrix  = targets(:,1:1000);
%
%     You can then call this function with the smaller targets matrix.
%       >> [x, y, color, valid_targets] = ExtractFromTargets(small_target_matrix);
%
% v1.1.0

function [x,y,color,valid_targets] = ExtractFromTargets( nlx_targets )

    % find how many recs and how many targets per record
    [max_targets,max_records] = size( nlx_targets );

    %output number of rec in file
    str = sprintf('There are %d records in this file',max_records);
    disp(str);
	
	%initialize the variables in case no targets are found
	x(1) = 0;
	y(1) = 0;
	color(1) = 0;
	valid_targets(1) = 0;
    % loop to extract all targets from each record
    for rec_index = 1:max_records

        %get record
        current_record = nlx_targets(:,rec_index);

        % loop to extract all targets from current record
        for target_index = 1:max_targets
            
            % if the bitfield is equal to zero, then we know there is no valid data for that
            % field, or the rest of the bitfields in the record.
            if current_record( target_index ) == 0
                break;
            end

            % extract the x and y positions and store them to be returned to matlab
            [x( rec_index, target_index ),y( rec_index, target_index )] = ExtractXY( current_record( target_index ) );
            [color( rec_index, target_index, 1:7 )] = ExtractColor( current_record( target_index ) );
            
        end  %end inner loop

        % record the number of targets within each rec that contain data
        valid_targets( rec_index ) = target_index - 1;
            

    end  %end loop
    
    disp('All records processed.');
end

%----------------------------------------------------------------------------------------------------------------------
%   This function extracts the x and y coordinates from the bitfield for a given target.  
%----------------------------------------------------------------------------------------------------------------------
function [x, y] = ExtractXY(target_value)  

	binary_target_string = dec2bin(target_value, 32);
	x = bin2dec(binary_target_string(21:32));
	y = bin2dec(binary_target_string(5:16));
  end
  
%----------------------------------------------------------------------------------------------------------------------
%	Extracts color information from a target
%----------------------------------------------------------------------------------------------------------------------
function [color] = ExtractColor(target_value)

	binary_target_string = dec2bin(target_value, 32);
	color(1) = bin2dec(binary_target_string(2)); %pure red
	color(2) = bin2dec(binary_target_string(3)); %pure green
	color(3) = bin2dec(binary_target_string(4)); %pure blue
	color(4) = bin2dec(binary_target_string(18)); %raw red
	color(5) = bin2dec(binary_target_string(19)); %raw green
	color(6) = bin2dec(binary_target_string(20)); %raw blue
	color(7) = bin2dec(binary_target_string(17)); %intensity

 end

