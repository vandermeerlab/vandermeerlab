% NLX2MATNRD Imports data from Neuralynx NRD files to Matlab variables.
%
%   [Timestamps, Samples, Header] = Nlx2MatNRD( Filename, ChannelNumber,
%                                  FieldSelectionFlags, HeaderExtractionFlag,
%                                  ExtractMode, ExtractionModeVector);
%
%   Version 5.0.1 
%
%   INPUT ARGUMENTS:
%   FileName: String containing either the complete ('C:\CheetahData\
%             RawData.nrd') or relative ('RawData.nrd') path of the file you wish
%             to import. 
%   ChannelNumber: The zero based channel number in the NRD file whose data
%                  will be imported(i.e. channel 0 is the first channel,
%                  channel 1 is the second channel, etc).
%   FieldSelectionFlags: Vector with each item being either a zero (excludes
%                        data) or a one (includes data) that determines which
%                        data will be returned for each record. The order of
%                        the items in the vector correspond to the following:
%                           FieldSelectionFlags(1): Timestamps
%                           FieldSelectionFlags(2): Samples
%                        EXAMPLE: [1 0] imports timestamp data from each record
%                        and excludes all other data.
%   HeaderExtractionFlag: Either a zero if you do not want to import the header
%                         or a one if header import is desired..
%   ExtractionMode: A number indicating how records will be processed during
%                   import. The numbers and their effect are described below:
%                      1 (Extract All): Extracts data from every record in
%                        the file.
%                      2 (Extract Record Index Range): Extracts every record
%                        whose index is within a range.
%                      3 (Extract Record Index List): Extracts a specific list
%                        of records based on record index.
%                      4 (Extract Timestamp Range): Extracts every record whose
%                        timestamp is within a range of timestamps.
%                      5 (Extract Timestamp List): Extracts a specific list of
%                        records based on their timestamp.
%   ExtractionModeVector: The contents of this vector varies based on the
%                         ExtractionMode. Each extraction mode is listed with
%                         a description of the ExtractionModeVector contents.
%                      1 (Extract All): The vector value is ignored.
%                      2 (Extract Record Index Range): A vector of two indices,
%                        in increasing order, indicating a range of records to
%                        extract. A record index is the number of the record in
%                        the file in temporal order (i.e. first record is index
%                        1, second is 2, etc.). This range is inclusive of the
%                        beginning and end indices. If the last record in the
%                        range is larger than the number of records in the
%                        file, all records until the end of the file will be
%                        extracted.
%                        EXAMPLE: [10 50] imports the 10th record through the
%                        50th record (total of 41 records) of the file.
%                      3 (Extract Record Index List): A vector of indices
%                        indicating individual records to extract. A record
%                        index is the number of the record in the file in
%                        temporal order (i.e. first record is index
%                        1, second is 2, etc.). Data will be extracted in the
%                        order specified by this vector. If an index in the
%                        vector is less than 1 or greater than the number of
%                        records in the file, the index will be skipped.
%                        EXAMPLE: [7 10 1] imports record 7 then 10 then 1,
%                        it is not sorted temporally
%                      4 (Extract Timestamp Range): A vector of two timestamps,
%                        in increasing order, indicating a range of time to use
%                        when extracting records. If either of the timestamps
%                        in the vector are not contained within the timeframe
%                        of the file, the range will be set to the closest
%                        valid timestamp (e.g. first or last). The range is
%                        inclusive of the beginning and end timestamps. If a
%                        specified timestamp occurs within a record, the entire
%                        record will be extracted. This means that the first
%                        record extracted may have a timestamp that occurs
%                        before the specified start time.
%                        EXAMPLE: [12500 25012] extracts all records that
%                        contain data that occurred between the timestamps
%                        12500 and 25012, inclusive of data at those times.
%                      5 (Extract Timestamp List): A vector of timestamps
%                        indicating individual records to extract. If a
%                        specified timestamp occurs within a record, the entire
%                        record will be extracted. This means that the a record
%                        extracted may have a timestamp that occurs before the
%                        specified timestamp. If there is no data available for
%                        a specified timestamp, the timestamp will be ignored.
%                        Data will be retrieved in the order specified by this
%                        vector.
%                        EXAMPLE: [45032 10125 75000] imports records that
%                        contain data that occurred at timestamp 45035 then
%                        10125 then 75000, it is not sorted temporally.
%
%   Notes on output data:
%   1. Each output variable's Nth element corresponds to the Nth element in
%      all the other output variables with the exception of the header output
%      variable.
%   2. The value of N in the output descriptions below is the total number of
%      records extracted.
%   3. For more information on Neuralynx records see:
%      http://www.neuralynx.com/static/software/NeuralynxDataFileFormats.pdf
%   4. Output data will always be assigned in the order indicated in the
%      FieldSelectionFlags. If data is not imported via a FieldSelectionFlags
%      index being 0, simply omit the output variable from the command.
%      EXAMPLE: FieldSelectionFlags = [1 0];
%      [Timestamps] = Nlx2MatNRD('test.nrd', 3, FieldSelectionFlags,0,1,[]);
%
%   OUTPUT VARIABLES:
%   Timestamps: A 1xN vector of timestamps. These timestamps  are the full
%               64 bit Cheetah timestamps.
%   Samples: A 1xN vector of the data points. These values are in AD counts.
%   Header: A Mx1 vector of all the text from the Neuralynx file header, where
%           M is the number of lines of text in the header.
%
%
%   EXAMPLE: [Timestamps, Samples, Header] = Nlx2MatNRD('test.nrd', 3,
%                                            [1 1], 1, 1, [] );
%   Uses extraction mode 1 to return all of the data from all of the records of
%   AD channel 3 in the file test.nrd.
%
%   EXAMPLE: [Timestamps, Header] = Nlx2MatNRD('test.nrd',[1 0],1, 2, [14 30]);
%   Uses extraction mode 2 to return the Timestamps between record index 14 and
%   30 as well as the complete file header.
%
