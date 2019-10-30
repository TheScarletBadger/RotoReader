clear all; clc
trackpath = uigetdir([],'Select location of track files');
outpath = uigetdir([],'Select location to save output');
outname = inputdlg('Please provide name for output file');

%Set analysis criteria-----------------------------------------------------
%Time [s] between pulses before wheel is automatically considered stopped
encoder_timeout = 1;
%Minimum duration [s] of event to analyse
duration_threshold = 2;
%Encoder pulses per rotation
Encoder_PPR = 100;
%--------------------------------------------------------------------------

%Define columns in raw_data------------------------------------------------
timestamps = 1;
centre_position_x = 2;
centre_position_y = 3;
distance_to_wheel_centre_point = 4;
clockwise_pulse = 5;
anticlockwise_pulse = 6;
in_wheel_true = 7;
%--------------------------------------------------------------------------

%For the love of baphomet please dont edit this part-----------------------
D = dir(trackpath);
D = D(arrayfun(@(x) ~strcmp(x.name(1),'.'),D));
for file = 1:length(D)
    fname = D(file).name;
    f2read = [trackpath filesep fname];
    track_data = csvread(f2read,1,0);
    %void arrays which are built up in loops
    encoder_data = []; runmat = [];
    %col_1 = track_data row indices
    encoder_data(:,1) = 1:numel(track_data(:,1));
    %col_2 = Anymaze timestamps
    encoder_data(:,2) = track_data(:,timestamps);
    %col_3 = clockwise pulses
    encoder_data(:,3) = track_data(:,clockwise_pulse);
    %col_4 = anticlockwise clicks
    encoder_data(:,4) = track_data(:,anticlockwise_pulse);
    %col_5 = animal inside wheel
    encoder_data(:,5) = track_data(:,in_wheel_true);
    %remove all time points except those with pulses on encoder
    encoder_data(find(sum(encoder_data(:,3:4),2)==0),:) = [];
    %remove all time points where animal is outside wheel
    encoder_data(find(encoder_data(:,5)==0),:) = [];
    for i=2:length(encoder_data)
        tnow = encoder_data(i,2);
        tlast = encoder_data(i-1,2);
        %returns 1 if ith pulse was clockwise -1 if anticlockwise
        this_pulse = encoder_data(i,3)-encoder_data(i,4);
        %returns 1 if i-1th pulse was clockwise -1 if anticlockwise
        previous_pulse = encoder_data(i-1,3)-encoder_data(i-1,4);
        %True if this pulse was from same run as last pulse
        encoder_data(i,6) = (tnow-tlast)<encoder_timeout...
            & this_pulse==previous_pulse;
    end
    %locate start and stop points of discrete runs in encoder_data
    tsig = encoder_data(:,6)';
    dsig = diff([0 tsig 0]);
    startindex = find(dsig > 0)-1;
    endindex = find(dsig < 0)-1;
    %iterate over discrete runs
    for run = 1:length(startindex)
        %timestamp corresponding to run start
        runmat(run,1) = encoder_data(startindex(run),2);
        %timestamp corresponding to run end
        runmat(run,2) = encoder_data(endindex(run),2);
        %run duration [s]
        runmat(run,3) = runmat(run,2) - runmat(run,1);
        %number of wheel rotations during run
        runmat(run,4) = sum(sum(encoder_data(startindex(run):...
            endindex(run),3:4)))/Encoder_PPR;
        %wheel speed [rpm]
        runmat(run,5) = (runmat(run,4)/runmat(run,3))*60;
        %row indices in track_data corresponding to start and end of run
        runmat(run,6:7) = [encoder_data(startindex(run),1)...
            encoder_data(endindex(run),1)];
    end
    %delete runs with less than half a rotation of the wheel
    runmat(find(runmat(:,4)<0.5),:) = [];
    %delete runs with duration less than duration_threshold
    runmat(find(runmat(:,3)<duration_threshold),:) = [];
    %iterate over remaining runs
    for run = 1:numel(runmat(:,1))
        %track_data row indices for this run
        anymaze_inds = runmat(run,6:7);
        %distance moved [cm] by body cenre of mass during run
        dist = sum(hypot(diff(track_data(anymaze_inds(1):...
            anymaze_inds(2),centre_position_x)),diff(track_data...
                (anymaze_inds(1):anymaze_inds(2),centre_position_y))))/10;
        %distance moved by body centre [cm] per rotation of the wheel.
        runmat(run,8) = dist / runmat(run,4);
        %mean distance [cm] from wheel centre point during run
        runmat(run,9) = mean(track_data(anymaze_inds(1):anymaze_inds(2),...
            distance_to_wheel_centre_point))*100;
    end
    testname{file,1} = fname;
    if ~isempty(runmat)   
        n_run_evts(file,1) = numel(runmat(:,1));
        total_run_duration(file,1) = sum(runmat(:,3));
        max_duration(file,1) = max(runmat(:,3));
        mean_duration(file,1) = mean(runmat(:,3));
        max_rotations(file,1) = max(runmat(:,4));
        mean_rotations(file,1) = mean(runmat(:,4));
        max_speed(file,1) = max(runmat(:,5));
        mean_speed(file,1) = mean(runmat(:,5));
        max_wobb(file,1) = max(runmat(:,8));
        mean_wobb(file,1) = mean(runmat(:,8));
        max_disp(file,1) = max(runmat(:,9));
        mean_disp(file,1) = mean(runmat(:,9));
    else
        n_run_evts(file,1) = 0;
        total_run_duration(file,1) = 0;
        max_duration(file,1) = 0;
        mean_duration(file,1) = 0;
        max_rotations(file,1) = 0;
        mean_rotations(file,1) = 0;
        max_speed(file,1) = 0;
        mean_speed(file,1) = 0;
        max_wobb(file,1) = 0;
        mean_wobb(file,1) = 0;
        max_disp(file,1) = 0;
        mean_disp(file,1) = 0;
    end
end
T = table(testname, n_run_evts, total_run_duration, mean_duration,...
    max_duration, mean_rotations, max_rotations, mean_speed, max_speed,...
    mean_wobb, max_wobb, mean_disp, max_disp);

writetable(T,[outpath filesep outname{:} '.xlsx']);
