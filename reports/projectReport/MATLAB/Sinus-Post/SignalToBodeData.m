%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%         Function to process data (magnitude, phase lag)       %%
%-----------------------------------------------------------------%
%                     File: SignalToBodeData.m                    %
%                       Author: Peter Munch                       %
%                   E-Mail: peterrmuench@aol.com                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                  %
function [ f0,fin,fout,w,a,p,c,e1,e2 ]...                         %
                       = SignalToBodeData( file,file2f,in,out )   %
                                                                  %
    % Read in data from file                                      %
    fid=fopen(file);                                              %
    inenum = 1;                                                   %
    C = textscan(fid, '%d', 1, 'delimiter', '\n',...              %
                        'headerlines', inenum-1);                 %
    c = cell2mat(C);                                              %
    M = dlmread(file,'\t',c+1,0);                                 %
                                                                  %
    % Extract from data time-, input- and output-data             %
    TIME = M(:,1); TIME = TIME-min(TIME);                         %
    IN   = M(:,in);                                               %
    OUT  = M(:,out);                                              %
                                                                  %
    % Perform FFT on in- and out-signal & ignore average value    %
    X=fft(IN); X =X(2:end);                                       %
    Y=fft(OUT);Y =Y(2:end);                                       %
                                                                  %
    % Determin the mode with the hightest amplitude               %
    [mag_x, idx_x] = max(abs(X));                                 %
    [mag_y, idx_y] = max(abs(Y));                                 %
                                                                  %
    % Determin the phase lag between these modes                  %
    px = angle(X(idx_x));                                         %
    py = angle(Y(idx_y));                                         %
                                                                  %
    phase_lag = (py - px)*180/pi;                                 %
                                                                  %
    if phase_lag>0                                                %
        phase_lag  = phase_lag-360;                               %
    end                                                           %
                                                                  %
    p = phase_lag;                                                %
                                                                  %
    % Determin the amplitude ratio of these modes                 %
    a = mag_y/mag_x;                                              %
                                                                  %
    % Determin the frequency: from the name of the input file     %
    f0 = file2f(file)^-1;                                         %
    w  = 2*pi*f0;                                                 %
                                                                  %
    % Determin the frequency: from the dominating input mode      %
    fs=(TIME(end)-TIME(1))^-1*(length(IN)-1);                     %
    fin = (idx_x)*fs/(length(IN)-1);                              %
                                                                  %
    % Determin the frequency: from the dominating output mode     %
    fout = (idx_y)*fs/(length(OUT)-1);                            %
                                                                  %
    % Determin coherence of input and outpu                       %
    ct = mscohere(IN,OUT);                                        %
    c = sum(ct)/length(ct);                                       %
                                                                  %
    % Determin error                                              %
    e1 = geterror(IN);                                            %
    e2 = geterror(OUT);                                           %
                                                                  %
    % Function to determin errror                                 %
    function e = geterror(vec)                                    %
        os=vec;                                                   %
        % Perfom FFT                                              %
        temp = fft(vec);                                          %
        % Filter signal                                           %
        [sortedX,sortingIndices] = sort(temp,'descend');          %
        temp = temp * 0;                                          %
        temp(sortingIndices(1:3))=sortedX(1:3);                   %
        % Perfom iFFT on filtered signal                          %
        ns=ifft(temp);                                            %
        % Determin error                                          %
        e=max((ns-os)./os);                                       %
    end                                                           %
                                                                  %
end                                                               %
                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%