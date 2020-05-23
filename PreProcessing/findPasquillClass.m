function [ classP ] = findPasquillClass( simCARt , Rsw )
%% findPasquillClass find the potential Pasquill classes for each time step
%   Inputs:
%   simCARt - hourly cloud cover data [-]
%   Rsw - hourly mean incoming shortwave radiation data [W m^-2]
%   Output:
%   classP - Pasquill classes
%% Stability parameter (Pasquill classes)
% Modified after Mohan and Siddiqui (1998).
for i=1:8760 % Assuming first and last hourly time step are at night time
    if i==1 || i==8760
        if simCARt(i)>0.875
            classP{i,1}=4; % Class D
        else
            classP{i,1}=[4 5 6]; % Classes D-F
        end
    else
        if Rsw(i)>0 % Day
            if Rsw(i)<50
                classP{i,1}=[3 4]; % Classes C-D
            elseif Rsw(i)>=50 && Rsw(i)<300
                classP{i,1}=[2 3 4];
            elseif Rsw(i)>=300 && Rsw(i)<600
                classP{i,1}=[1 2 3 4];
            else
                classP{i,1}=[1 2 3];
            end
        elseif Rsw(i)==0 && Rsw(i-1)==0 && Rsw(i+1)==0 % Night
            if simCARt(i)>0.875
                classP{i,1}=4;
            else
                classP{i,1}=[4 5 6];
            end
        else % Transition between night and day (1 h)
            classP{i,1}=4; % Class D
        end
    end
end
end