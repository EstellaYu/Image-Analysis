function [xsmooth, ysmooth] = Remove_redundancies(x,y)
%
% by W.D. Ristenpart
% March 31, 2006
%
% INPUTS 
%  x, y = vectors describing x and y coordinates of detected drop edge
%
% OUTPUT
%  xsmooth, ysmooth = vectors with redundant coordiantes removed
%  
% example: if on the detected curve there are three adjacent coordinates
% with different y values but the same x value (i.e., a straight line
% segment), this function returns a single coordinate located at x and the
% mean of the three y positions
%


% remove x redundancies
xnew=[];
ynew=[];
k=1;
while k < length(x)
    xcurrent = x(k);
    indrange=k;
    for m = k+1 : length(x)
        if x(m)==xcurrent && m<length(x)
            indrange = [indrange m];
        elseif x(m)~=xcurrent
            xnew =[xnew; xcurrent];
            ynew =[ynew; mean(y(indrange))];
            k = m;
            break;
        elseif x(m)==xcurrent && m==length(x)
            xnew =[xnew; xcurrent];
            ynew =[ynew; mean(y(indrange))];
            k = m;
            break;
        end
    end
end

% remove y redundancies
xnew2=[];
ynew2=[];
k=1;
while k <= length(ynew)
    ycurrent = ynew(k);
    indrange=k;
    if k>512
        2;
    end
    if k < length(ynew)
        for m = k+1 : length(ynew)
            if ynew(m)==ycurrent && m<length(ynew)
                indrange = [indrange m];
            elseif ynew(m)~=ycurrent
                xnew2 =[xnew2; mean(xnew(indrange))];
                ynew2 =[ynew2; ycurrent];
                k = m;
                break;
            elseif ynew(m)==ycurrent && m==length(ynew)
                xnew2 =[xnew2; mean(xnew(indrange))];
                ynew2 =[ynew2; ycurrent];
                k = m;
                break;
            end
        end
    elseif k == length(ynew)
        xnew2 =[xnew2; xnew(k)];
        ynew2 =[ynew2; ycurrent];
        k=k+1;
        break;
    end
end

% rename to return smoothed values
xsmooth = xnew2;
ysmooth = ynew2;

