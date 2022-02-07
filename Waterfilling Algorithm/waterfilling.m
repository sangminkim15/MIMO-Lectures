function [wf] = waterfilling (p, noiselevel)

temp = length(noiselevel);
waterlevel = (p.Nt + sum(noiselevel)) / temp;
wf = waterlevel - noiselevel;

while true
    if wf >= zeros(1, length(noiselevel))
        break
    else
        wf(temp) = 0;
        temp0 = waterlevel * temp - noiselevel(temp);
        temp = temp - 1;
        waterlevel = temp0 / temp;
        
        for idx = 1 : temp
            wf(idx) = waterlevel - noiselevel(idx);
        end
    end
end

end