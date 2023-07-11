
function z = SimModel(modelname, AllCoh, AllW, NumIter)
load I_Params.mat
%%
rng(1);
z = 0;
if modelname == 1  % TD
    for cc=1:size(AllCoh,2)
        coh = AllCoh(:,cc);
        for ww=1:size(AllW,2)
            [Conf,time,choice] = TD_model(AllW(ww), [Params_I, ConfParams_I],coh);
            Allconf(cc,ww) = mean([Conf]);
            Decision(cc,ww) = size(find(choice),1)/NumIter;
            AllRT(cc,ww) = mean([time]);
        end
        disp(['Simulation is finished for coherence = ', num2str(cc)]);
    end
    save('TD','Allconf','Decision','AllRT');
    z = 1;

end


if modelname == 2  % Bound
    for cc=1:size(AllCoh,2)
        coh = AllCoh(:,cc);
        for ww=1:size(AllW,2)
            [Conf,time,choice] = BoundModel(AllW(ww), [Params_I, ConfParams_I],coh);
            Allconf(cc,ww) = mean([Conf]);
            Decision(cc,ww) = size(find(choice),1)/NumIter;
            AllRT(cc,ww) = mean([time]);
        end
        disp(['Simulation is finished for coherence = ', num2str(cc)]);
    end
    save('Bound','Allconf','Decision','AllRT');
    z = 1;
end

if modelname == 3  % NDT
    for cc=1:size(AllCoh,2)
        coh = AllCoh(:,cc);
        for ww=1:size(AllW,2)
            [Conf,time,choice] = NDT_model(AllW(ww), [Params_I, ConfParams_I],coh);
            Allconf(cc,ww) = mean([Conf]);
            Decision(cc,ww) = size(find(choice),1)/NumIter;
            AllRT(cc,ww) = mean([time]);
        end
        disp(['Simulation is finished for coherence = ', num2str(cc)]);
    end
    save('NDT','Allconf','Decision','AllRT');
end

if modelname == 4   % Mu
    for cc=1:size(AllCoh,2)
        coh = AllCoh(:,cc);
        for ww=1:size(AllW,2)
            [Conf,time,choice] = Mu_model(AllW(ww), [Params_I, ConfParams_I],coh);
            Allconf(cc,ww) = mean([Conf]);
            Decision(cc,ww) = size(find(choice),1)/NumIter;
            AllRT(cc,ww) = mean([time]);
        end
        disp(['Simulation is finished for coherence = ', num2str(cc)]);
    end
    save('Mu','Allconf','Decision','AllRT');
    z = 1;
end


end
