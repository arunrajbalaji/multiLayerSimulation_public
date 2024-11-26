folderList = {'/fastscratch/abalaji/N117_i0_01_V_2/'; ...
'/fastscratch/abalaji/allSpeciesNoMembrane_i0_0.1/'; ...
'/fastscratch/abalaji/allSpeciesNoMembraneWithReactions/'};

fileName = 'input.json';
parpool(length(folderList));
parfor ii = 1:length(folderList)
    main(folderList{ii}, fileName)
end
