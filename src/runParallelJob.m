function runParallelJob(inputFilePath, taskIndex, jobNumber, loopArgs)

newFolderPathList = generateDataFolders(inputFilePath, taskIndex, jobNumber, loopArgs{:});

while (~exist([newFolderPathList{taskIndex} 'inputFile_' num2str(jobNumber) '.json'], 'file'))
end

main(newFolderPathList{taskIndex}, ['inputFile_' num2str(jobNumber) '.json'])

end
