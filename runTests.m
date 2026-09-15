function exitCode = runTests()
    addpath('.');
    addpath('tests');
    result = runtests('tests');
    disp(result);
    if all([result.Passed])
        exitCode = 0;
    else
        exitCode = 1;
    end
end
