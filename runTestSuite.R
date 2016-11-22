require(RUnit)

testsuite.Kaphi <- defineTestSuite('Kaphi',
    dirs=file.path(getwd(), 'tests'),
    testFileRegexp='^test_.+\\.R',
    testFuncRegexp='^test.+'
)
testResult <- runTestSuite(testsuite.Kaphi)
printTextProtocol(testResult)
