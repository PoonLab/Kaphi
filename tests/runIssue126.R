require(RUnit, quietly=TRUE)

testResult <- runTestFile("tests/test_issue126.R")
printTextProtocol(testResult)
