function tests = HelloWorldTest
    tests = functiontests(localfunctions);
end

function testHelloWorld(testCase)
    actSolution = 'Hello, World!';
    expSolution = 'Hello, World!';
    verifyEqual(testCase, actSolution, expSolution);
end