classdef TestSuite < matlab.unittest.TestSuite
    methods (Static)
        function tests = suite
            tests = [tests, TestSuite.fromClass(?HelloWorldTest)];
        end
    end
end