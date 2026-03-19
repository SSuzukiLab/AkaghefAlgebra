function Z=VectAlgExamples(type,param)
    % VectAlgExamples generate various examples for experiments
    % first, specify type name and parameters following.
    arguments
        type (1,1) string {mustBeMember(type, ...
            ["cyc","KP","uqsl2+"])}
    end
    arguments(Repeating)
        param
    end
    switch type
        case "cyc"
            Z=CyclicGroupAlg.getGenerator(param{1});
        case "KP"
            Z=VectKPAlg.getGenerator;
        case "uqsl2+"
            r=param{1};
            if round(r)==r, r=1/r; end
            Z=Uqsl2BorelSmall.getGenerator(r,"L");
    end
end