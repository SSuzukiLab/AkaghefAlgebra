classdef (Abstract) HopfAlg
    % Abstract class representing a Hopf algebra interface

    methods (Abstract)
        % Abstract method for multiplication
        ret = mtimes(arg1,arg2)

        % Abstract method for comultiplication
        ret = Delta(arg)

        % Abstract method for antipode
        ret = antipode(arg, x)

        % Abstract method for unit
        ret = unit(arg)

        % Abstract method for counit
        ret = counit(arg)

        ret=algfun(arg,fun,unit)
    end
    methods
        function result = DeltaN(input, Rank)
            % Recursive function to compute the Nth coproduct (Delta)
            if Rank == 0
                % Base case: counit for N = 0
                result = counit(input);
            elseif Rank == 1
                % Base case: return the input for N = 1
                result = input;
            else
                % Ensure numSteps is an integer
                mustBeInteger(Rank);

                % Recursive call to compute DeltaN for the previous step
                splitParts = split(DeltaN(input, Rank - 1));

                % Compute Delta on the last column
                combinedResults = arrayfun(@Delta, splitParts(:, end));

                % Combine results across columns in reverse order
                for colIndex = Rank - 2:-1:1
                    combinedResults = arrayfun(@or, splitParts(:, colIndex), combinedResults);
                end

                % Sum up the results to produce the final result
                result = sum(combinedResults);
            end
        end

    end
    methods
        function ret=HP(i1,i2)
            % HP: Hopf pairing

            % disp(i1)
            % disp(i2)
            if i1.term*i2.term>1
                x1=split(i1);
                x2=split(i2);
                T=combinations(x1,x2);
                X=arrayfun(@HP,T.x1,T.x2,UniformOutput=false);
                ret=sum(vertcat(X{:}));
            elseif i1.term*i2.term==0
                ret=0;
            else
                p1=i1.pw{1};
                p2=i2.pw{1};
                if length(p1)==0
                    ret=i1.cf*counit(i2);
                elseif length(p2)==0
                    ret=i2.cf*counit(i1);
                elseif length(p1)>1||length(p2)>1
                    if length(p1)>1
                        X12=split(Delta(i2));
                        X3=i1;
                    else
                        X12=split(Delta(i1));
                        X3=i2;
                    end
                    hp1=arrayfun(@(y)HP(X3.set_cp(X3.cf,{X3.pw{1}(1)},{X3.bs{1}(1)}),y),X12(:,1),UniformOutput=false);
                    hp2=arrayfun(@(y)HP(X3.set_cp(ones(size(X3.cf)),{X3.pw{1}(2:end)},{X3.bs{1}(2:end)}),y),X12(:,2),UniformOutput=false);
                    % hp3=[X12(:,1).cf].';
                    ret=sum(prod([vertcat(hp1{:}),vertcat(hp2{:})],2));

                else
                    ret=i1.cf*i2.cf*HPGen(i1,i2);
                end
            end
            if 0
            disp("HP="+string(ret))
            disp(i1)
            disp(i2)
            end
        end
        function ret=HPGen(i1,i2)
            % 生成元に対するhopf pairngの値を出力する(通常これをオーバーライドする)
            ret=nan;
        end
        function ret=LRA(u,f)
            f12=split(Delta(f));
            X=arrayfun(@(f1,f2)HP(u,f2)*f1,f12(:,1),f12(:,2));
            ret=sum(X);

        end
    end
end
