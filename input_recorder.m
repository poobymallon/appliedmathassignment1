%this class exists to store the sequence of guesses made by a root finder
%how to use:
%Step 1: Create an instance of the input_recorder:
% my_recorder = input_recorder();
%Step 2: Use input_recorder to generate a version of the test function
% that records the input after every iteration:
% If test_fun is defined as an anonymous function:
% f_record = my_recorder.generate_recorder_fun(test_fun);
% If test_fun is defined using function keyword
% f_record = my_recorder.generate_recorder_fun(@test_fun);
%Step 3: Call your root finder using the recording function:
% x_root = newtons_method(f_record, x_guess);
%Step 4: See what input values were used when f_record was called:
% input_list = my_recorder.get_input_list();
%Step 5: Clear the recorder for the next trial:
% my_recorder.clear_input_list();
%After this, you can circle back to step 3 and repeat for another trial
classdef input_recorder < handle
    properties
        %list that store function inputs
        input_list;
    end
    methods
        function obj=input_recorder()
            obj.input_list = [];
        end
        %resets input_list
        function clear_input_list(obj)
            obj.input_list = [];
        end
        %returns input_list
        function input_list = get_input_list(obj)
            input_list = obj.input_list;
        end
        %creates a version of function fun(x) that will store
        %x into the input_list whenever fun is called
        function f_record = generate_recorder_fun(obj,fun)
            f_record = @(x) obj.record_wrapper(fun,x);
        end
        function varargout = record_wrapper(obj,fun,x)
            nOutputs = nargout;
            varargout = cell(1,nOutputs);
            if nargout == 1
                a = fun(x);
                varargout{1} = a;
            end
            if nargout == 2
                [a, b] = fun(x);
                varargout{1} = a;
                varargout{2} = b;
            end
            if nargout == 3
                [a, b, c] = fun(x);
                varargout{1} = a;
                varargout{2} = b;
                varargout{3} = c;
            end
            obj.input_list(end+1) = x;
        end
    end
end
