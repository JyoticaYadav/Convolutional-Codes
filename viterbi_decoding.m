% Generating convolution code for rate 1/4 i.e. (4,1,4) code

clc;
clear all;
close all;

% Initializations - m1, m2, m3 and m4 are the memory blocks and input m1 is the input bit

indexBSC = 1;
N = 100;                              % Number of input bits generated
u = rand(N,1)>0.5;                    % Source generating bits at the probability of 0.5
u(N+1:N+4) = 0;                       % 4 bits have to be added at the end to take the state back to 0
inputM1 = 0;
outM1 = 0;                            % Available value at the output of the memory block
outM2 = 0;
outM3 = 0;
outM4 = 0;
v = zeros(N,4);

%% Convolutional Encoding

for i = 1:N+4

   inputM1 = u(i,1);

   temp2 = 0;
   conn1 = [inputM1 outM3 outM4];
   for k =1: 3                                            % v(i,1) = bitxor(inputM1, outM3, outM4);
      temp1 = conn1(1,k);
      temp2 = bitxor(temp1, temp2);                      
   end
   v(i,1) = temp2;    

   temp2 = 0;
   conn2 = [inputM1 outM1 outM2 outM4];
   for k =1: 4                                            % v(i,2) = bitxor(inputM1, outM1, outM2, outM4);
       temp1 = conn2(1,k);
       temp2 = bitxor(temp1, temp2);                      
   end
   v(i,2) = temp2;
   temp2 = 0;
   conn3 = [inputM1 outM2 outM3 outM4]; 
   for k =1: 4                                            % v(i,3) = bitxor(inputM1, outM2, outM3, outM4);                 
       temp1 = conn3(1,k);
       temp2 = bitxor(temp1, temp2);                      
   end
   v(i,3) = temp2;
   temp2 = 0;
   conn4 = [inputM1 outM1 outM3 outM4]; 
   for k =1: 4                                            % v(i,4) = bitxor(inputM1, outM1, outM3, outM4);                                                             
       temp1 = conn4(1,k);
       temp2 = bitxor(temp1, temp2);                      
   end
   v(i,4) = temp2; 
   outM4 = outM3;                                         % Shifting operation
   outM3 = outM2;
   outM2 = outM1;
   outM1 = inputM1;
end


for pBSC = 0: 0.05: 0.8
    %% Binary Symmetric channel modulation (BSC Channel)

    error = rand(size(v)) < pBSC;                             % pBSC is the probability of BER in BSC channel      
    rxdV = v;                                                 
    rxdV(error) = 1 - rxdV(error);                            % Shuffle some bits according to the probability of error specified

    %% Trellis table formation

    conn = [1, 0, 0, 1, 1;1, 1, 1, 0, 1;1, 0, 1, 1, 1;1, 1, 0, 1, 1];
    initState = dec2bin(0:2^4-1) - '0';
    table(:,2:5) = kron(initState,ones(2,1));
    for i= 2:2:32
       table(i,1)=1; 
    end

    for i = 1:32 
        for j = 1:4
            inter1 = table(i,1:5).*conn(j,:);
            table(i,10:13) = table(i,1:4);
            temp2 = 0;
            for k =1: 5                                            % v(i,1) = bitxor(inputM1, outM3, outM4);
                temp1 = inter1(1,k);
                temp2 = bitxor(temp1, temp2);                      
            end
            table(i,5+j) = temp2; 
        end
    end

    %% Viterbi Decoding for Convolutional Code

    pathMetric  = zeros(16,1);                                      % Total states are 16
    branchMetric = zeros(16,N);
    initState = [0 0 0 0];

    for iter1 = 1: N+4                                                    % Same process has to be repeated input number of times
        rowRxdV = rxdV(iter1,:);                                          % Extracting the received vector for the first iteration
        for iter2 = 1:16                                                  % For loop for number of states 
            binStateNum = decimalToBinaryVector(iter2-1,4);
            [~,indInitState1] = ismember(binStateNum,table(:,10:13),'rows');    % First row which has an output state as iter2-1
            indInitState2 = indInitState1+2;
            trackTable(iter2,1:4) = table(indInitState1,2:5);
            trackTable(iter2,5:8) = table(indInitState2,2:5);                   % initial states saved

            output1 = table(indInitState1,6:9);                                 % Hamming distance calculations
            output2 = table(indInitState2,6:9);
            hammDist1 = sum(xor(rowRxdV, output1),2);
            hammDist2 = sum(xor(rowRxdV, output2),2);
            trackTable(iter2,9) = hammDist1;
            trackTable(iter2,10) = hammDist2;

            inter1 = trackTable(iter2,1)*8 + trackTable(iter2,2)*4 + trackTable(iter2,3)*2 +trackTable(iter2,4)*1 +1;
            inter2 = trackTable(iter2,5)*8 + trackTable(iter2,6)*4 + trackTable(iter2,7)*2 +trackTable(iter2,8)*1 +1;
            bm1 = branchMetric(inter1,iter1)+ hammDist1;
            bm2 = branchMetric(inter2,iter1)+ hammDist2;
            trackTable(iter2, 11) = bm1;
            trackTable(iter2, 12) = bm2;

            if iter1 == 1                                               % Talking about the first input bit
                if iter2 == 1 || iter2 == 9
                    pathMetric = [1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0];
                    branchMetric(iter2,iter1+1) = bm1;
                else
                end
            end
            if iter1 == 2
                if iter2 == 1 || iter2 == 5 || iter2 == 9 || iter2 == 13
                   pathMetric = [1;0;0;0;9;0;0;0;1;0;0;0;9;0;0;0];  
                   branchMetric(iter2,iter1+1) = bm1;
                else
                end
            else
            end
            if  iter1 == 3
                 if iter2 == 1 || iter2 == 3 || iter2 == 5 || iter2 == 7|| iter2 == 9 || iter2 == 11|| iter2 == 13|| iter2== 15
                    pathMetric = [1;0;5;0;9;0;13;0;1;0;5;0;9;0;13;0];  
                    branchMetric(iter2,iter1+1) = bm1;
                 else
                 end 
            else
            end

            if iter1 == 4
                pathMetric = [1;3;5;7;9;11;13;15;1;3;5;7;9;11;13;15];
                branchMetric(iter2,iter1+1) = bm1;
            else
            end

            if iter1 > 4   
                if bm1 < bm2
                   pathMetric(iter2,1)= inter1;
                   branchMetric(iter2,iter1+1)=bm1;
                else
                   pathMetric(iter2,1)= inter2;  
                   branchMetric(iter2,iter1+1)=bm2;
                end
            end
        end
        pathMetric_N(:,iter1) = pathMetric(:,1);
    end

    % trace back unit
    finalState = 1;                                           % current state is 1
    for var = N+4:-1:1
        initState = pathMetric_N(finalState,var);
        row = 2*initState-1;
        fStVec = table(row,10:13);
        decEqFinalState = fStVec(1,1)*8+ fStVec(1,2)*4+ fStVec(1,3)*2+fStVec(1,4)*1;
        if (decEqFinalState + 1) == finalState
            input(var,1)=0;
        else
            input(var,1)=1;
        end
        finalState = initState;
    end    
        
    errors = abs(input - u);
    BERConv(indexBSC,1) = sum(errors);
    
    %% Viterbi decoding for punctured convolutional code

    pathMetricPunct  = zeros(16,1);                                   % Total states are 16
    branchMetricPunct = zeros(16,N);
    initState = [0 0 0 0];

    var_check = 0;

    for iter1 = 1: N+4                                                    % Same process has to be repeated input number of times

        for iter2 = 1:16                                                  % For loop for number of states 
            binStateNum = decimalToBinaryVector(iter2-1,4);
            [~,indInitState1] = ismember(binStateNum,table(:,10:13),'rows');     % First row which has an output state as iter2-1
            indInitState2 = indInitState1+2;
            trackTable(iter2,1:4) = table(indInitState1,2:5);
            trackTable(iter2,5:8) = table(indInitState2,2:5);                   % initial states saved

            if var_check < 2
                rowRxdV = rxdV(iter1,1:3);
                output1 = table(indInitState1,6:8);                                 % Hamming distance calculations
                output2 = table(indInitState2,6:8);
                hammDist1 = sum(xor(rowRxdV, output1),2);
                hammDist2 = sum(xor(rowRxdV, output2),2);
                trackTable(iter2,9) = hammDist1;
                trackTable(iter2,10) = hammDist2;

                inter1 = trackTable(iter2,1)*8 + trackTable(iter2,2)*4 + trackTable(iter2,3)*2 +trackTable(iter2,4)*1 +1;
                inter2 = trackTable(iter2,5)*8 + trackTable(iter2,6)*4 + trackTable(iter2,7)*2 +trackTable(iter2,8)*1 +1;
                bm1 = branchMetricPunct(inter1,iter1)+ hammDist1;
                bm2 = branchMetricPunct(inter2,iter1)+ hammDist2;
                trackTable(iter2, 11) = bm1;
                trackTable(iter2, 12) = bm2;

                if iter1 == 1                                               % Talking about the first input bit
                    if iter2 == 1 || iter2 == 9
                        pathMetricPunct = [1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0];
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end
                end
                if iter1 == 2
                    if iter2 == 1 || iter2 == 5 || iter2 == 9 || iter2 == 13
                        pathMetricPunct = [1;0;0;0;9;0;0;0;1;0;0;0;9;0;0;0];  
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end
                else
                end
                if  iter1 == 3
                    if iter2 == 1 || iter2 == 3 || iter2 == 5 || iter2 == 7|| iter2 == 9 || iter2 == 11|| iter2 == 13|| iter2== 15
                        pathMetricPunct = [1;0;5;0;9;0;13;0;1;0;5;0;9;0;13;0];  
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end 
                else
                end

                if iter1 == 4
                    pathMetricPunct = [1;3;5;7;9;11;13;15;1;3;5;7;9;11;13;15];
                    branchMetricPunct(iter2,iter1+1) = bm1;
                else
                end

                if iter1 > 4   
                    if bm1 < bm2
                        pathMetricPunct(iter2,1)= inter1;
                        branchMetricPunct(iter2,iter1+1)=bm1;
                    else
                        pathMetricPunct(iter2,1)= inter2;  
                        branchMetricPunct(iter2,iter1+1)=bm2;
                    end
                end
                var_check = var_check + 1;
            end
            if (var_check >= 2   && var_check < 4)
                rowRxdV = rxdV(iter1,1:2);
                output1 = table(indInitState1,6:7);                                 % Hamming distance calculations
                output2 = table(indInitState2,6:7);
                hammDist1 = sum(xor(rowRxdV, output1),2);
                hammDist2 = sum(xor(rowRxdV, output2),2);
                trackTable(iter2,9) = hammDist1;
                trackTable(iter2,10) = hammDist2;

                inter1 = trackTable(iter2,1)*8 + trackTable(iter2,2)*4 + trackTable(iter2,3)*2 +trackTable(iter2,4)*1 +1;
                inter2 = trackTable(iter2,5)*8 + trackTable(iter2,6)*4 + trackTable(iter2,7)*2 +trackTable(iter2,8)*1 +1;
                bm1 = branchMetricPunct(inter1,iter1)+ hammDist1;
                bm2 = branchMetricPunct(inter2,iter1)+ hammDist2;
                trackTable(iter2, 11) = bm1;
                trackTable(iter2, 12) = bm2;

                if iter1 == 1                                               % Talking about the first input bit
                    if iter2 == 1 || iter2 == 9
                        pathMetricPunct = [1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0];
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end
                end
                if iter1 == 2
                    if iter2 == 1 || iter2 == 5 || iter2 == 9 || iter2 == 13
                        pathMetricPunct = [1;0;0;0;9;0;0;0;1;0;0;0;9;0;0;0];  
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end
                else
                end
                if  iter1 == 3
                    if iter2 == 1 || iter2 == 3 || iter2 == 5 || iter2 == 7|| iter2 == 9 || iter2 == 11|| iter2 == 13|| iter2== 15
                        pathMetricPunct = [1;0;5;0;9;0;13;0;1;0;5;0;9;0;13;0];  
                        branchMetricPunct(iter2,iter1+1) = bm1;
                    else
                    end 
                else
                end

                if iter1 == 4
                    pathMetricPunct = [1;3;5;7;9;11;13;15;1;3;5;7;9;11;13;15];
                    branchMetricPunct(iter2,iter1+1) = bm1;
                else
                end

                if iter1 > 4   
                    if bm1 < bm2
                        pathMetricPunct(iter2,1)= inter1;
                        branchMetricPunct(iter2,iter1+1)=bm1;
                    else
                        pathMetricPunct(iter2,1)= inter2;  
                        branchMetricPunct(iter2,iter1+1)=bm2;
                    end
                end
                var_check = var_check + 1;
            end
            if var_check == 4
                var_check = 0;
            end
        end
        pathMetric_N_Punct(:,iter1) = pathMetricPunct(:,1);
    end

    % trace back unit
    finalState = 1;                                           % current state is 1
    for var = N+4:-1:1
        initState = pathMetric_N_Punct(finalState,var);
        row = 2*initState-1;
        fStVec = table(row,10:13);
        decEqFinalState = fStVec(1,1)*8+ fStVec(1,2)*4+ fStVec(1,3)*2+fStVec(1,4)*1;
        if (decEqFinalState + 1) == finalState
            input_Punct(var,1)=0;
        else
            input_Punct(var,1)=1;
        end
        finalState = initState;
    end    
    
    errors = abs(input_Punct - u);
    BERPunc(indexBSC,1) = sum(errors);
    indexBSC = indexBSC + 1;
end

%% Plotting 

Conv_Viterbi = log10(BERConv/N);                     % Viterbi decoding BER
Punct_Viterbi = log10(BERPunc/N);

pBSC = 0:0.05:0.8;
figure
plot(pBSC, BERConv/N)
hold on
grid on
plot(pBSC, BERPunc/N)
xlabel('pBSC');
ylabel('Bit Error Rate');
title('BER vs bit flipping probability')
legend('Convolutional Codes','Punctured Convolutional Codes','Location','northwest')

