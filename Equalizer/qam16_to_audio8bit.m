function xsint=qam16_to_audio8bit(qam_data)
% qam16_to_audio8bit.m: convert QAM16 4-bit codewords to 8-bit audio
%       samples. Each audio sample is built from a pair of input codewords
%       as (LSB,MSB), with LSB=qam_data(n), MSB=qam_data(n+1).
% Inputs:
%       qam_data: sequence of (input) QAM16 codewords
% Outputs:
%       xsint: output 8-bit audio signal (double format)

LS4=qam_data(1:2:end);
MS4=qam_data(2:2:end);
x8sint=16*MS4+LS4;
xsint=double(x8sint)-128;