function qam_data=audio8bit_to_qam16(x)
% audio8bit_to_qam16.m: convert 8-bit audio samples to QAM16 codewords.
%       Each audio sample yields a pair of input 4-bit codewords.
% Inputs:
%       x: output 8-bit audio signal (double format)
% Outputs:
%       qam_data: sequence of (input) QAM16 codewords

x8=uint8(x+128);
MS4_org=bitand(x8,uint8(240))/16;
LS4_org=bitand(x8,uint8(15));
qam_data(1:2:2*length(LS4_org))=double(LS4_org);
qam_data(2:2:2*length(LS4_org))=double(MS4_org);
qam_data=qam_data';