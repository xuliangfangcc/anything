% function phase = phase_compensation(cent_freq, DL, sample_rate, scs)

%% STO
function [phase, a_angel] = STO_phase_compensation(DL)
scs = 2;%1--15k,2---30k
FIX_COE = pow2(14);
T0_symb_table = ...
[
[320, 4704, 9088, 13472, 17856, 22240, 26624, 31040, 35424, 39808, 44192, 48576, 52960, 57344];%15k
[352, 4736, 9120, 13504, 17888, 22272, 26656, 31040, 35424, 39808, 44192, 48576, 52960, 57344];%30k
];

cent_freq = 620.4;%MHz
sample_rate = 122.880;%MHz
% sita = cent_freq/sample_rate*T0_symb_table(2,:)
% phasecomple= (cos(-2*pi*sita) + 1j*sin(-2*pi*sita))
% DL = 1;
if DL == 1
   a_internal = - 2 * pi * cent_freq / sample_rate;
else
   a_internal =  2 * pi * cent_freq / sample_rate;
end
phase = zeros(1,14);
a_angel = zeros(1,14);
for symIdx = 1:14
    a_angel(symIdx) = a_internal * T0_symb_table(scs,symIdx);

    phase(symIdx) = exp(1j*a_angel(symIdx));
    % sin 定点值放到高16bit, cos定点值放到低16bit
%     a_sin(symIdx) = sin(a_angel(symIdx));
%     a_cos(symIdx) = cos(a_angel(symIdx));
 
%     phase_high16bit = bitand(bitshift(round(a_sin(symIdx) * FIX_COE), 16 ), 0xffff0000);
%     phase_low16bit  = bitand(round(a_cos(symIdx) * FIX_COE), 0x0000ffff,'uint32');
%     phase(symIdx) = bitor(phase_high16bit, phase_low16bit);

end

end


