function floatdata = fix2float(fixdata, decimal_width, data_type)

if 0 == data_type
	fixdata_hex = dec2hex(fixdata, 8);
	fixdata_imag = hex2dec(fixdata_hex(:, 1 : 4 ));
	fixdata_real   = hex2dec(fixdata_hex(:, 5 : 8));
else
	fixdata_imag = int32(fixdata(:, 1 ));
	fixdata_real   = int32(fixdata(:, 2 ));
end

for idx = 1:size(fixdata, 1)
	if (fixdata_imag(idx) > 32767 )
		fixdata_imag(idx) = fixdata_imag(idx) - 65536;
	end

	if (fixdata_real(idx) > 32767 )
		fixdata_real(idx) = fixdata_real(idx) - 65536;
	end

end

format
fixdata_imag = double(fixdata_imag) ./ (2^ decimal_width);
fixdata_real = double(fixdata_real) ./ (2^ decimal_width);
floatdata =  complex(fixdata_real, fixdata_imag);

end

