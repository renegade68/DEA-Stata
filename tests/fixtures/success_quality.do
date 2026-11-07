run "@ADO@/_dea_simplex.ado"
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125,256\8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25)
b=(7340032\6.67572021484375e-06\7\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0081_result.bin",(R.status,R.objval,R.iters));savebin("s0081_x.bin",R.x);savebin("s0081_audit.bin",_rsm_audit_read())
fh=fopen("s0081_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0081_eqmat.bin",E)
savebin("s0081_eqrhs.bin",b)
savebin("s0081_ineqmat.bin",U)
savebin("s0081_ineqrhs.bin",u)
savebin("s0081_cost.bin",c)
end
mata: savebin("s0081_success.bin",_rsm_success_read());savebin("s0081_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322,262144.00024414062\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\8.0000000037252903,4095.9999990463257,49152.000015258789,384,327679.99993896484,32,262144\8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322,262144.00024414062)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0084_result.bin",(R.status,R.objval,R.iters));savebin("s0084_x.bin",R.x);savebin("s0084_audit.bin",_rsm_audit_read())
fh=fopen("s0084_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0084_eqmat.bin",E)
savebin("s0084_eqrhs.bin",b)
savebin("s0084_ineqmat.bin",U)
savebin("s0084_ineqrhs.bin",u)
savebin("s0084_cost.bin",c)
end
mata: savebin("s0084_success.bin",_rsm_success_read());savebin("s0084_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999990463257,12582912.00390625,1.4999999997671694,327680,8192.0000076293945,1024.0000009536743\2.9103830456733704e-11,3.814697265625e-06,0.01171875,1.3969838619232178e-09,0.00030517578125,7.62939453125e-06,9.5367431640625e-07\0.031250000014551915,4095.9999990463257,12582912.00390625,1.5,327679.99993896484,8192,1024\0.03125,4095.9999990463257,12582912.00390625,1.4999999997671694,327680,8192.0000076293945,1024.0000009536743)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0087_result.bin",(R.status,R.objval,R.iters));savebin("s0087_x.bin",R.x);savebin("s0087_audit.bin",_rsm_audit_read())
fh=fopen("s0087_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0087_eqmat.bin",E)
savebin("s0087_eqrhs.bin",b)
savebin("s0087_ineqmat.bin",U)
savebin("s0087_ineqrhs.bin",u)
savebin("s0087_cost.bin",c)
end
mata: savebin("s0087_success.bin",_rsm_success_read());savebin("s0087_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578\0.000244140625,7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08\256,0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125\268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578)
b=(7340032\6.67572021484375e-06\7\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0090_result.bin",(R.status,R.objval,R.iters));savebin("s0090_x.bin",R.x);savebin("s0090_audit.bin",_rsm_audit_read())
fh=fopen("s0090_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0090_eqmat.bin",E)
savebin("s0090_eqrhs.bin",b)
savebin("s0090_ineqmat.bin",U)
savebin("s0090_ineqrhs.bin",u)
savebin("s0090_cost.bin",c)
end
mata: savebin("s0090_success.bin",_rsm_success_read());savebin("s0090_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322\0.000244140625,7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08\262144,8.0000000037252903,4095.9999990463257,49152.000015258789,384,327679.99993896484,32\262144.00024414062,8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0093_result.bin",(R.status,R.objval,R.iters));savebin("s0093_x.bin",R.x);savebin("s0093_audit.bin",_rsm_audit_read())
fh=fopen("s0093_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0093_eqmat.bin",E)
savebin("s0093_eqrhs.bin",b)
savebin("s0093_ineqmat.bin",U)
savebin("s0093_ineqrhs.bin",u)
savebin("s0093_cost.bin",c)
end
mata: savebin("s0093_success.bin",_rsm_success_read());savebin("s0093_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,8,1048575.9997558594,192.00000005960464,383.99999994039536,83886080,0.12500000011641532\9.5367431640625e-07,7.4505805969238281e-09,0.0009765625,1.7881393432617188e-07,3.5762786865234375e-07,0.078125,1.1641532182693481e-10\1024,8.0000000037252903,1048575.9997558594,192.00000005960464,384,83886079.984375,0.125\1024.0000009536743,8,1048575.9997558594,192.00000005960464,383.99999994039536,83886080,0.12500000011641532)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0096_result.bin",(R.status,R.objval,R.iters));savebin("s0096_x.bin",R.x);savebin("s0096_audit.bin",_rsm_audit_read())
fh=fopen("s0096_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0096_eqmat.bin",E)
savebin("s0096_eqrhs.bin",b)
savebin("s0096_ineqmat.bin",U)
savebin("s0096_ineqrhs.bin",u)
savebin("s0096_cost.bin",c)
end
mata: savebin("s0096_success.bin",_rsm_success_read());savebin("s0096_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192\0.000244140625,2.9802322387695312e-08,0.00030517578125,3.5762786865234375e-07,4.57763671875e-05,3.814697265625e-06,7.4505805969238281e-09\256,0.03125,319.99999994039536,0.375,48.000000014901161,3.9999999990686774,0.0078125000036379788\268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192)
b=(7340032\6.67572021484375e-06\7\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0099_result.bin",(R.status,R.objval,R.iters));savebin("s0099_x.bin",R.x);savebin("s0099_audit.bin",_rsm_audit_read())
fh=fopen("s0099_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0099_eqmat.bin",E)
savebin("s0099_eqrhs.bin",b)
savebin("s0099_ineqmat.bin",U)
savebin("s0099_ineqrhs.bin",u)
savebin("s0099_cost.bin",c)
end
mata: savebin("s0099_success.bin",_rsm_success_read());savebin("s0099_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,32.000000029802322,327680,383.99999994039536,49152.000015258789,4095.9999990463257,8\0.000244140625,2.9802322387695312e-08,0.00030517578125,3.5762786865234375e-07,4.57763671875e-05,3.814697265625e-06,7.4505805969238281e-09\262144,32,327679.99993896484,384,49152.000015258789,4095.9999990463257,8.0000000037252903\262144.00024414062,32.000000029802322,327680,383.99999994039536,49152.000015258789,4095.9999990463257,8)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0102_result.bin",(R.status,R.objval,R.iters));savebin("s0102_x.bin",R.x);savebin("s0102_audit.bin",_rsm_audit_read())
fh=fopen("s0102_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0102_eqmat.bin",E)
savebin("s0102_eqrhs.bin",b)
savebin("s0102_ineqmat.bin",U)
savebin("s0102_ineqrhs.bin",u)
savebin("s0102_cost.bin",c)
end
mata: savebin("s0102_success.bin",_rsm_success_read());savebin("s0102_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,32.000000029802322,83886080,1.4999999997671694,49152.000015258789,1048575.9997558594,0.03125\9.5367431640625e-07,2.9802322387695312e-08,0.078125,1.3969838619232178e-09,4.57763671875e-05,0.0009765625,2.9103830456733704e-11\1024,32,83886079.984375,1.5,49152.000015258789,1048575.9997558594,0.031250000014551915\1024.0000009536743,32.000000029802322,83886080,1.4999999997671694,49152.000015258789,1048575.9997558594,0.03125)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0105_result.bin",(R.status,R.objval,R.iters));savebin("s0105_x.bin",R.x);savebin("s0105_audit.bin",_rsm_audit_read())
fh=fopen("s0105_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0105_eqmat.bin",E)
savebin("s0105_eqrhs.bin",b)
savebin("s0105_ineqmat.bin",U)
savebin("s0105_ineqrhs.bin",u)
savebin("s0105_cost.bin",c)
end
mata: savebin("s0105_success.bin",_rsm_success_read());savebin("s0105_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25\8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125,256)
b=(7340032\7340032\6.67572021484375e-06\7)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0108_result.bin",(R.status,R.objval,R.iters));savebin("s0108_x.bin",R.x);savebin("s0108_audit.bin",_rsm_audit_read())
fh=fopen("s0108_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0108_eqmat.bin",E)
savebin("s0108_eqrhs.bin",b)
savebin("s0108_ineqmat.bin",U)
savebin("s0108_ineqrhs.bin",u)
savebin("s0108_cost.bin",c)
end
mata: savebin("s0108_success.bin",_rsm_success_read());savebin("s0108_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322,262144.00024414062\8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25\7.62939453125e-06,0.00390625,0.046875,0.0003662109375,0.3125,3.0517578125e-05,0.25\7.6293945348027137e-06,0.0039062499990905053,0.046875000014551915,0.0003662109375,0.31249999994179234,3.0517578125e-05,0.25)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0111_result.bin",(R.status,R.objval,R.iters));savebin("s0111_x.bin",R.x);savebin("s0111_audit.bin",_rsm_audit_read())
fh=fopen("s0111_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0111_eqmat.bin",E)
savebin("s0111_eqrhs.bin",b)
savebin("s0111_ineqmat.bin",U)
savebin("s0111_ineqrhs.bin",u)
savebin("s0111_cost.bin",c)
end
mata: savebin("s0111_success.bin",_rsm_success_read());savebin("s0111_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999990463257,12582912.00390625,1.4999999997671694,327680,8192.0000076293945,1024.0000009536743\32,4194303.9990234375,12884901892,1535.9999997615814,335544320,8388608.0078125,1048576.0009765625\2.9802322387695312e-08,0.00390625,12,1.430511474609375e-06,0.3125,0.0078125,0.0009765625\2.98023224015731e-08,0.0039062499990905053,12.00000000372529,1.430511474609375e-06,0.31249999994179234,0.0078125,0.0009765625)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0114_result.bin",(R.status,R.objval,R.iters));savebin("s0114_x.bin",R.x);savebin("s0114_audit.bin",_rsm_audit_read())
fh=fopen("s0114_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0114_eqmat.bin",E)
savebin("s0114_eqrhs.bin",b)
savebin("s0114_ineqmat.bin",U)
savebin("s0114_ineqrhs.bin",u)
savebin("s0114_cost.bin",c)
end
mata: savebin("s0114_success.bin",_rsm_success_read());savebin("s0114_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578\268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578\0.000244140625,7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08\256,0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125)
b=(7340032\7340032\6.67572021484375e-06\7)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0117_result.bin",(R.status,R.objval,R.iters));savebin("s0117_x.bin",R.x);savebin("s0117_audit.bin",_rsm_audit_read())
fh=fopen("s0117_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0117_eqmat.bin",E)
savebin("s0117_eqrhs.bin",b)
savebin("s0117_ineqmat.bin",U)
savebin("s0117_ineqrhs.bin",u)
savebin("s0117_cost.bin",c)
end
mata: savebin("s0117_success.bin",_rsm_success_read());savebin("s0117_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322\268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578\0.25,7.62939453125e-06,0.00390625,0.046875,0.0003662109375,0.3125,3.0517578125e-05\0.25,7.6293945348027137e-06,0.0039062499990905053,0.046875000014551915,0.0003662109375,0.31249999994179234,3.0517578125e-05)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0120_result.bin",(R.status,R.objval,R.iters));savebin("s0120_x.bin",R.x);savebin("s0120_audit.bin",_rsm_audit_read())
fh=fopen("s0120_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0120_eqmat.bin",E)
savebin("s0120_eqrhs.bin",b)
savebin("s0120_ineqmat.bin",U)
savebin("s0120_ineqrhs.bin",u)
savebin("s0120_cost.bin",c)
end
mata: savebin("s0120_success.bin",_rsm_success_read());savebin("s0120_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,8,1048575.9997558594,192.00000005960464,383.99999994039536,83886080,0.12500000011641532\1048576.0009765625,8192,1073741823.75,196608.00006103516,393215.99993896484,85899345920,128.00000011920929\0.0009765625,7.62939453125e-06,1,0.00018310546875,0.0003662109375,80,1.1920928955078125e-07\0.0009765625,7.6293945348027137e-06,0.99999999976716936,0.00018310546880684342,0.0003662109375,79.999999985098839,1.1920928955078125e-07)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0123_result.bin",(R.status,R.objval,R.iters));savebin("s0123_x.bin",R.x);savebin("s0123_audit.bin",_rsm_audit_read())
fh=fopen("s0123_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0123_eqmat.bin",E)
savebin("s0123_eqrhs.bin",b)
savebin("s0123_ineqmat.bin",U)
savebin("s0123_ineqrhs.bin",u)
savebin("s0123_cost.bin",c)
end
mata: savebin("s0123_success.bin",_rsm_success_read());savebin("s0123_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192\268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192\0.000244140625,2.9802322387695312e-08,0.00030517578125,3.5762786865234375e-07,4.57763671875e-05,3.814697265625e-06,7.4505805969238281e-09\256,0.03125,319.99999994039536,0.375,48.000000014901161,3.9999999990686774,0.0078125000036379788)
b=(7340032\7340032\6.67572021484375e-06\7)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0126_result.bin",(R.status,R.objval,R.iters));savebin("s0126_x.bin",R.x);savebin("s0126_audit.bin",_rsm_audit_read())
fh=fopen("s0126_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0126_eqmat.bin",E)
savebin("s0126_eqrhs.bin",b)
savebin("s0126_ineqmat.bin",U)
savebin("s0126_ineqrhs.bin",u)
savebin("s0126_cost.bin",c)
end
mata: savebin("s0126_success.bin",_rsm_success_read());savebin("s0126_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,32.000000029802322,327680,383.99999994039536,49152.000015258789,4095.9999990463257,8\268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192\0.25,3.0517578125e-05,0.3125,0.0003662109375,0.046875,0.00390625,7.62939453125e-06\0.25,3.0517578125e-05,0.31249999994179234,0.0003662109375,0.046875000014551915,0.0039062499990905053,7.6293945348027137e-06)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0129_result.bin",(R.status,R.objval,R.iters));savebin("s0129_x.bin",R.x);savebin("s0129_audit.bin",_rsm_audit_read())
fh=fopen("s0129_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0129_eqmat.bin",E)
savebin("s0129_eqrhs.bin",b)
savebin("s0129_ineqmat.bin",U)
savebin("s0129_ineqrhs.bin",u)
savebin("s0129_cost.bin",c)
end
mata: savebin("s0129_success.bin",_rsm_success_read());savebin("s0129_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,32.000000029802322,83886080,1.4999999997671694,49152.000015258789,1048575.9997558594,0.03125\1048576.0009765625,32768.000030517578,85899345920,1535.9999997615814,50331648.015625,1073741823.75,32\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\0.0009765625,3.0517578125e-05,79.999999985098839,1.430511474609375e-06,0.046875000014551915,0.99999999976716936,2.98023224015731e-08)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0132_result.bin",(R.status,R.objval,R.iters));savebin("s0132_x.bin",R.x);savebin("s0132_audit.bin",_rsm_audit_read())
fh=fopen("s0132_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0132_eqmat.bin",E)
savebin("s0132_eqrhs.bin",b)
savebin("s0132_ineqmat.bin",U)
savebin("s0132_ineqrhs.bin",u)
savebin("s0132_cost.bin",c)
end
mata: savebin("s0132_success.bin",_rsm_success_read());savebin("s0132_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25\0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125,256\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578,268435456.25)
b=(7340032\7\6.67572021484375e-06\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0135_result.bin",(R.status,R.objval,R.iters));savebin("s0135_x.bin",R.x);savebin("s0135_audit.bin",_rsm_audit_read())
fh=fopen("s0135_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0135_eqmat.bin",E)
savebin("s0135_eqrhs.bin",b)
savebin("s0135_ineqmat.bin",U)
savebin("s0135_ineqrhs.bin",u)
savebin("s0135_cost.bin",c)
end
mata: savebin("s0135_success.bin",_rsm_success_read());savebin("s0135_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322,262144.00024414062\0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125,256\7.62939453125e-06,0.00390625,0.046875,0.0003662109375,0.3125,3.0517578125e-05,0.25\8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322,262144.00024414062)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0138_result.bin",(R.status,R.objval,R.iters));savebin("s0138_x.bin",R.x);savebin("s0138_audit.bin",_rsm_audit_read())
fh=fopen("s0138_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0138_eqmat.bin",E)
savebin("s0138_eqrhs.bin",b)
savebin("s0138_ineqmat.bin",U)
savebin("s0138_ineqrhs.bin",u)
savebin("s0138_cost.bin",c)
end
mata: savebin("s0138_success.bin",_rsm_success_read());savebin("s0138_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999990463257,12582912.00390625,1.4999999997671694,327680,8192.0000076293945,1024.0000009536743\3.0517578139210855e-05,3.9999999990686774,12288.000003814697,0.00146484375,319.99999994039536,8,1\2.9802322387695312e-08,0.00390625,12,1.430511474609375e-06,0.3125,0.0078125,0.0009765625\0.03125,4095.9999990463257,12582912.00390625,1.4999999997671694,327680,8192.0000076293945,1024.0000009536743)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0141_result.bin",(R.status,R.objval,R.iters));savebin("s0141_x.bin",R.x);savebin("s0141_audit.bin",_rsm_audit_read())
fh=fopen("s0141_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0141_eqmat.bin",E)
savebin("s0141_eqrhs.bin",b)
savebin("s0141_ineqmat.bin",U)
savebin("s0141_ineqrhs.bin",u)
savebin("s0141_cost.bin",c)
end
mata: savebin("s0141_success.bin",_rsm_success_read());savebin("s0141_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578\256,0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125\0.000244140625,7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08\268435456.25,8192,4194303.9990234375,50331648.015625,393215.99993896484,335544320,32768.000030517578)
b=(7340032\7\6.67572021484375e-06\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0144_result.bin",(R.status,R.objval,R.iters));savebin("s0144_x.bin",R.x);savebin("s0144_audit.bin",_rsm_audit_read())
fh=fopen("s0144_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0144_eqmat.bin",E)
savebin("s0144_eqrhs.bin",b)
savebin("s0144_ineqmat.bin",U)
savebin("s0144_ineqrhs.bin",u)
savebin("s0144_cost.bin",c)
end
mata: savebin("s0144_success.bin",_rsm_success_read());savebin("s0144_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322\256,0.0078125000036379788,3.9999999990686774,48.000000014901161,0.375,319.99999994039536,0.03125\0.25,7.62939453125e-06,0.00390625,0.046875,0.0003662109375,0.3125,3.0517578125e-05\262144.00024414062,8,4095.9999990463257,49152.000015258789,383.99999994039536,327680,32.000000029802322)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0147_result.bin",(R.status,R.objval,R.iters));savebin("s0147_x.bin",R.x);savebin("s0147_audit.bin",_rsm_audit_read())
fh=fopen("s0147_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0147_eqmat.bin",E)
savebin("s0147_eqrhs.bin",b)
savebin("s0147_ineqmat.bin",U)
savebin("s0147_ineqrhs.bin",u)
savebin("s0147_cost.bin",c)
end
mata: savebin("s0147_success.bin",_rsm_success_read());savebin("s0147_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,8,1048575.9997558594,192.00000005960464,383.99999994039536,83886080,0.12500000011641532\1,0.0078125000036379788,1023.9999997615814,0.18750000005820766,0.375,81919.999984741211,0.0001220703125\0.0009765625,7.62939453125e-06,1,0.00018310546875,0.0003662109375,80,1.1920928955078125e-07\1024.0000009536743,8,1048575.9997558594,192.00000005960464,383.99999994039536,83886080,0.12500000011641532)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0150_result.bin",(R.status,R.objval,R.iters));savebin("s0150_x.bin",R.x);savebin("s0150_audit.bin",_rsm_audit_read())
fh=fopen("s0150_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0150_eqmat.bin",E)
savebin("s0150_eqrhs.bin",b)
savebin("s0150_ineqmat.bin",U)
savebin("s0150_ineqrhs.bin",u)
savebin("s0150_cost.bin",c)
end
mata: savebin("s0150_success.bin",_rsm_success_read());savebin("s0150_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192\256,0.03125,319.99999994039536,0.375,48.000000014901161,3.9999999990686774,0.0078125000036379788\0.000244140625,2.9802322387695312e-08,0.00030517578125,3.5762786865234375e-07,4.57763671875e-05,3.814697265625e-06,7.4505805969238281e-09\268435456.25,32768.000030517578,335544320,393215.99993896484,50331648.015625,4194303.9990234375,8192)
b=(7340032\7\6.67572021484375e-06\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0153_result.bin",(R.status,R.objval,R.iters));savebin("s0153_x.bin",R.x);savebin("s0153_audit.bin",_rsm_audit_read())
fh=fopen("s0153_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0153_eqmat.bin",E)
savebin("s0153_eqrhs.bin",b)
savebin("s0153_ineqmat.bin",U)
savebin("s0153_ineqrhs.bin",u)
savebin("s0153_cost.bin",c)
end
mata: savebin("s0153_success.bin",_rsm_success_read());savebin("s0153_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(262144.00024414062,32.000000029802322,327680,383.99999994039536,49152.000015258789,4095.9999990463257,8\256,0.03125,319.99999994039536,0.375,48.000000014901161,3.9999999990686774,0.0078125000036379788\0.25,3.0517578125e-05,0.3125,0.0003662109375,0.046875,0.00390625,7.62939453125e-06\262144.00024414062,32.000000029802322,327680,383.99999994039536,49152.000015258789,4095.9999990463257,8)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.046875,320,0.375,64,3,0.0078125)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0156_result.bin",(R.status,R.objval,R.iters));savebin("s0156_x.bin",R.x);savebin("s0156_audit.bin",_rsm_audit_read())
fh=fopen("s0156_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0156_eqmat.bin",E)
savebin("s0156_eqrhs.bin",b)
savebin("s0156_ineqmat.bin",U)
savebin("s0156_ineqrhs.bin",u)
savebin("s0156_cost.bin",c)
end
mata: savebin("s0156_success.bin",_rsm_success_read());savebin("s0156_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000009536743,32.000000029802322,83886080,1.4999999997671694,49152.000015258789,1048575.9997558594,0.03125\1,0.03125,81919.999984741211,0.00146484375,48.000000014901161,1023.9999997615814,3.0517578139210855e-05\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\1024.0000009536743,32.000000029802322,83886080,1.4999999997671694,49152.000015258789,1048575.9997558594,0.03125)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0159_result.bin",(R.status,R.objval,R.iters));savebin("s0159_x.bin",R.x);savebin("s0159_audit.bin",_rsm_audit_read())
fh=fopen("s0159_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0159_eqmat.bin",E)
savebin("s0159_eqrhs.bin",b)
savebin("s0159_ineqmat.bin",U)
savebin("s0159_ineqrhs.bin",u)
savebin("s0159_cost.bin",c)
end
mata: savebin("s0159_success.bin",_rsm_success_read());savebin("s0159_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8192,4194303.9997558594,50331648.00390625,393215.99998474121,335544320,32768.000007629395,268435456.0625\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\0.0078125000009094947,3.9999999997671694,48.00000000372529,0.375,319.99999998509884,0.03125,256\8192,4194303.9997558594,50331648.00390625,393215.99998474121,335544320,32768.000007629395,268435456.0625)
b=(7340032\6.67572021484375e-06\7\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0324_result.bin",(R.status,R.objval,R.iters));savebin("s0324_x.bin",R.x);savebin("s0324_audit.bin",_rsm_audit_read())
fh=fopen("s0324_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0324_eqmat.bin",E)
savebin("s0324_eqrhs.bin",b)
savebin("s0324_ineqmat.bin",U)
savebin("s0324_ineqrhs.bin",u)
savebin("s0324_cost.bin",c)
end
mata: savebin("s0324_success.bin",_rsm_success_read());savebin("s0324_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8,4095.9999997615814,49152.000003814697,383.99999998509884,327680,32.000000007450581,262144.00006103516\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\8.0000000009313226,4095.9999997615814,49152.000003814697,384,327679.99998474121,32,262144\8,4095.9999997615814,49152.000003814697,383.99999998509884,327680,32.000000007450581,262144.00006103516)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0327_result.bin",(R.status,R.objval,R.iters));savebin("s0327_x.bin",R.x);savebin("s0327_audit.bin",_rsm_audit_read())
fh=fopen("s0327_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0327_eqmat.bin",E)
savebin("s0327_eqrhs.bin",b)
savebin("s0327_ineqmat.bin",U)
savebin("s0327_ineqrhs.bin",u)
savebin("s0327_cost.bin",c)
end
mata: savebin("s0327_success.bin",_rsm_success_read());savebin("s0327_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999997615814,12582912.000976562,1.4999999999417923,327680,8192.0000019073486,1024.0000002384186\2.9103830456733704e-11,3.814697265625e-06,0.01171875,1.3969838619232178e-09,0.00030517578125,7.62939453125e-06,9.5367431640625e-07\0.031250000003637979,4095.9999997615814,12582912.000976562,1.5,327679.99998474121,8192,1024\0.03125,4095.9999997615814,12582912.000976562,1.4999999999417923,327680,8192.0000019073486,1024.0000002384186)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0330_result.bin",(R.status,R.objval,R.iters));savebin("s0330_x.bin",R.x);savebin("s0330_audit.bin",_rsm_audit_read())
fh=fopen("s0330_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0330_eqmat.bin",E)
savebin("s0330_eqrhs.bin",b)
savebin("s0330_ineqmat.bin",U)
savebin("s0330_ineqrhs.bin",u)
savebin("s0330_cost.bin",c)
end
mata: savebin("s0330_success.bin",_rsm_success_read());savebin("s0330_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,8,1048575.9999389648,192.00000001490116,383.99999998509884,83886080,0.12500000002910383\9.5367431640625e-07,7.4505805969238281e-09,0.0009765625,1.7881393432617188e-07,3.5762786865234375e-07,0.078125,1.1641532182693481e-10\1024,8.0000000009313226,1048575.9999389648,192.00000001490116,384,83886079.99609375,0.125\1024.0000002384186,8,1048575.9999389648,192.00000001490116,383.99999998509884,83886080,0.12500000002910383)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0339_result.bin",(R.status,R.objval,R.iters));savebin("s0339_x.bin",R.x);savebin("s0339_audit.bin",_rsm_audit_read())
fh=fopen("s0339_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0339_eqmat.bin",E)
savebin("s0339_eqrhs.bin",b)
savebin("s0339_ineqmat.bin",U)
savebin("s0339_ineqrhs.bin",u)
savebin("s0339_cost.bin",c)
end
mata: savebin("s0339_success.bin",_rsm_success_read());savebin("s0339_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\9.5367431640625e-07,2.9802322387695312e-08,0.078125,1.3969838619232178e-09,4.57763671875e-05,0.0009765625,2.9103830456733704e-11\1024,32,83886079.99609375,1.5,49152.000003814697,1048575.9999389648,0.031250000003637979\1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0348_result.bin",(R.status,R.objval,R.iters));savebin("s0348_x.bin",R.x);savebin("s0348_audit.bin",_rsm_audit_read())
fh=fopen("s0348_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0348_eqmat.bin",E)
savebin("s0348_eqrhs.bin",b)
savebin("s0348_ineqmat.bin",U)
savebin("s0348_ineqrhs.bin",u)
savebin("s0348_cost.bin",c)
end
mata: savebin("s0348_success.bin",_rsm_success_read());savebin("s0348_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\9.5367431640625e-07,2.9802322387695312e-08,0.078125,1.3969838619232178e-09,4.57763671875e-05,0.0009765625,2.9103830456733704e-11\1024,32,83886079.99609375,1.5,49152.000003814697,1048575.9999389648,0.031250000003637979\1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125)
b=(7168\6.67572021484375e-06\7168\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-10,"")
savebin("s0349_result.bin",(R.status,R.objval,R.iters));savebin("s0349_x.bin",R.x);savebin("s0349_audit.bin",_rsm_audit_read())
fh=fopen("s0349_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0349_eqmat.bin",E)
savebin("s0349_eqrhs.bin",b)
savebin("s0349_ineqmat.bin",U)
savebin("s0349_ineqrhs.bin",u)
savebin("s0349_cost.bin",c)
end
mata: savebin("s0349_success.bin",_rsm_success_read());savebin("s0349_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8192,4194303.9997558594,50331648.00390625,393215.99998474121,335544320,32768.000007629395,268435456.0625\8192,4194303.9997558594,50331648.00390625,393215.99998474121,335544320,32768.000007629395,268435456.0625\7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08,0.000244140625\0.0078125000009094947,3.9999999997671694,48.00000000372529,0.375,319.99999998509884,0.03125,256)
b=(7340032\7340032\6.67572021484375e-06\7)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0351_result.bin",(R.status,R.objval,R.iters));savebin("s0351_x.bin",R.x);savebin("s0351_audit.bin",_rsm_audit_read())
fh=fopen("s0351_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0351_eqmat.bin",E)
savebin("s0351_eqrhs.bin",b)
savebin("s0351_ineqmat.bin",U)
savebin("s0351_ineqrhs.bin",u)
savebin("s0351_cost.bin",c)
end
mata: savebin("s0351_success.bin",_rsm_success_read());savebin("s0351_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(8,4095.9999997615814,49152.000003814697,383.99999998509884,327680,32.000000007450581,262144.00006103516\8192,4194303.9997558594,50331648.00390625,393215.99998474121,335544320,32768.000007629395,268435456.0625\7.62939453125e-06,0.00390625,0.046875,0.0003662109375,0.3125,3.0517578125e-05,0.25\7.6293945321381784e-06,0.0039062499997726263,0.046875000003637979,0.0003662109375,0.31249999998544808,3.0517578125e-05,0.25)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(0.0078125,3,64,0.375,320,0.046875,256)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0354_result.bin",(R.status,R.objval,R.iters));savebin("s0354_x.bin",R.x);savebin("s0354_audit.bin",_rsm_audit_read())
fh=fopen("s0354_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0354_eqmat.bin",E)
savebin("s0354_eqrhs.bin",b)
savebin("s0354_ineqmat.bin",U)
savebin("s0354_ineqrhs.bin",u)
savebin("s0354_cost.bin",c)
end
mata: savebin("s0354_success.bin",_rsm_success_read());savebin("s0354_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999997615814,12582912.000976562,1.4999999999417923,327680,8192.0000019073486,1024.0000002384186\32,4194303.9997558594,12884901889,1535.9999999403954,335544320,8388608.001953125,1048576.0002441406\2.9802322387695312e-08,0.00390625,12,1.430511474609375e-06,0.3125,0.0078125,0.0009765625\2.9802322391164759e-08,0.0039062499997726263,12.000000000931323,1.430511474609375e-06,0.31249999998544808,0.0078125,0.0009765625)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0357_result.bin",(R.status,R.objval,R.iters));savebin("s0357_x.bin",R.x);savebin("s0357_audit.bin",_rsm_audit_read())
fh=fopen("s0357_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0357_eqmat.bin",E)
savebin("s0357_eqrhs.bin",b)
savebin("s0357_ineqmat.bin",U)
savebin("s0357_ineqrhs.bin",u)
savebin("s0357_cost.bin",c)
end
mata: savebin("s0357_success.bin",_rsm_success_read());savebin("s0357_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,8,1048575.9999389648,192.00000001490116,383.99999998509884,83886080,0.12500000002910383\1048576.0002441406,8192,1073741823.9375,196608.00001525879,393215.99998474121,85899345920,128.00000002980232\0.0009765625,7.62939453125e-06,1,0.00018310546875,0.0003662109375,80,1.1920928955078125e-07\0.0009765625,7.6293945321381784e-06,0.99999999994179234,0.00018310546876421085,0.0003662109375,79.99999999627471,1.1920928955078125e-07)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0366_result.bin",(R.status,R.objval,R.iters));savebin("s0366_x.bin",R.x);savebin("s0366_audit.bin",_rsm_audit_read())
fh=fopen("s0366_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0366_eqmat.bin",E)
savebin("s0366_eqrhs.bin",b)
savebin("s0366_ineqmat.bin",U)
savebin("s0366_ineqrhs.bin",u)
savebin("s0366_cost.bin",c)
end
mata: savebin("s0366_success.bin",_rsm_success_read());savebin("s0366_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\1048576.0002441406,32768.000007629395,85899345920,1535.9999999403954,50331648.00390625,1073741823.9375,32\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\0.0009765625,3.0517578125e-05,79.99999999627471,1.430511474609375e-06,0.046875000003637979,0.99999999994179234,2.9802322391164759e-08)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0375_result.bin",(R.status,R.objval,R.iters));savebin("s0375_x.bin",R.x);savebin("s0375_audit.bin",_rsm_audit_read())
fh=fopen("s0375_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0375_eqmat.bin",E)
savebin("s0375_eqrhs.bin",b)
savebin("s0375_ineqmat.bin",U)
savebin("s0375_ineqrhs.bin",u)
savebin("s0375_cost.bin",c)
end
mata: savebin("s0375_success.bin",_rsm_success_read());savebin("s0375_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\1048576.0002441406,32768.000007629395,85899345920,1535.9999999403954,50331648.00390625,1073741823.9375,32\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\0.0009765625,3.0517578125e-05,79.99999999627471,1.430511474609375e-06,0.046875000003637979,0.99999999994179234,2.9802322391164759e-08)
b=(7168\7340032\0.0068359375\0.0068359375)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-10,"")
savebin("s0376_result.bin",(R.status,R.objval,R.iters));savebin("s0376_x.bin",R.x);savebin("s0376_audit.bin",_rsm_audit_read())
fh=fopen("s0376_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0376_eqmat.bin",E)
savebin("s0376_eqrhs.bin",b)
savebin("s0376_ineqmat.bin",U)
savebin("s0376_ineqrhs.bin",u)
savebin("s0376_cost.bin",c)
end
mata: savebin("s0376_success.bin",_rsm_success_read());savebin("s0376_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.03125,4095.9999997615814,12582912.000976562,1.4999999999417923,327680,8192.0000019073486,1024.0000002384186\3.0517578128552714e-05,3.9999999997671694,12288.000000953674,0.00146484375,319.99999998509884,8,1\2.9802322387695312e-08,0.00390625,12,1.430511474609375e-06,0.3125,0.0078125,0.0009765625\0.03125,4095.9999997615814,12582912.000976562,1.4999999999417923,327680,8192.0000019073486,1024.0000002384186)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(3.0517578125e-05,3,16384,0.00146484375,320,12,1)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0384_result.bin",(R.status,R.objval,R.iters));savebin("s0384_x.bin",R.x);savebin("s0384_audit.bin",_rsm_audit_read())
fh=fopen("s0384_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0384_eqmat.bin",E)
savebin("s0384_eqrhs.bin",b)
savebin("s0384_ineqmat.bin",U)
savebin("s0384_ineqrhs.bin",u)
savebin("s0384_cost.bin",c)
end
mata: savebin("s0384_success.bin",_rsm_success_read());savebin("s0384_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,8,1048575.9999389648,192.00000001490116,383.99999998509884,83886080,0.12500000002910383\1,0.0078125000009094947,1023.9999999403954,0.18750000001455192,0.375,81919.999996185303,0.0001220703125\0.0009765625,7.62939453125e-06,1,0.00018310546875,0.0003662109375,80,1.1920928955078125e-07\1024.0000002384186,8,1048575.9999389648,192.00000001490116,383.99999998509884,83886080,0.12500000002910383)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.0078125,768,0.25,0.375,81920,0.00018310546875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0393_result.bin",(R.status,R.objval,R.iters));savebin("s0393_x.bin",R.x);savebin("s0393_audit.bin",_rsm_audit_read())
fh=fopen("s0393_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0393_eqmat.bin",E)
savebin("s0393_eqrhs.bin",b)
savebin("s0393_ineqmat.bin",U)
savebin("s0393_ineqrhs.bin",u)
savebin("s0393_cost.bin",c)
end
mata: savebin("s0393_success.bin",_rsm_success_read());savebin("s0393_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\1,0.03125,81919.999996185303,0.00146484375,48.00000000372529,1023.9999999403954,3.0517578128552714e-05\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125)
b=(7168\7\0.0068359375\7168)
U=J(0,7,0)
u=J(0,1,0)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("s0402_result.bin",(R.status,R.objval,R.iters));savebin("s0402_x.bin",R.x);savebin("s0402_audit.bin",_rsm_audit_read())
fh=fopen("s0402_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("s0402_eqmat.bin",E)
savebin("s0402_eqrhs.bin",b)
savebin("s0402_ineqmat.bin",U)
savebin("s0402_ineqrhs.bin",u)
savebin("s0402_cost.bin",c)
end
mata: savebin("s0402_success.bin",_rsm_success_read());savebin("s0402_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(0.000244140625,7.4505805969238281e-09,3.814697265625e-06,4.57763671875e-05,3.5762786865234375e-07,0.00030517578125,2.9802322387695312e-08\256,0.0078125001164153218,3.9999999701976776,48.000000476837158,0.375,319.99999809265137,0.03125\268435464,8192,4194303.96875,50331648.5,393215.998046875,335544320,32768.0009765625)
b=(6.67572021484375e-06\7\7340032)
U=J(0,7,0)
u=J(0,1,0)
c=(256,0.0078125,3,64,0.375,320,0.046875)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-08,"")
savebin("e015_result.bin",(R.status,R.objval,R.iters));savebin("e015_x.bin",R.x);savebin("e015_audit.bin",_rsm_audit_read())
fh=fopen("e015_message.txt","w");fput(fh,R.msg);fclose(fh)
savebin("e015_eqmat.bin",E)
savebin("e015_eqrhs.bin",b)
savebin("e015_ineqmat.bin",U)
savebin("e015_ineqrhs.bin",u)
savebin("e015_cost.bin",c)
end
mata: savebin("e015_success.bin",_rsm_success_read());savebin("e015_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
1 2 3 "" .
1 2 4 "=" 7
1 2.03125 3.96875 "=" 7
1.03125 2 3.96875 "=" 7
end
mata: savebin("w191_none_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("w191_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w191_none_x.bin",st_matrix("r(x)"));savebin("w191_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w191_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w191_none_success.bin",_rsm_success_read());savebin("w191_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
1 2 3 "" .
1 2 4 "=" 7
1 2.03125 3.96875 "=" 7
1.03125 2 3.96875 "=" 7
end
mata: savebin("w191_row_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("w191_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w191_row_x.bin",st_matrix("r(x)"));savebin("w191_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w191_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w191_row_success.bin",_rsm_success_read());savebin("w191_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
1 2 3 "" .
1 2 4 "=" 7
1 2.0000000298023224 3.9999999701976776 "=" 7
1.0000000298023224 2 3.9999999701976776 "=" 7
end
mata: savebin("w264_none_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(none) tolerance(1e-10)
mata: savebin("w264_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w264_none_x.bin",st_matrix("r(x)"));savebin("w264_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w264_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w264_none_success.bin",_rsm_success_read());savebin("w264_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
1 2 3 "" .
1 2 4 "=" 7
1 2.0000000298023224 3.9999999701976776 "=" 7
1.0000000298023224 2 3.9999999701976776 "=" 7
end
mata: savebin("w264_row_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(row) tolerance(1e-10)
mata: savebin("w264_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w264_row_x.bin",st_matrix("r(x)"));savebin("w264_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w264_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w264_row_success.bin",_rsm_success_read());savebin("w264_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
2 3 1 "" .
2 3.9999999701976776 1.0000000298023224 "=" 7
2 4 1 "=" 7
2.0000000298023224 3.9999999701976776 1 "=" 7
end
mata: savebin("w267_none_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(none) tolerance(1e-10)
mata: savebin("w267_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w267_none_x.bin",st_matrix("r(x)"));savebin("w267_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w267_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w267_none_success.bin",_rsm_success_read());savebin("w267_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 str2 rel double rhs
2 3 1 "" .
2 3.9999999701976776 1.0000000298023224 "=" 7
2 4 1 "=" 7
2.0000000298023224 3.9999999701976776 1 "=" 7
end
mata: savebin("w267_row_input.bin",st_data(.,tokens("x1 x2 x3 rhs")))
lp x1 x2 x3,min engine(simplex) scaling(row) tolerance(1e-10)
mata: savebin("w267_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("w267_row_x.bin",st_matrix("r(x)"));savebin("w267_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("w267_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("w267_row_success.bin",_rsm_success_read());savebin("w267_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t0_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("e015_t0_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t0_none_x.bin",st_matrix("r(x)"));savebin("e015_t0_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t0_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t0_none_success.bin",_rsm_success_read());savebin("e015_t0_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t0_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("e015_t0_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t0_row_x.bin",st_matrix("r(x)"));savebin("e015_t0_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t0_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t0_row_success.bin",_rsm_success_read());savebin("e015_t0_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t1_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-10)
mata: savebin("e015_t1_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t1_none_x.bin",st_matrix("r(x)"));savebin("e015_t1_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t1_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t1_none_success.bin",_rsm_success_read());savebin("e015_t1_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t1_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-10)
mata: savebin("e015_t1_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t1_row_x.bin",st_matrix("r(x)"));savebin("e015_t1_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t1_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t1_row_success.bin",_rsm_success_read());savebin("e015_t1_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t2_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(9.9999999999999998e-13)
mata: savebin("e015_t2_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t2_none_x.bin",st_matrix("r(x)"));savebin("e015_t2_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t2_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t2_none_success.bin",_rsm_success_read());savebin("e015_t2_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125001164153218 3.9999999701976776 48.000000476837158 0.375 319.99999809265137 0.03125 "=" 7
268435464 8192 4194303.96875 50331648.5 393215.998046875 335544320 32768.0009765625 "=" 7340032
end
mata: savebin("e015_t2_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(9.9999999999999998e-13)
mata: savebin("e015_t2_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("e015_t2_row_x.bin",st_matrix("r(x)"));savebin("e015_t2_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("e015_t2_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("e015_t2_row_success.bin",_rsm_success_read());savebin("e015_t2_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
end
mata: savebin("s0081_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0081_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0081_public_none_x.bin",st_matrix("r(x)"));savebin("s0081_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0081_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0081_public_none_success.bin",_rsm_success_read());savebin("s0081_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
end
mata: savebin("s0081_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0081_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0081_public_row_x.bin",st_matrix("r(x)"));savebin("s0081_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0081_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0081_public_row_success.bin",_rsm_success_read());savebin("s0081_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8.0000000037252903 4095.9999990463257 49152.000015258789 384 327679.99993896484 32 262144 "=" 7168
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
end
mata: savebin("s0084_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0084_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0084_public_none_x.bin",st_matrix("r(x)"));savebin("s0084_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0084_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0084_public_none_success.bin",_rsm_success_read());savebin("s0084_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8.0000000037252903 4095.9999990463257 49152.000015258789 384 327679.99993896484 32 262144 "=" 7168
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
end
mata: savebin("s0084_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0084_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0084_public_row_x.bin",st_matrix("r(x)"));savebin("s0084_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0084_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0084_public_row_success.bin",_rsm_success_read());savebin("s0084_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
2.9103830456733704e-11 3.814697265625e-06 0.01171875 1.3969838619232178e-09 0.00030517578125 7.62939453125e-06 9.5367431640625e-07 "=" 6.67572021484375e-06
0.031250000014551915 4095.9999990463257 12582912.00390625 1.5 327679.99993896484 8192 1024 "=" 7168
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
end
mata: savebin("s0087_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0087_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0087_public_none_x.bin",st_matrix("r(x)"));savebin("s0087_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0087_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0087_public_none_success.bin",_rsm_success_read());savebin("s0087_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
2.9103830456733704e-11 3.814697265625e-06 0.01171875 1.3969838619232178e-09 0.00030517578125 7.62939453125e-06 9.5367431640625e-07 "=" 6.67572021484375e-06
0.031250000014551915 4095.9999990463257 12582912.00390625 1.5 327679.99993896484 8192 1024 "=" 7168
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
end
mata: savebin("s0087_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0087_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0087_public_row_x.bin",st_matrix("r(x)"));savebin("s0087_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0087_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0087_public_row_success.bin",_rsm_success_read());savebin("s0087_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
end
mata: savebin("s0090_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0090_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0090_public_none_x.bin",st_matrix("r(x)"));savebin("s0090_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0090_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0090_public_none_success.bin",_rsm_success_read());savebin("s0090_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
end
mata: savebin("s0090_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0090_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0090_public_row_x.bin",st_matrix("r(x)"));savebin("s0090_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0090_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0090_public_row_success.bin",_rsm_success_read());savebin("s0090_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
262144 8.0000000037252903 4095.9999990463257 49152.000015258789 384 327679.99993896484 32 "=" 7168
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
end
mata: savebin("s0093_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0093_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0093_public_none_x.bin",st_matrix("r(x)"));savebin("s0093_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0093_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0093_public_none_success.bin",_rsm_success_read());savebin("s0093_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
262144 8.0000000037252903 4095.9999990463257 49152.000015258789 384 327679.99993896484 32 "=" 7168
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
end
mata: savebin("s0093_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0093_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0093_public_row_x.bin",st_matrix("r(x)"));savebin("s0093_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0093_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0093_public_row_success.bin",_rsm_success_read());savebin("s0093_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
9.5367431640625e-07 7.4505805969238281e-09 0.0009765625 1.7881393432617188e-07 3.5762786865234375e-07 0.078125 1.1641532182693481e-10 "=" 6.67572021484375e-06
1024 8.0000000037252903 1048575.9997558594 192.00000005960464 384 83886079.984375 0.125 "=" 7168
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
end
mata: savebin("s0096_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0096_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0096_public_none_x.bin",st_matrix("r(x)"));savebin("s0096_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0096_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0096_public_none_success.bin",_rsm_success_read());savebin("s0096_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
9.5367431640625e-07 7.4505805969238281e-09 0.0009765625 1.7881393432617188e-07 3.5762786865234375e-07 0.078125 1.1641532182693481e-10 "=" 6.67572021484375e-06
1024 8.0000000037252903 1048575.9997558594 192.00000005960464 384 83886079.984375 0.125 "=" 7168
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
end
mata: savebin("s0096_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0096_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0096_public_row_x.bin",st_matrix("r(x)"));savebin("s0096_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0096_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0096_public_row_success.bin",_rsm_success_read());savebin("s0096_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
end
mata: savebin("s0099_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0099_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0099_public_none_x.bin",st_matrix("r(x)"));savebin("s0099_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0099_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0099_public_none_success.bin",_rsm_success_read());savebin("s0099_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
end
mata: savebin("s0099_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0099_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0099_public_row_x.bin",st_matrix("r(x)"));savebin("s0099_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0099_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0099_public_row_success.bin",_rsm_success_read());savebin("s0099_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
262144 32 327679.99993896484 384 49152.000015258789 4095.9999990463257 8.0000000037252903 "=" 7168
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
end
mata: savebin("s0102_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0102_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0102_public_none_x.bin",st_matrix("r(x)"));savebin("s0102_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0102_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0102_public_none_success.bin",_rsm_success_read());savebin("s0102_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
262144 32 327679.99993896484 384 49152.000015258789 4095.9999990463257 8.0000000037252903 "=" 7168
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
end
mata: savebin("s0102_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0102_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0102_public_row_x.bin",st_matrix("r(x)"));savebin("s0102_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0102_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0102_public_row_success.bin",_rsm_success_read());savebin("s0102_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.984375 1.5 49152.000015258789 1048575.9997558594 0.031250000014551915 "=" 7168
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
end
mata: savebin("s0105_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0105_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0105_public_none_x.bin",st_matrix("r(x)"));savebin("s0105_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0105_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0105_public_none_success.bin",_rsm_success_read());savebin("s0105_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.984375 1.5 49152.000015258789 1048575.9997558594 0.031250000014551915 "=" 7168
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
end
mata: savebin("s0105_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0105_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0105_public_row_x.bin",st_matrix("r(x)"));savebin("s0105_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0105_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0105_public_row_success.bin",_rsm_success_read());savebin("s0105_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
end
mata: savebin("s0108_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0108_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0108_public_none_x.bin",st_matrix("r(x)"));savebin("s0108_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0108_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0108_public_none_success.bin",_rsm_success_read());savebin("s0108_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
end
mata: savebin("s0108_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0108_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0108_public_row_x.bin",st_matrix("r(x)"));savebin("s0108_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0108_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0108_public_row_success.bin",_rsm_success_read());savebin("s0108_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
7.6293945348027137e-06 0.0039062499990905053 0.046875000014551915 0.0003662109375 0.31249999994179234 3.0517578125e-05 0.25 "=" 0.0068359375
end
mata: savebin("s0111_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0111_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0111_public_none_x.bin",st_matrix("r(x)"));savebin("s0111_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0111_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0111_public_none_success.bin",_rsm_success_read());savebin("s0111_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
7.6293945348027137e-06 0.0039062499990905053 0.046875000014551915 0.0003662109375 0.31249999994179234 3.0517578125e-05 0.25 "=" 0.0068359375
end
mata: savebin("s0111_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0111_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0111_public_row_x.bin",st_matrix("r(x)"));savebin("s0111_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0111_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0111_public_row_success.bin",_rsm_success_read());savebin("s0111_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
32 4194303.9990234375 12884901892 1535.9999997615814 335544320 8388608.0078125 1048576.0009765625 "=" 7340032
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
2.98023224015731e-08 0.0039062499990905053 12.00000000372529 1.430511474609375e-06 0.31249999994179234 0.0078125 0.0009765625 "=" 0.0068359375
end
mata: savebin("s0114_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0114_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0114_public_none_x.bin",st_matrix("r(x)"));savebin("s0114_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0114_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0114_public_none_success.bin",_rsm_success_read());savebin("s0114_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
32 4194303.9990234375 12884901892 1535.9999997615814 335544320 8388608.0078125 1048576.0009765625 "=" 7340032
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
2.98023224015731e-08 0.0039062499990905053 12.00000000372529 1.430511474609375e-06 0.31249999994179234 0.0078125 0.0009765625 "=" 0.0068359375
end
mata: savebin("s0114_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0114_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0114_public_row_x.bin",st_matrix("r(x)"));savebin("s0114_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0114_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0114_public_row_success.bin",_rsm_success_read());savebin("s0114_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
end
mata: savebin("s0117_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0117_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0117_public_none_x.bin",st_matrix("r(x)"));savebin("s0117_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0117_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0117_public_none_success.bin",_rsm_success_read());savebin("s0117_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
end
mata: savebin("s0117_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0117_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0117_public_row_x.bin",st_matrix("r(x)"));savebin("s0117_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0117_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0117_public_row_success.bin",_rsm_success_read());savebin("s0117_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.25 7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 "=" 0.0068359375
0.25 7.6293945348027137e-06 0.0039062499990905053 0.046875000014551915 0.0003662109375 0.31249999994179234 3.0517578125e-05 "=" 0.0068359375
end
mata: savebin("s0120_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0120_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0120_public_none_x.bin",st_matrix("r(x)"));savebin("s0120_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0120_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0120_public_none_success.bin",_rsm_success_read());savebin("s0120_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
0.25 7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 "=" 0.0068359375
0.25 7.6293945348027137e-06 0.0039062499990905053 0.046875000014551915 0.0003662109375 0.31249999994179234 3.0517578125e-05 "=" 0.0068359375
end
mata: savebin("s0120_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0120_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0120_public_row_x.bin",st_matrix("r(x)"));savebin("s0120_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0120_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0120_public_row_success.bin",_rsm_success_read());savebin("s0120_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
1048576.0009765625 8192 1073741823.75 196608.00006103516 393215.99993896484 85899345920 128.00000011920929 "=" 7340032
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
0.0009765625 7.6293945348027137e-06 0.99999999976716936 0.00018310546880684342 0.0003662109375 79.999999985098839 1.1920928955078125e-07 "=" 0.0068359375
end
mata: savebin("s0123_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0123_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0123_public_none_x.bin",st_matrix("r(x)"));savebin("s0123_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0123_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0123_public_none_success.bin",_rsm_success_read());savebin("s0123_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
1048576.0009765625 8192 1073741823.75 196608.00006103516 393215.99993896484 85899345920 128.00000011920929 "=" 7340032
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
0.0009765625 7.6293945348027137e-06 0.99999999976716936 0.00018310546880684342 0.0003662109375 79.999999985098839 1.1920928955078125e-07 "=" 0.0068359375
end
mata: savebin("s0123_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0123_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0123_public_row_x.bin",st_matrix("r(x)"));savebin("s0123_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0123_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0123_public_row_success.bin",_rsm_success_read());savebin("s0123_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
end
mata: savebin("s0126_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0126_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0126_public_none_x.bin",st_matrix("r(x)"));savebin("s0126_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0126_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0126_public_none_success.bin",_rsm_success_read());savebin("s0126_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
end
mata: savebin("s0126_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0126_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0126_public_row_x.bin",st_matrix("r(x)"));savebin("s0126_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0126_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0126_public_row_success.bin",_rsm_success_read());savebin("s0126_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.25 3.0517578125e-05 0.3125 0.0003662109375 0.046875 0.00390625 7.62939453125e-06 "=" 0.0068359375
0.25 3.0517578125e-05 0.31249999994179234 0.0003662109375 0.046875000014551915 0.0039062499990905053 7.6293945348027137e-06 "=" 0.0068359375
end
mata: savebin("s0129_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0129_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0129_public_none_x.bin",st_matrix("r(x)"));savebin("s0129_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0129_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0129_public_none_success.bin",_rsm_success_read());savebin("s0129_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
0.25 3.0517578125e-05 0.3125 0.0003662109375 0.046875 0.00390625 7.62939453125e-06 "=" 0.0068359375
0.25 3.0517578125e-05 0.31249999994179234 0.0003662109375 0.046875000014551915 0.0039062499990905053 7.6293945348027137e-06 "=" 0.0068359375
end
mata: savebin("s0129_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0129_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0129_public_row_x.bin",st_matrix("r(x)"));savebin("s0129_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0129_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0129_public_row_success.bin",_rsm_success_read());savebin("s0129_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
1048576.0009765625 32768.000030517578 85899345920 1535.9999997615814 50331648.015625 1073741823.75 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.999999985098839 1.430511474609375e-06 0.046875000014551915 0.99999999976716936 2.98023224015731e-08 "=" 0.0068359375
end
mata: savebin("s0132_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0132_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0132_public_none_x.bin",st_matrix("r(x)"));savebin("s0132_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0132_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0132_public_none_success.bin",_rsm_success_read());savebin("s0132_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
1048576.0009765625 32768.000030517578 85899345920 1535.9999997615814 50331648.015625 1073741823.75 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.999999985098839 1.430511474609375e-06 0.046875000014551915 0.99999999976716936 2.98023224015731e-08 "=" 0.0068359375
end
mata: savebin("s0132_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0132_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0132_public_row_x.bin",st_matrix("r(x)"));savebin("s0132_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0132_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0132_public_row_success.bin",_rsm_success_read());savebin("s0132_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
end
mata: savebin("s0135_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0135_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0135_public_none_x.bin",st_matrix("r(x)"));savebin("s0135_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0135_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0135_public_none_success.bin",_rsm_success_read());savebin("s0135_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 268435456.25 "=" 7340032
end
mata: savebin("s0135_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0135_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0135_public_row_x.bin",st_matrix("r(x)"));savebin("s0135_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0135_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0135_public_row_success.bin",_rsm_success_read());savebin("s0135_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
end
mata: savebin("s0138_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0138_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0138_public_none_x.bin",st_matrix("r(x)"));savebin("s0138_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0138_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0138_public_none_success.bin",_rsm_success_read());savebin("s0138_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 256 "=" 7
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 262144.00024414062 "=" 7168
end
mata: savebin("s0138_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0138_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0138_public_row_x.bin",st_matrix("r(x)"));savebin("s0138_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0138_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0138_public_row_success.bin",_rsm_success_read());savebin("s0138_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
3.0517578139210855e-05 3.9999999990686774 12288.000003814697 0.00146484375 319.99999994039536 8 1 "=" 7
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
end
mata: savebin("s0141_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0141_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0141_public_none_x.bin",st_matrix("r(x)"));savebin("s0141_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0141_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0141_public_none_success.bin",_rsm_success_read());savebin("s0141_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
3.0517578139210855e-05 3.9999999990686774 12288.000003814697 0.00146484375 319.99999994039536 8 1 "=" 7
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
0.03125 4095.9999990463257 12582912.00390625 1.4999999997671694 327680 8192.0000076293945 1024.0000009536743 "=" 7168
end
mata: savebin("s0141_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0141_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0141_public_row_x.bin",st_matrix("r(x)"));savebin("s0141_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0141_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0141_public_row_success.bin",_rsm_success_read());savebin("s0141_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
end
mata: savebin("s0144_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0144_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0144_public_none_x.bin",st_matrix("r(x)"));savebin("s0144_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0144_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0144_public_none_success.bin",_rsm_success_read());savebin("s0144_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
0.000244140625 7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 "=" 6.67572021484375e-06
268435456.25 8192 4194303.9990234375 50331648.015625 393215.99993896484 335544320 32768.000030517578 "=" 7340032
end
mata: savebin("s0144_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0144_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0144_public_row_x.bin",st_matrix("r(x)"));savebin("s0144_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0144_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0144_public_row_success.bin",_rsm_success_read());savebin("s0144_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
0.25 7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 "=" 0.0068359375
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
end
mata: savebin("s0147_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0147_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0147_public_none_x.bin",st_matrix("r(x)"));savebin("s0147_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0147_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0147_public_none_success.bin",_rsm_success_read());savebin("s0147_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.0078125 3 64 0.375 320 0.046875 "" .
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
256 0.0078125000036379788 3.9999999990686774 48.000000014901161 0.375 319.99999994039536 0.03125 "=" 7
0.25 7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 "=" 0.0068359375
262144.00024414062 8 4095.9999990463257 49152.000015258789 383.99999994039536 327680 32.000000029802322 "=" 7168
end
mata: savebin("s0147_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0147_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0147_public_row_x.bin",st_matrix("r(x)"));savebin("s0147_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0147_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0147_public_row_success.bin",_rsm_success_read());savebin("s0147_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
1 0.0078125000036379788 1023.9999997615814 0.18750000005820766 0.375 81919.999984741211 0.0001220703125 "=" 7
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
end
mata: savebin("s0150_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0150_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0150_public_none_x.bin",st_matrix("r(x)"));savebin("s0150_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0150_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0150_public_none_success.bin",_rsm_success_read());savebin("s0150_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
1 0.0078125000036379788 1023.9999997615814 0.18750000005820766 0.375 81919.999984741211 0.0001220703125 "=" 7
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
1024.0000009536743 8 1048575.9997558594 192.00000005960464 383.99999994039536 83886080 0.12500000011641532 "=" 7168
end
mata: savebin("s0150_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0150_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0150_public_row_x.bin",st_matrix("r(x)"));savebin("s0150_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0150_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0150_public_row_success.bin",_rsm_success_read());savebin("s0150_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
end
mata: savebin("s0153_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0153_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0153_public_none_x.bin",st_matrix("r(x)"));savebin("s0153_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0153_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0153_public_none_success.bin",_rsm_success_read());savebin("s0153_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
0.000244140625 2.9802322387695312e-08 0.00030517578125 3.5762786865234375e-07 4.57763671875e-05 3.814697265625e-06 7.4505805969238281e-09 "=" 6.67572021484375e-06
268435456.25 32768.000030517578 335544320 393215.99993896484 50331648.015625 4194303.9990234375 8192 "=" 7340032
end
mata: savebin("s0153_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0153_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0153_public_row_x.bin",st_matrix("r(x)"));savebin("s0153_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0153_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0153_public_row_success.bin",_rsm_success_read());savebin("s0153_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
0.25 3.0517578125e-05 0.3125 0.0003662109375 0.046875 0.00390625 7.62939453125e-06 "=" 0.0068359375
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
end
mata: savebin("s0156_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0156_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0156_public_none_x.bin",st_matrix("r(x)"));savebin("s0156_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0156_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0156_public_none_success.bin",_rsm_success_read());savebin("s0156_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
256 0.046875 320 0.375 64 3 0.0078125 "" .
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
256 0.03125 319.99999994039536 0.375 48.000000014901161 3.9999999990686774 0.0078125000036379788 "=" 7
0.25 3.0517578125e-05 0.3125 0.0003662109375 0.046875 0.00390625 7.62939453125e-06 "=" 0.0068359375
262144.00024414062 32.000000029802322 327680 383.99999994039536 49152.000015258789 4095.9999990463257 8 "=" 7168
end
mata: savebin("s0156_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0156_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0156_public_row_x.bin",st_matrix("r(x)"));savebin("s0156_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0156_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0156_public_row_success.bin",_rsm_success_read());savebin("s0156_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
1 0.03125 81919.999984741211 0.00146484375 48.000000014901161 1023.9999997615814 3.0517578139210855e-05 "=" 7
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
end
mata: savebin("s0159_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0159_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0159_public_none_x.bin",st_matrix("r(x)"));savebin("s0159_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0159_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0159_public_none_success.bin",_rsm_success_read());savebin("s0159_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
1 0.03125 81919.999984741211 0.00146484375 48.000000014901161 1023.9999997615814 3.0517578139210855e-05 "=" 7
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
1024.0000009536743 32.000000029802322 83886080 1.4999999997671694 49152.000015258789 1048575.9997558594 0.03125 "=" 7168
end
mata: savebin("s0159_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0159_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0159_public_row_x.bin",st_matrix("r(x)"));savebin("s0159_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0159_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0159_public_row_success.bin",_rsm_success_read());savebin("s0159_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000009094947 3.9999999997671694 48.00000000372529 0.375 319.99999998509884 0.03125 256 "=" 7
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
end
mata: savebin("s0324_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0324_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0324_public_none_x.bin",st_matrix("r(x)"));savebin("s0324_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0324_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0324_public_none_success.bin",_rsm_success_read());savebin("s0324_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000009094947 3.9999999997671694 48.00000000372529 0.375 319.99999998509884 0.03125 256 "=" 7
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
end
mata: savebin("s0324_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0324_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0324_public_row_x.bin",st_matrix("r(x)"));savebin("s0324_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0324_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0324_public_row_success.bin",_rsm_success_read());savebin("s0324_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8.0000000009313226 4095.9999997615814 49152.000003814697 384 327679.99998474121 32 262144 "=" 7168
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
end
mata: savebin("s0327_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0327_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0327_public_none_x.bin",st_matrix("r(x)"));savebin("s0327_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0327_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0327_public_none_success.bin",_rsm_success_read());savebin("s0327_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
8.0000000009313226 4095.9999997615814 49152.000003814697 384 327679.99998474121 32 262144 "=" 7168
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
end
mata: savebin("s0327_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0327_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0327_public_row_x.bin",st_matrix("r(x)"));savebin("s0327_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0327_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0327_public_row_success.bin",_rsm_success_read());savebin("s0327_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
2.9103830456733704e-11 3.814697265625e-06 0.01171875 1.3969838619232178e-09 0.00030517578125 7.62939453125e-06 9.5367431640625e-07 "=" 6.67572021484375e-06
0.031250000003637979 4095.9999997615814 12582912.000976562 1.5 327679.99998474121 8192 1024 "=" 7168
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
end
mata: savebin("s0330_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0330_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0330_public_none_x.bin",st_matrix("r(x)"));savebin("s0330_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0330_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0330_public_none_success.bin",_rsm_success_read());savebin("s0330_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
2.9103830456733704e-11 3.814697265625e-06 0.01171875 1.3969838619232178e-09 0.00030517578125 7.62939453125e-06 9.5367431640625e-07 "=" 6.67572021484375e-06
0.031250000003637979 4095.9999997615814 12582912.000976562 1.5 327679.99998474121 8192 1024 "=" 7168
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
end
mata: savebin("s0330_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0330_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0330_public_row_x.bin",st_matrix("r(x)"));savebin("s0330_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0330_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0330_public_row_success.bin",_rsm_success_read());savebin("s0330_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
9.5367431640625e-07 7.4505805969238281e-09 0.0009765625 1.7881393432617188e-07 3.5762786865234375e-07 0.078125 1.1641532182693481e-10 "=" 6.67572021484375e-06
1024 8.0000000009313226 1048575.9999389648 192.00000001490116 384 83886079.99609375 0.125 "=" 7168
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
end
mata: savebin("s0339_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0339_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0339_public_none_x.bin",st_matrix("r(x)"));savebin("s0339_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0339_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0339_public_none_success.bin",_rsm_success_read());savebin("s0339_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
9.5367431640625e-07 7.4505805969238281e-09 0.0009765625 1.7881393432617188e-07 3.5762786865234375e-07 0.078125 1.1641532182693481e-10 "=" 6.67572021484375e-06
1024 8.0000000009313226 1048575.9999389648 192.00000001490116 384 83886079.99609375 0.125 "=" 7168
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
end
mata: savebin("s0339_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0339_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0339_public_row_x.bin",st_matrix("r(x)"));savebin("s0339_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0339_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0339_public_row_success.bin",_rsm_success_read());savebin("s0339_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.99609375 1.5 49152.000003814697 1048575.9999389648 0.031250000003637979 "=" 7168
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0348_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0348_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0348_public_none_x.bin",st_matrix("r(x)"));savebin("s0348_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0348_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0348_public_none_success.bin",_rsm_success_read());savebin("s0348_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.99609375 1.5 49152.000003814697 1048575.9999389648 0.031250000003637979 "=" 7168
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0348_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0348_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0348_public_row_x.bin",st_matrix("r(x)"));savebin("s0348_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0348_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0348_public_row_success.bin",_rsm_success_read());savebin("s0348_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.99609375 1.5 49152.000003814697 1048575.9999389648 0.031250000003637979 "=" 7168
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0349_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-10)
mata: savebin("s0349_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0349_public_none_x.bin",st_matrix("r(x)"));savebin("s0349_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0349_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0349_public_none_success.bin",_rsm_success_read());savebin("s0349_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
9.5367431640625e-07 2.9802322387695312e-08 0.078125 1.3969838619232178e-09 4.57763671875e-05 0.0009765625 2.9103830456733704e-11 "=" 6.67572021484375e-06
1024 32 83886079.99609375 1.5 49152.000003814697 1048575.9999389648 0.031250000003637979 "=" 7168
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0349_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-10)
mata: savebin("s0349_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0349_public_row_x.bin",st_matrix("r(x)"));savebin("s0349_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0349_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0349_public_row_success.bin",_rsm_success_read());savebin("s0349_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000009094947 3.9999999997671694 48.00000000372529 0.375 319.99999998509884 0.03125 256 "=" 7
end
mata: savebin("s0351_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0351_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0351_public_none_x.bin",st_matrix("r(x)"));savebin("s0351_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0351_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0351_public_none_success.bin",_rsm_success_read());savebin("s0351_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.4505805969238281e-09 3.814697265625e-06 4.57763671875e-05 3.5762786865234375e-07 0.00030517578125 2.9802322387695312e-08 0.000244140625 "=" 6.67572021484375e-06
0.0078125000009094947 3.9999999997671694 48.00000000372529 0.375 319.99999998509884 0.03125 256 "=" 7
end
mata: savebin("s0351_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0351_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0351_public_row_x.bin",st_matrix("r(x)"));savebin("s0351_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0351_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0351_public_row_success.bin",_rsm_success_read());savebin("s0351_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
7.6293945321381784e-06 0.0039062499997726263 0.046875000003637979 0.0003662109375 0.31249999998544808 3.0517578125e-05 0.25 "=" 0.0068359375
end
mata: savebin("s0354_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0354_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0354_public_none_x.bin",st_matrix("r(x)"));savebin("s0354_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0354_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0354_public_none_success.bin",_rsm_success_read());savebin("s0354_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
0.0078125 3 64 0.375 320 0.046875 256 "" .
8 4095.9999997615814 49152.000003814697 383.99999998509884 327680 32.000000007450581 262144.00006103516 "=" 7168
8192 4194303.9997558594 50331648.00390625 393215.99998474121 335544320 32768.000007629395 268435456.0625 "=" 7340032
7.62939453125e-06 0.00390625 0.046875 0.0003662109375 0.3125 3.0517578125e-05 0.25 "=" 0.0068359375
7.6293945321381784e-06 0.0039062499997726263 0.046875000003637979 0.0003662109375 0.31249999998544808 3.0517578125e-05 0.25 "=" 0.0068359375
end
mata: savebin("s0354_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0354_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0354_public_row_x.bin",st_matrix("r(x)"));savebin("s0354_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0354_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0354_public_row_success.bin",_rsm_success_read());savebin("s0354_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
32 4194303.9997558594 12884901889 1535.9999999403954 335544320 8388608.001953125 1048576.0002441406 "=" 7340032
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
2.9802322391164759e-08 0.0039062499997726263 12.000000000931323 1.430511474609375e-06 0.31249999998544808 0.0078125 0.0009765625 "=" 0.0068359375
end
mata: savebin("s0357_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0357_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0357_public_none_x.bin",st_matrix("r(x)"));savebin("s0357_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0357_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0357_public_none_success.bin",_rsm_success_read());savebin("s0357_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
32 4194303.9997558594 12884901889 1535.9999999403954 335544320 8388608.001953125 1048576.0002441406 "=" 7340032
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
2.9802322391164759e-08 0.0039062499997726263 12.000000000931323 1.430511474609375e-06 0.31249999998544808 0.0078125 0.0009765625 "=" 0.0068359375
end
mata: savebin("s0357_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0357_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0357_public_row_x.bin",st_matrix("r(x)"));savebin("s0357_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0357_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0357_public_row_success.bin",_rsm_success_read());savebin("s0357_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
1048576.0002441406 8192 1073741823.9375 196608.00001525879 393215.99998474121 85899345920 128.00000002980232 "=" 7340032
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
0.0009765625 7.6293945321381784e-06 0.99999999994179234 0.00018310546876421085 0.0003662109375 79.99999999627471 1.1920928955078125e-07 "=" 0.0068359375
end
mata: savebin("s0366_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0366_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0366_public_none_x.bin",st_matrix("r(x)"));savebin("s0366_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0366_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0366_public_none_success.bin",_rsm_success_read());savebin("s0366_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
1048576.0002441406 8192 1073741823.9375 196608.00001525879 393215.99998474121 85899345920 128.00000002980232 "=" 7340032
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
0.0009765625 7.6293945321381784e-06 0.99999999994179234 0.00018310546876421085 0.0003662109375 79.99999999627471 1.1920928955078125e-07 "=" 0.0068359375
end
mata: savebin("s0366_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0366_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0366_public_row_x.bin",st_matrix("r(x)"));savebin("s0366_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0366_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0366_public_row_success.bin",_rsm_success_read());savebin("s0366_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1048576.0002441406 32768.000007629395 85899345920 1535.9999999403954 50331648.00390625 1073741823.9375 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.99999999627471 1.430511474609375e-06 0.046875000003637979 0.99999999994179234 2.9802322391164759e-08 "=" 0.0068359375
end
mata: savebin("s0375_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0375_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0375_public_none_x.bin",st_matrix("r(x)"));savebin("s0375_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0375_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0375_public_none_success.bin",_rsm_success_read());savebin("s0375_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1048576.0002441406 32768.000007629395 85899345920 1535.9999999403954 50331648.00390625 1073741823.9375 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.99999999627471 1.430511474609375e-06 0.046875000003637979 0.99999999994179234 2.9802322391164759e-08 "=" 0.0068359375
end
mata: savebin("s0375_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0375_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0375_public_row_x.bin",st_matrix("r(x)"));savebin("s0375_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0375_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0375_public_row_success.bin",_rsm_success_read());savebin("s0375_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1048576.0002441406 32768.000007629395 85899345920 1535.9999999403954 50331648.00390625 1073741823.9375 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.99999999627471 1.430511474609375e-06 0.046875000003637979 0.99999999994179234 2.9802322391164759e-08 "=" 0.0068359375
end
mata: savebin("s0376_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-10)
mata: savebin("s0376_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0376_public_none_x.bin",st_matrix("r(x)"));savebin("s0376_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0376_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0376_public_none_success.bin",_rsm_success_read());savebin("s0376_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1048576.0002441406 32768.000007629395 85899345920 1535.9999999403954 50331648.00390625 1073741823.9375 32 "=" 7340032
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
0.0009765625 3.0517578125e-05 79.99999999627471 1.430511474609375e-06 0.046875000003637979 0.99999999994179234 2.9802322391164759e-08 "=" 0.0068359375
end
mata: savebin("s0376_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-10)
mata: savebin("s0376_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0376_public_row_x.bin",st_matrix("r(x)"));savebin("s0376_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0376_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0376_public_row_success.bin",_rsm_success_read());savebin("s0376_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
3.0517578128552714e-05 3.9999999997671694 12288.000000953674 0.00146484375 319.99999998509884 8 1 "=" 7
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
end
mata: savebin("s0384_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0384_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0384_public_none_x.bin",st_matrix("r(x)"));savebin("s0384_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0384_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0384_public_none_success.bin",_rsm_success_read());savebin("s0384_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
3.0517578125e-05 3 16384 0.00146484375 320 12 1 "" .
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
3.0517578128552714e-05 3.9999999997671694 12288.000000953674 0.00146484375 319.99999998509884 8 1 "=" 7
2.9802322387695312e-08 0.00390625 12 1.430511474609375e-06 0.3125 0.0078125 0.0009765625 "=" 0.0068359375
0.03125 4095.9999997615814 12582912.000976562 1.4999999999417923 327680 8192.0000019073486 1024.0000002384186 "=" 7168
end
mata: savebin("s0384_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0384_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0384_public_row_x.bin",st_matrix("r(x)"));savebin("s0384_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0384_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0384_public_row_success.bin",_rsm_success_read());savebin("s0384_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
1 0.0078125000009094947 1023.9999999403954 0.18750000001455192 0.375 81919.999996185303 0.0001220703125 "=" 7
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
end
mata: savebin("s0393_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0393_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0393_public_none_x.bin",st_matrix("r(x)"));savebin("s0393_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0393_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0393_public_none_success.bin",_rsm_success_read());savebin("s0393_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.0078125 768 0.25 0.375 81920 0.00018310546875 "" .
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
1 0.0078125000009094947 1023.9999999403954 0.18750000001455192 0.375 81919.999996185303 0.0001220703125 "=" 7
0.0009765625 7.62939453125e-06 1 0.00018310546875 0.0003662109375 80 1.1920928955078125e-07 "=" 0.0068359375
1024.0000002384186 8 1048575.9999389648 192.00000001490116 383.99999998509884 83886080 0.12500000002910383 "=" 7168
end
mata: savebin("s0393_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0393_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0393_public_row_x.bin",st_matrix("r(x)"));savebin("s0393_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0393_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0393_public_row_success.bin",_rsm_success_read());savebin("s0393_public_row_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1 0.03125 81919.999996185303 0.00146484375 48.00000000372529 1023.9999999403954 3.0517578128552714e-05 "=" 7
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0402_public_none_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(none) tolerance(1e-08)
mata: savebin("s0402_public_none_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0402_public_none_x.bin",st_matrix("r(x)"));savebin("s0402_public_none_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0402_public_none_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0402_public_none_success.bin",_rsm_success_read());savebin("s0402_public_none_counts.bin",_rsm_success_counts_read())
mata: _rsm_success_reset(@POLICY@,1)
clear
input double x1 double x2 double x3 double x4 double x5 double x6 double x7 str2 rel double rhs
1 0.046875 81920 0.00146484375 64 768 3.0517578125e-05 "" .
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
1 0.03125 81919.999996185303 0.00146484375 48.00000000372529 1023.9999999403954 3.0517578128552714e-05 "=" 7
0.0009765625 3.0517578125e-05 80 1.430511474609375e-06 0.046875 1 2.9802322387695312e-08 "=" 0.0068359375
1024.0000002384186 32.000000007450581 83886080 1.4999999999417923 49152.000003814697 1048575.9999389648 0.03125 "=" 7168
end
mata: savebin("s0402_public_row_input.bin",st_data(.,tokens("x1 x2 x3 x4 x5 x6 x7 rhs")))
lp x1 x2 x3 x4 x5 x6 x7,min engine(simplex) scaling(row) tolerance(1e-08)
mata: savebin("s0402_public_row_result.bin",(st_numscalar("r(status)"),st_numscalar("r(objective)"),st_numscalar("r(iterations)")))
mata: savebin("s0402_public_row_x.bin",st_matrix("r(x)"));savebin("s0402_public_row_attempts.bin",st_matrix("r(attempts)"))
mata: savebin("s0402_public_row_public_residual.bin",(st_numscalar("r(max_eq_violation)"),st_numscalar("r(max_scaled_eq_violation)"),st_numscalar("r(objective_scale)"),st_numscalar("r(row_scale_min)"),st_numscalar("r(row_scale_max)")))
mata: savebin("s0402_public_row_success.bin",_rsm_success_read());savebin("s0402_public_row_counts.bin",_rsm_success_counts_read())

mata: _rsm_success_reset(@POLICY@,1)
mata:
E=(1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125\1,0.03125,81919.999996185303,0.00146484375,48.00000000372529,1023.9999999403954,3.0517578128552714e-05\0.0009765625,3.0517578125e-05,80,1.430511474609375e-06,0.046875,1,2.9802322387695312e-08\1024.0000002384186,32.000000007450581,83886080,1.4999999999417923,49152.000003814697,1048575.9999389648,0.03125)
b=(7168\7\0.0068359375\7168)
c=(1,0.046875,81920,0.00146484375,64,768,3.0517578125e-05)
P=RSM_Prob()
_rsm_std_form(P,J(0,cols(E),0),J(0,1,0),E,b,c)
_rsm_audit_reset(1)
R=_rsm_two_phase(P,1e-10,"")
savebin("boundary_result.bin",(R.status,R.objval,R.iters));savebin("boundary_x.bin",R.x)
end

mata: savebin("boundary_success.bin",_rsm_success_read());savebin("boundary_counts.bin",_rsm_success_counts_read())
