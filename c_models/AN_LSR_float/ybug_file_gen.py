def input_gen(chips=1,segments=1,cores=16):

#	address_in=['0x607eeac0','0x60b4c010','0x60ea9560','0x61206ab0','0x61564000','0x618c1550',
#		     '0x61c1eaa0','0x61f7bff0','0x622d9540','0x62636a90','0x62993fe0','0x62cf1530',
#			'0x6304ea80','0x633abfd0','0x63709520','0x63a66a70','0x63dc3fc0']
	address_in=['0x607eeac0','0x60b4d4c4','0x60eabec8','0x6120a8cc','0x615692d0','0x618c7cd4',
			     '0x61c266d8','0x61f850dc','0x622e3ae0','0x626424e4','0x629a0ee8','0x62cff8ec',
				'0x6305e2f0','0x633bccf4','0x6371b6f8','0x63a7a0fc']

	f=open('./timed_input_write','w')
	#f.write('test.\n')
	for j in range(segments):
		for k in range(chips):
			f.write("sp {} {}\n".format(k>>1,k%2))
			for i in range(cores):
				line="sload ./load_files/load{}_{} {}\n".format((cores*k)+i+1,j+1,address_in[i])
				f.write(line)
		if j==0:	
			f.write("app_sig all 20 sync0\n")
	#		f.write("sleep 0.5\n")
		else:
			f.write("sleep\n")
	f.close()

def dump(chips,segment_length,num_seg,num_fibres=10,cores=16):

#	address_out=['0x60640018','0x6099d568','0x60cfaab8','0x61058008','0x613b5558','0x61712aa8',
#			     '0x61a6fff8','0x61dcd548','0x6212aa98','0x62487fe8','0x627e5538','0x62b42a88',
#				'0x62e9ffd8','0x631fd528','0x6355aa78','0x638b7fc8','0x63c15518']

	address_out=['0x60640018','0x6099ea1c','0x60cfd420','0x6105be24','0x613ba828','0x6171922c',
			     '0x61a77c30','0x61dd6634','0x62135038','0x62493a3c','0x627f2440','0x62b50e44',
				'0x62eaf848','0x6320e24c','0x6356cc50','0x638cb654']

	profile_out=['0x6099d568','0x60cfbf6c','0x6105a970','0x613b9374','0x61717d78','0x61a7677c',
			'0x61dd5180','0x62133b84','0x62492588','0x627f0f8c','0x62b4f990','0x62eae394',
				'0x6320cd98','0x6356b79c','0x638ca1a0','0x63c28ba4']
	bytes=segment_length*num_seg*num_fibres*4 
	profile_bytes=3*num_seg*4
	f=open('./RAM_dump','w')
	for k in range(chips):
		f.write("sp {} {}\n".format(k>>1,k%2))
		for i in range(cores):	
			line="sdump ./dump_files/dump{}.bin {} {}\n".format((cores*k)+i+1,address_out[i],hex(bytes))
			f.write(line)
			#profile dump
			f.write("sdump ./dump_files/profile{}.bin {} {}\n".format((cores*k)+i+1,profile_out[i],hex(profile_bytes)))	
	f.close()
		
def test_script_gen(segments,chips=1,dur=7,boot_string="boot",cores=16,segment_length=100,num_fibres=10,app_no=20):

	#generate up to date timed_input_write file
	input_gen(chips,1,cores)#1 segment inputs hard coded
	#generate up to date dump file
	dump(chips,segment_length,segments,num_fibres,cores)

	f=open('./test','w')
#	f.write("boot scamp.boot no_wdog.conf\n")	
	f.write("{}\n".format(boot_string))
	f.write("app_stop {}\n".format(app_no))	
	for i in range(chips):
		#switch chip
		f.write("sp {} {}\n".format(i>>1,i%2))
		#load application 
		f.write("@ application_load_single_chip\n")
	f.write("@ timed_input_write\n")
	f.write("sleep {}\n".format(dur))	
	f.write("@ RAM_dump\n")





