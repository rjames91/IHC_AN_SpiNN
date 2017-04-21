def input_gen(chips=1,segments=1,cores=16):


	address_in=['0x60cfaaa0','0x613e1af4','0x61ac8b48','0x621afb9c','0x62896bf0','0x62f7dc44',
			     '0x63664c98','0x63d4bcec','0x64432d40','0x64b19d94','0x65200de8','0x658e7e3c',
				'0x65fcee90','0x666b5ee4','0x66d9cf38','0x67483f8c']


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

	address_out=['0x60640018','0x60d2706c','0x6140e0c0','0x61af5114','0x621dc168','0x628c31bc',
			     '0x62faa210','0x63691264','0x63d782b8','0x6445f30c','0x64b46360','0x6522d3b4',
				'0x65914408','0x65ffb45c','0x666e24b0','0x66dc9504']


	profile_out=['0x60d25bb8','0x6140cc0c','0x61af3c60','0x621dacb4','0x628c1d08','0x62fa8d5c',
			'0x6368fdb0','0x63d76e04','0x6445de58','0x64b44eac','0x6522bf00','0x65912f54',
				'0x65ff9fa8','0x666e0ffc','0x66dc8050','0x674af0a4']

	bytes=segment_length*num_seg*num_fibres*4 
	profile_bytes=3*num_seg*4
	f=open('./RAM_dump','w')
	for k in range(chips):
		f.write("sp {} {}\n".format(k>>1,k%2))
		for i in range(cores):	
			line="sdump ./dump_files/dump{}.bin {} {}\n".format((cores*k)+i+1,address_out[i],hex(bytes))
			f.write(line)
			#profile dump
			#f.write("sdump ./dump_files/profile{}.bin {} {}\n".format((cores*k)+i+1,profile_out[i],hex(profile_bytes)))	
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





