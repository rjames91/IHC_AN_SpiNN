address_out=['0x60640018','0x60718a2c','0x607f1440','0x608c9e54','0x609a2868','0x60a7b27c',
		     '0x60b53c90','0x60c2c6a4','0x60d050b8','0x60dddacc','0x60eb64e0','0x60f8eef4',
			'0x61067908','0x6114031c','0x61218d30','0x612f1744']

address_in=['0x606ec460','0x607c4e74','0x6089d888','0x6097629c','0x60a4ecb0','0x60b276c4',
		     '0x60c000d8','0x60cd8aec','0x60db1500','0x60e89f14','0x60f62928','0x6103b33c',
			'0x61113d50','0x611ec764','0x612c5178','0x6139db8c']

profile_out=['0x60717578','0x607eff8c','0x608c89a0','0x609a13b4','0x60a79dc8','0x60b527dc',
		'0x60c2b1f0','0x60d03c04','0x60ddc618','0x60eb502c','0x60f8da40','0x61066454',
			'0x6113ee68','0x6121787c','0x612f0290','0x613c8ca4']

chip_list=['0 0','1 0','1 1','0 1','2 0','2 1','2 2','1 2','0 2','3 0','3 1','3 2','4 3','3 3','2 3','1 3','4 0','4 1','4 2','5 3','5 4','4 4','3 4','2 4','0 3','5 1','5 2','6 3','6 4','6 5','6 6','5 5','4 6','4 5','3 5','1 4','6 2','7 3','7 4','7 5','7 6','7 7','6 7','5 6','5 7','4 7','3 6','2 5']

#print len(chip_list)

def input_gen(chips=1,segments=1,cores=16):

	f=open('./timed_input_write','w')
	#f.write('test.\n')
	for j in range(segments):
		for k in range(chips):
			#f.write("sp {} {}\n".format(k>>1,k%2))
			f.write("sp {}\n".format(chip_list[k]))
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



	bytes=segment_length*num_seg*num_fibres*4 
	profile_bytes=3*num_seg*4
	f=open('./RAM_dump','w')
	for k in range(chips):
		#f.write("sp {} {}\n".format(k>>1,k%7))#2))
		f.write("sp {}\n".format(chip_list[k]))
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
		#f.write("sp {} {}\n".format(i>>1,i%2))
		f.write("sp {}\n".format(chip_list[i]))
		#load application 
		f.write("@ application_load_single_chip\n")
	f.write("@ timed_input_write\n")
	f.write("sleep {}\n".format(dur))	
	f.write("@ RAM_dump\n")


test_script_gen(segments=200,chips=48,dur=20,boot_string="boot scamp.boot no_wdog.conf",cores=16,segment_length=100//5,num_fibres=20,app_no=20)


