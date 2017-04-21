address_out=['0x60640018','0x606721cc','0x606a4380','0x606d6534','0x607086e8','0x6073a89c',
		     '0x6076ca50','0x6079ec04','0x607d0db8','0x60802f6c','0x60835120','0x608672d4',
			'0x60899488','0x608cb63c','0x608fd7f0','0x6092f9a4']

address_in=['0x60645c00','0x60677db4','0x606a9f68','0x606dc11c','0x6070e2d0','0x60740484',
		     '0x60772638','0x607a47ec','0x607d69a0','0x60808b54','0x6083ad08','0x6086cebc',
			'0x6089f070','0x608d1224','0x609033d8','0x6093558c']

profile_out=['0x60670d18','0x606a2ecc','0x606d5080','0x60707234','0x607393e8','0x6076b59c',
		'0x6079d750','0x607cf904','0x60801ab8','0x60833c6c','0x60865e20','0x60897fd4',
			'0x608ca188','0x608fc33c','0x6092e4f0','0x609606a4']

chip_list=['0 0','1 0','1 1','0 1','2 0','2 1','2 2','1 2','0 2','3 0','3 1','3 2','4 3','3 3','2 3','1 3','4 0','4 1','4 2','5 3','5 4','4 4','3 4','2 4','0 3','5 1','5 2','6 3','6 4','6 5','6 6','5 5','4 6','4 5','3 5','1 4','6 2','7 3','7 4','7 5','7 6','7 7','6 7','5 6','5 7','4 7','3 6','2 5']

#print len(chip_list)

def input_gen(chips=1,segments=1,cores=16):

	f=open('./timed_input_write','w')
	#f.write('test.\n')
	for j in range(segments):
		for k in range(chips):
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
		f.write("sp {}\n".format(chip_list[i]))
		#load application 
		f.write("@ application_load_single_chip\n")
	f.write("@ timed_input_write\n")
	f.write("sleep {}\n".format(dur))	
	f.write("@ RAM_dump\n")


test_script_gen(segments=200,chips=4,dur=20,boot_string="boot scamp.boot no_wdog.conf",cores=16,segment_length=100//30,num_fibres=4,app_no=20)


