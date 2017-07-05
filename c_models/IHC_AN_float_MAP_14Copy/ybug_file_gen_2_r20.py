address_out=['0x60640018','0x60670ad4','0x606a1590','0x606d204c','0x60702b08','0x607335c4',
		     '0x60764080','0x60794b3c','0x607c55f8','0x607f60b4','0x60826b70','0x6085762c',
			'0x608880e8','0x608b8ba4','0x608e9660','0x6091a11c']

address_in=['0x60644508','0x60674fc4','0x606a5a80','0x606d653c','0x60706ff8','0x60737ab4',
		     '0x60768570','0x6079902c','0x607c9ae8','0x607fa5a4','0x6082b060','0x6085bb1c',
			'0x6088c5d8','0x608bd094','0x608edb50','0x6091e60c']

profile_out=['0x6066f620','0x606a00dc','0x606d0b98','0x60701654','0x60732110','0x60762bcc',
		'0x60793688','0x607c4144','0x607f4c00','0x608256bc','0x60856178','0x60886c34',
			'0x608b76f0','0x608e81ac','0x60918c68','0x60949724']


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
				#single channel hack
				#line="sload ./load_files/load{}_{} {}\n".format(1,j+1,address_in[i])
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


test_script_gen(segments=240,chips=1,dur=10,boot_string="boot scamp.boot no_wdog.conf",cores=16,segment_length=100//20,num_fibres=2,app_no=20)


