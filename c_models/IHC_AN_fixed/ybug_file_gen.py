def input_gen(chips=1,segments=1,cores=16):

	address_in=['0x60696240','0x60744ff0','0x607f3da0','0x608a2b50','0x60951900','0x60a006b0',
			     '0x60aaf460','0x60b5e210','0x60c0cfc0','0x60cbbd70','0x60d6ab20','0x60e198d0',
				'0x60ec8680','0x60f77430','0x610261e0','0x610d4f90']

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

def dump(chips,segment_length,num_seg,num_fibres=10,cores=16,word_length=64):


	address_out=['0x60640018','0x606eedc8','0x6079db78','0x6084c928','0x608fb6d8','0x609aa488',
			     '0x60a59238','0x60b07fe8','0x60bb6d98','0x60c65b48','0x60d148f8','0x60dc36a8',
				'0x60e72458','0x60f21208','0x60fcffb8','0x6107ed68']


	profile_out=['0x606ec468','0x6079b218','0x60849fc8','0x608f8d78','0x609a7b28','0x60a568d8',
			'0x60b05688','0x60bb4438','0x60c631e8','0x60d11f98','0x60dc0d48','0x60e6faf8',
				'0x60f1e8a8','0x60fcd658','0x6107c408','0x6112b1b8']
	bytes=segment_length*num_seg*num_fibres*(word_length/8) 
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





