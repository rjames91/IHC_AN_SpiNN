address_out=['0x60640018','0x60674560','0x606a8aa8','0x606dcff0','0x60711538','0x60745a80',
		     '0x60779fc8','0x607ae510','0x607e2a58','0x60816fa0','0x6084b4e8','0x6087fa30',
			'0x608b3f78','0x608e84c0','0x6091ca08','0x60950f50']

address_in=['0x606489f0','0x6067cf38','0x606b1480','0x606e59c8','0x60719f10','0x6074e458',
		     '0x607829a0','0x607b6ee8','0x607eb430','0x6081f978','0x60853ec0','0x60888408',
			'0x608bc950','0x608f0e98','0x609253e0','0x60959928']

profile_out=['0x60673b08','0x606a8050','0x606dc598','0x60710ae0','0x60745028','0x60779570',
		'0x607adab8','0x607e2000','0x60816548','0x6084aa90','0x6087efd8','0x608b3520',
			'0x608e7a68','0x6091bfb0','0x609504f8','0x60984a40']


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


test_script_gen(segments=100,chips=4,dur=20,boot_string="boot scamp.boot no_wdog.conf",cores=16,segment_length=200//20,num_fibres=4,app_no=20)


