def input_gen(chips=1,segments=1,cores=15):

	address_in=['0x607eeac0','0x609c9b34','0x60ba4ba8','0x60d7fc1c','0x60f5ac90','0x61135d04',
			     '0x61310d78','0x614ebdec','0x616c6e60','0x618a1ed4','0x61a7cf48','0x61c57fbc',
				'0x61e33030','0x6200e0a4','0x621e9118','']

	f=open('./timed_input_write','w')
	#f.write('test.\n')
	for j in range(segments):
		for k in range(chips):
			f.write("sp {} {}\n".format(k>>1,k%2))
			for i in range(cores):
				#line="sload ./load_files/load{}_{} {}\n".format((cores*k)+i+1,j+1,address_in[i])
				line="sload ./load_files/load{}_{} {}\n".format(((cores/5)*k)+(i/5)+1,j+1,address_in[i])
				f.write(line)
		if j==0:	
			f.write("app_sig all 20 sync0\n")
	#		f.write("sleep 0.5\n")
		else:
			f.write("sleep\n")
	f.close()

def dump(chips,segment_length,num_seg,num_fibres=2,cores=15):


	address_out=['0x60640018','0x6081b08c','0x609f6100','0x60bd1174','0x60dac1e8','0x60f8725c',
			     '0x611622d0','0x6133d344','0x615183b8','0x616f342c','0x618ce4a0','0x61aa9514',
				'0x61c84588','0x61e5f5fc','0x6203a670','']


	profile_out=['0x60819bd8','0x609f4c4c','0x60bcfcc0','0x60daad34','0x60f85da8','0x61160e1c',
			'0x6133be90','0x61516f04','0x616f1f78','0x618ccfec','0x61aa8060','0x61c830d4',
				'0x61e5e148','0x620391bc','0x62214230','']
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
		
def test_script_gen(segments,chips=1,dur=7,boot_string="boot",cores=15,segment_length=100,num_fibres=2,app_no=20):

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





