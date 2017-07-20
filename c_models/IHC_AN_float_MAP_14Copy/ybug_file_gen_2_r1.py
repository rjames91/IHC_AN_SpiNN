address_out=['0x60640018','0x606c280c','0x60745000','0x607c77f4','0x60849fe8','0x608cc7dc',
             '0x6094efd0','0x609d17c4','0x60a53fb8','0x60ad67ac','0x60b58fa0','0x60bdb794',
            '0x60c5df88','0x60ce077c','0x60d62f70','0x60de5764']

address_in=['0x60696240','0x60718a34','0x6079b228','0x6081da1c','0x608a0210','0x60922a04',
             '0x609a51f8','0x60a279ec','0x60aaa1e0','0x60b2c9d4','0x60baf1c8','0x60c319bc',
            '0x60cb41b0','0x60d369a4','0x60db9198','0x60e3b98c']

profile_out=['0x606c1358','0x60743b4c','0x607c6340','0x60848b34','0x608cb328','0x6094db1c',
        '0x609d0310','0x60a52b04','0x60ad52f8','0x60b57aec','0x60bda2e0','0x60c5cad4',
            '0x60cdf2c8','0x60d61abc','0x60de42b0','0x60e66aa4']

chip_list=['0 0','1 0','1 1','0 1','2 0','2 1','2 2','1 2','0 2','3 0','3 1','3 2','4 3','3 3','2 3','1 3','4 0','4 1','4 2','5 3','5 4',
           '4 4','3 4','2 4','0 3','5 1','5 2','6 3','6 4','6 5','6 6','5 5','4 6','4 5','3 5','1 4','6 2','7 3','7 4','7 5','7 6','7 7',
           '6 7','5 6','5 7','4 7','3 6','2 5']

#print len(chip_list)

def input_gen(chips=1,segments=1,cores=16):

    f=open('./timed_input_write','w')
    #f.write('test.\n')
    for j in range(segments):
        for k in range(chips):
            f.write("sp {}\n".format(chip_list[k]))
            for i in range(cores):
                #line="sload ./load_files/load{}_{} {}\n".format((cores*k)+i+1,j+1,address_in[i])
                #single channel hack
                line="sload ./load_files/load{}_{} {}\n".format(1,j+1,address_in[i])
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


test_script_gen(segments=61,chips=48,dur=5,boot_string="boot scamp.boot no_wdog.conf",cores=16,segment_length=100,num_fibres=2,app_no=20)


