
startupVars generateStartupVars(void)
{
	startupVars out;
	REAL Gu0,kt0,IHCV,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,tauca,CaCurrLSR,CaCurrMSR,CaCurrHSR,kt0LSR,kt0MSR,kt0HSR,
	ANCleftLSR,ANCleftMSR,ANCleftHSR,ANAvailLSR,ANAvailMSR,ANAvailHSR,ANReproLSR,ANReproMSR,
	ANReproHSR,z,power,y,x,l,r,MLSR,MMSR,MHSR;

	//======Calculate Gu0===========//
	gmax=6e-9;//(float)Cilia.Gmax/pow(2,30);
	u0=5e-9;//(float)Cilia.u0/pow(2,30);
	u1=1e-9;//(float)Cilia.u1/pow(2,30);
	s1=1e-9;//(float)Cilia.recips1/pow(2,1);
	s0=3e-8;//(float)Cilia.recips0/pow(2,1);
	ga=8e-10;//(float)Cilia.Ga/pow(2,31);
	Gu0=ga+gmax/(1+(exp(u0/s0)*(1+exp(u1/s1))));
	float test=exp(u0/s0);
	float test1=exp(u1/s1);
	//======Calculate IHCVnow=========//
	gk=2e-8;//(float)Cilia.Gk/pow(2,30);
	et=0.1;//(float)Cilia.Et/pow(2,30);
	ek=-0.08;//(float)Cilia.Ek/pow(2,30);
	rpc=0.04;//(float)Cilia.Rpc/pow(2,30);
	ekp=ek+rpc*et;
	IHCV=(gk*ekp+Gu0*et)/(Gu0+gk);

	//=======Calculate mICaCurrent========//
	gamma=130;
	beta=400;
	mICaCurr=1/(1+exp(-gamma*IHCV)*(1/beta));

	//=======Calculate CaCurrent========//
	gmaxcalsr=1e-11;
	gmaxcamsr=7e-11;
	gmaxcahsr=5e-10;
	eca=0.066;
	tauca=1e-4;
	CaCurrLSR=gmaxcalsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;
	CaCurrMSR=gmaxcamsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;
	CaCurrHSR=gmaxcahsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;

	//========Calculate AN==========//
	z=1e38;//reduced from 1e40 fit as single float type
	power=3;
	y=10;
	x=30;
	MLSR=13;
	MMSR=14;
	MHSR=8;
	l=2580;
	r=6580;

	kt0LSR=z*pow(CaCurrLSR,power)*100;//added *100
	kt0MSR=z*pow(CaCurrMSR,power)*100;
	kt0HSR=z*pow(CaCurrHSR,power)*100;
	ANCleftLSR=(kt0LSR*y*MLSR)/(y*(l+r)+kt0LSR*l);
	ANCleftMSR=(kt0MSR*y*MMSR)/(y*(l+r)+kt0MSR*l);
	ANCleftHSR=(kt0HSR*y*MHSR)/(y*(l+r)+kt0HSR*l);
	ANAvailLSR=round((ANCleftLSR*(l+r))/kt0LSR);
	ANAvailMSR=round((ANCleftMSR*(l+r))/kt0MSR);
	ANAvailHSR=round((ANCleftHSR*(l+r))/kt0HSR);
	ANReproLSR=(ANCleftLSR*r)/x;
	ANReproMSR=(ANCleftMSR*r)/x;
	ANReproHSR=(ANCleftHSR*r)/x;

	out.IHCVnow0=IHCV;
	out.Ekp=ekp;
	out.mICaCurr0=mICaCurr;
	out.CaCurrLSR0=CaCurrLSR;
	out.CaCurrMSR0=CaCurrMSR;
	out.CaCurrHSR0=CaCurrHSR;
	out.ANCleftLSR0=ANCleftLSR;
	out.ANCleftMSR0=ANCleftMSR;
	out.ANCleftHSR0=ANCleftHSR;
	out.ANAvailLSR0=ANAvailLSR;
	out.ANAvailMSR0=ANAvailMSR;
	out.ANAvailHSR0=ANAvailHSR;
	out.ANReproLSR0=ANReproLSR;
	out.ANReproMSR0=ANReproMSR;
	out.ANReproHSR0=ANReproHSR;

	return out;
}