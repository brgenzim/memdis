from lxml import etree
from os import popen
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
from tensorflow.keras.models import load_model



def readmod(installDir):

	MODELI=load_model(installDir+"/models/ModelInside0610")
	MODELi=load_model(installDir+"/models/Modelinside0610")
	MODELo=load_model(installDir+"/models/Modeloutside0610")
	MODELO=load_model(installDir+"/models/ModelOutside0610")
	MODELlstm=load_model(installDir+"/models/ModelBiDirLSTM0610")

	return [MODELI,MODELi,MODELo,MODELO,MODELlstm]

def readfasta(x):
	infile=open(x)
	s=""
	header=""
	while 1:
		l=infile.readline()
		if l=="":
			break
		if l[0]==">":
			header=l[1:].strip()
		else:
			s=s+l.strip()
	infile.close()
	return [header,s]

def readcctop(x):
	t=""
	doc=etree.parse(x)
	for elt in doc.getiterator():
		if elt.tag=="Topology":
			t=""
			for elt2 in elt.getiterator():
				if elt2.tag=="Region":
					for i in range (int(elt2.attrib["from"]),int(elt2.attrib["to"])+1):
						t=t+elt2.attrib["loc"]
	return t

def readtopology(f):
	t=""
	if f[-3:]=="xml":
		t=readcctop(f)
	else:
		[h,t]=readfasta(f)
	return t

def readnetsurfp(x):
	infile=open(x)
	res=[]
	while 1:
		l=infile.readline()
		if l=="":
			break
		if l[0]!="#":
			res.append([0,0,0,0])
			l=l.split()
			res[len(res)-1][0]=l[4]
			res[len(res)-1][1]=l[7]
			res[len(res)-1][2]=l[8]
			res[len(res)-1][3]=l[9]
	infile.close()
	return res
def readpsiblast(x):
	ABC="RKHDESTYFWALVIMGNQPC"
	infile=open(x)
	l=infile.readline()
	l=infile.readline()
	l=infile.readline()
	ranges=l.split()
	consi=[]
	hit=[]
	mx={}
	mx2={}
	n=0
	while 1:

		l=infile.readline()
		l=l.split()
		if len(l)<20:
			break
		mx[n]={}
		total=0
		cons=0
		for i in range (22,42):
			total=total+int(l[i])
			mx[n][ranges[i-2]]=int(l[i])
		hit.append(total)

		for key in mx[n]:
			if cons>0:
				if mx[n][key]/total>cons:
					cons=mx[n][key]/total
			if total>0:
				mx[n][key]=mx[n][key]/total
			else:
				mx[n][key]=0
		consi.append(cons)
		n+=1

	for i in range (0, len(mx)):
		mx2[i]=[]
		for j in range (0, len(ABC)):
			mx2[i].append(mx[i][ABC[j]])

	infile.close()
	return [hit,consi,mx2]

def readseg(x):
        data2=""
        data=""
        a=open(x)
        while 1:
                al=a.readline()
                if al=="":
                        break
                if al[:12]=="{'proteins':":
                        data=al.strip().replace("'",'"')
                elif al.find("sequence")!=-1:
                        data2=al.strip().replace("'",'"')
        s=[]

        proc_dat = eval(data)
        proc_dat2 = eval(data2)
        for key in proc_dat2:
                if key=="sequence":
                        for i in range (0, len(proc_dat2[key])):
                                s.append("0")
        for key in proc_dat:
                if key=="proteins":
                        for key2 in proc_dat[key][0]:
                                if key2=="SEG":

                                        for j in range (0, len(proc_dat[key][0][key2])):

                                                for k in range (int(proc_dat[key][0][key2][j][0])-1,int(proc_dat[key][0][key2][j][1])):
                                                        s[k]="1"
        return s

def calctop(top):
	tail=[]
	numtm=[]
	segment=[]
	dist=[]
	ntm=0
	for i in range (0, len(top)-1):
		if top[i]=="M" and top[i+1]!="M":
			ntm+=1
	for i in range (0, len(top)):
		pos=i
		negative=99999
		positive=99999
		boolmem=[False,False]
		while 1:
			if pos==len(top):
				break
			if top[pos]=="M":
				boolmem[0]=True
				break

			else:
				pos+=1
		positive=pos-i
		pos=i
		while 1:
			if pos==-1:
				break
			if top[pos]=="M":
				boolmem[1]=True
				break

			else:
				pos-=1
		negative=i-pos
		if boolmem[0]==True and boolmem[1]==True:
			tail.append(0)
		else:
			tail.append(1)

		if positive+negative<100:
			segment.append((positive+negative)/100)
		else:
			segment.append(1)
		if ntm>1:
			numtm.append(0)
		else:
			numtm.append(1)
		dist.append(min(positive,negative))
	memleng=[]
	for j in range (0, len(dist)):
		memleng.append(0)
	addl=0
	saver=0
	bool1=False
	for j in range (0, len(dist)):

		if dist[j]==0:
			addl+=1
			bool1=True

		else:
			if bool1==True:
				for k in range (j,j+10):
					try:
						if dist[k]!=0:
							memleng[k]=(addl-10)/20
					except IndexError:
						pass
				bool1=False
				addl=0
	addl=0
	bool1=False
	for j in range (len(dist)-1,-1,-1):

		if dist[j]==0:
			addl+=1
			bool1=True
		else:
			if bool1==True:
				for k in range (j-1,j-10,-1):
					try:
						if dist[k]!=0:
							if memleng[k]==0:
								memleng[k]=(addl-10)/20
					except IndexError:
						pass
				bool1=False
				addl=0
	return [numtm,tail,segment,dist,memleng]

def cut(seq):
	res={}
	for i in range (0, len(seq)):
		res[i]=""
		for j in range (i-5,i+5):
			if 0<=j<len(seq):
				res[i]=res[i]+seq[j]
	return res

def protparam(x):
	we=[]
	iso=[]
	instab=[]
	for i in range (0, len(x)):
		analysed_seq = ProteinAnalysis(x[i])
		mw=str(analysed_seq.molecular_weight())
		ie=str(min(analysed_seq.isoelectric_point(),1000000))
		inst=str(min(analysed_seq.instability_index(),1000000))
		mw=str(analysed_seq.molecular_weight()/2000)
		ie=str(min(analysed_seq.isoelectric_point()/10,1))
		inst=str(min((analysed_seq.instability_index()+100)/500,1))
		we.append(mw)
		iso.append(ie)
		instab.append(inst)
	return [we,iso,instab]
def readAAIndex(x):
	shift=0
	shift2=0
	ABC="ARNDCQEGHILKMFPSTWYV"
	aain={}
	infile=open(x)
	while 1:
		l=infile.readline()
		if l=="":
			infile.close()
			break
		l=l.split()
		if l[0]=="I":
			shift+=1
			aain[shift]={}
			shift2=0
		else:
			for i in range (0, len(l)):
				aain[shift][ABC[shift2]]=float(l[i])
				shift2+=1

	infile.close()
	return aain
def indexing(x,y):
	aain=[]
	for i in range (0, len(x)):

		aa1=0
		aa2=0
		aa3=0
		aa4=0
		aa5=0
		for j in range (0, len(x[i])):
			aa1=aa1+y[1][x[i][j]]
			aa2=aa2+y[2][x[i][j]]
			aa3=aa3+y[3][x[i][j]]
			aa4=aa4+y[4][x[i][j]]
			aa5=aa5+y[5][x[i][j]]
		aain.append([aa1,aa2,aa3,aa4,aa5])
	return aain

def checktop(x):
	if x.find("M")==-1:
		return False
	else:
		return True
def castmx(hit,cons,psmx,top,tmreg,tailreg,segmentreg,dist,thick,surf,mw,ie,ins,aa,seg):
	MX=[]
	for i in range (0, len(top)):
		line=""
		line=line+" #"
		line=line+" "+top[i]
		line=line+" "+str(dist[i])
		line=line+" "+str(thick[i])
		for j in range (0, len(psmx[i])):
			line=line+" "+str(psmx[i][j])
		line=line+" "+str(cons[i])
		line=line+" "+str(min(1,float(hit[i])/100))
		line=line+" "+seg[i]
		line=line+" "+str(tailreg[i])
		line=line+" "+str(segmentreg[i])
		line=line+" "+str(tmreg[i])
		line=line+" "+str(surf[i][1])
		line=line+" "+str(surf[i][2])
		line=line+" "+str(surf[i][0])
		line=line+" "+str(mw[i])
		line=line+" "+str(ie[i])
		line=line+" "+str(ins[i])
		aa[i][0]=min(abs((aa[i][0]+15)/40),1)
		aa[i][1]=min(abs((aa[i][1]+5)/1000),1)
		aa[i][2]=min(abs((aa[i][2]+0)/7500),1)
		aa[i][3]=min(abs((aa[i][3]+0)/3500),1)
		aa[i][4]=min(abs((aa[i][4]+0)/70),1)
		line=line+" "+str(aa[i][0])
		line=line+" "+str(aa[i][1])
		line=line+" "+str(aa[i][2])
		line=line+" "+str(aa[i][3])
		line=line+" "+str(aa[i][4])
		MX.append(line.split())
	return MX

def cnn(INS,ins,outs,OUTS,seq,mx):
	top_seq=[]
	raw_res=[]
	for i in range (0, len(seq)):
		curr=[]
		posi=0
		for j in range (i-5,i+6):
			curr.append([])

			if j>=0 and j<len(mx)-1:
				for k in range (3, len(mx[j])):
					curr[posi].append(float(mx[j][k]))
			else:
				for k in range (3, len(mx[0])):
					curr[posi].append(0)

			posi+=1
		if mx[i][1]=="I" and int(mx[i][2])<15:
			x=np.array([curr])
			x = x.reshape(x.shape[0], 11, 38,1)
			x = x.astype('float32')
			res=ins.predict(x)
			raw_res.append(res[0][1])
			top_seq.append(4)
		elif mx[i][1]=="I" and int(mx[i][2])>=15:
			x=np.array([curr])
			x = x.reshape(x.shape[0], 11, 38,1)#
			x = x.astype('float32')
			res=INS.predict(x)
			raw_res.append(res[0][1])
			top_seq.append(3)
		elif mx[i][1]=="O" and int(mx[i][2])<15:
			x=np.array([curr])
			x = x.reshape(x.shape[0], 11, 38,1)
			x = x.astype('float32')
			res=outs.predict(x)
			raw_res.append(res[0][1])
			top_seq.append(2)
		elif mx[i][1]=="O" and int(mx[i][2])>=15:
			x=np.array([curr])
			x = x.reshape(x.shape[0], 11, 38,1)#
			x = x.astype('float32')
			res=OUTS.predict(x)
			raw_res.append(res[0][1])
			top_seq.append(1)
		else:
			raw_res.append(0)
			top_seq.append(0)
	return [raw_res,top_seq]
def magic(value):
	if value<=0.65:
		value=value*0.76923
	else:
		value=1-((1-value)*1.42857)
	return value
def rescale(raw):
	smooth=[]
	for i in range (0, len(raw)):
		summa=0
		cases=0
		for j in range (i-7,i+7):
			try:
				summa=summa+float(raw[j])
				cases+=1
			except Exception:
				pass
		if raw[i]=="-":
			smooth.append("-")
		else:
			smooth.append(magic(summa/cases))
	return smooth

def bidirlstm(bd,x,y):
	feat=[]
	raw_res=[]
	for i in range (0, len(x)):
		feat.append([])
		for k in range (i-11,i+1):
			if k>=0:
				feat[len(feat)-1].append([x[k],y[k]])
			else:
				feat[len(feat)-1].append([0,0])
	for i in range (0, len(feat)):
		xx=bd.predict(np.array([feat[i]]))
		raw_res.append(xx[0][0])
	return raw_res

def rescale2(raw):
	smooth=[]
	for i in range (0, len(raw)):
		summa=0
		cases=0
		for j in range (i-4,i+4):
			try:
				summa=summa+float(raw[j])
				cases+=1
			except Exception:
				pass
		if raw[i]=="-":
			smooth.append("-")
		else:
			smooth.append(summa/cases)
	return smooth

def fasprint(x,y,z):
	print("#number residue topology score")
	for i in range (0, len(x)):
		print(str(i)+" "+x[i]+" "+y[i]+" "+str(z[i]))

