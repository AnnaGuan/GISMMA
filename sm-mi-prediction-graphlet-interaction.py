import numpy as np 
import xlrd

fpredict=open('prediction.txt','w')
fnsm=xlrd.open_workbook('sm_no.xlsx','r')
fnmi=xlrd.open_workbook('mi_no.xlsx','r')


def filetomat(filename):
  f=open(filename,'r')
  array=f.readlines()
  lines=[]
  for line in array:
    line=line.strip('\n')
    linelist=[float(a) for a in line.split()]
    lines.append(linelist)
  matrix=np.array(lines)
  f.close()
  return matrix

def ni0(a,i,j):
    ni0=a[i,j]
    return ni0
def ni1(a,i,j):
    ni1=(1-a[i,j])*(np.dot(a[i,:],a[j,:]))
    return ni1
def ni2(a,b,i,j):
    ni2=a[i,j]*(np.dot(a[i,:],b[j,:]))
    return ni2
def ni3(a,b,i,j):
    ni3=a[i,j]*(np.dot(b[i,:],a[j,:]))
    return ni3
def ni4(a,b,i,j):
    ni4=a[i,j]*(np.dot(a[i,:],a[j,:]))
    return ni4
  
def ni5(a,b,i,j):
    ni5=a[i,j]*np.dot((b[i,:]*a[j,:]),np.dot((b[i,:]*b[j,:]),a))
    return ni5
def ni6(a,b,i,j):
    ni6=(1-a[i,j])*np.dot((a[i,:]*a[j,:]),np.dot((b[i,:]*a[j,:]),b))
    return ni6
def ni7(a,b,i,j):
    ni7=(1-a[i,j])*np.dot((a[i,:]*b[j,:]),np.dot((b[i,:]*a[j,:]),a))
    return ni7
def ni8(a,b,i,j):
    ni8=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((b[i,:]*b[j,:]),a))
    return ni8
def ni9(a,b,i,j):
    ni9=(1-a[i,j])*np.dot((a[i,:]*b[j,:]),np.dot((a[i,:]*a[j,:]),b))
    return ni9
def ni10(a,b,i,j):
    ni10=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((b[i,:]*a[j,:]),b))
    return ni10
def ni11(a,b,i,j):
    ni11=(1-a[i,j])*np.dot((a[i,:]*a[j,:]),np.dot((b[i,:]*b[j,:]),a))
    return ni11
def ni12(a,b,i,j):
    ni12=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((a[i,:]*b[j,:]),b))
    return ni12
def ni13(a,b,i,j):
    ni13=a[i,j]*np.dot((b[i,:]*a[j,:]),np.dot((b[i,:]*a[j,:]),b))
    return ni13
def ni14(a,b,i,j):
    ni14=(1-a[i,j])*np.dot((a[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),b))
    return ni14
def ni15(a,b,i,j):
    ni15=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((b[i,:]*a[j,:]),a))
    return ni15
def ni16(a,b,i,j):
    ni16=a[i,j]*np.dot((b[i,:]*b[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni16
def ni17(a,b,i,j):
    ni17=a[i,j]*np.dot((a[i,:]*a[j,:]),np.dot((b[i,:]*b[j,:]),a))
    return ni17
def ni18(a,b,i,j):
    ni18=a[i,j]*np.dot((a[i,:]*a[j,:]),np.dot((b[i,:]*a[j,:]),b))
    return ni18
def ni19(a,b,i,j):
    ni19=(1-a[i,j])*np.dot((a[i,:]*b[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni19
def ni20(a,b,i,j):
    ni20=a[i,j]*np.dot((a[i,:]*a[j,:]),np.dot((a[i,:]*b[j,:]),b))
    return ni20
def ni21(a,b,i,j):
    ni21=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((a[i,:]*b[j,:]),a))
    return ni21
def ni22(a,b,i,j):
    ni22=(1-a[i,j])*np.dot((b[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni22
def ni23(a,b,i,j):
    ni23=(1-a[i,j])*np.dot((a[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni23
def ni24(a,b,i,j):
    ni24=a[i,j]*np.dot((a[i,:]*b[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni24
def ni25(a,b,i,j):
    ni25=a[i,j]*np.dot((a[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),b))
    return ni25
def ni26(a,b,i,j):
    ni26=a[i,j]*np.dot((b[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni26
def ni27(a,b,i,j):
    ni27=a[i,j]*np.dot((a[i,:]*a[j,:]),np.dot((a[i,:]*a[j,:]),a))
    return ni27

r1=541
r2=831
mirna=np.zeros([r1,r1])
smallmo=np.zeros([r2,r2])
assoarr=np.zeros([664,2])
mirna=filetomat('miRNA similarity matrix.txt')  
smallmo=filetomat('SM similarity matrix.txt')  
assoarr=filetomat('SM-miRNA-association.txt') 

adj=np.zeros((r2,r1))
for i in range(assoarr.shape[0]):
    adj[(int(assoarr[i,0])-1),(int(assoarr[i,1])-1)]=1        
    
sm=mirna.copy()
sd=smallmo.copy()
          
eye1=np.eye(r1,r1)
eye2=np.eye(r2,r2)
reye1=1-eye1
reye2=1-eye2
smv=reye1*(1-sm)
sdv=reye2*(1-sd)

GI=np.zeros([r1,r1,28])
sumi=np.zeros([r1,28])

for i in range(r1):
  for j in range(r1):
    if j!=i:
      GI[i,j,0]=ni0(sm,i,j)
      GI[i,j,1]=ni1(sm,i,j)
      GI[i,j,2]=ni2(sm,smv,i,j)
      GI[i,j,3]=ni3(sm,smv,i,j)
      GI[i,j,4]=ni4(sm,smv,i,j)
      GI[i,j,5]=ni5(sm,smv,i,j)
      GI[i,j,6]=ni6(sm,smv,i,j)
      GI[i,j,7]=ni7(sm,smv,i,j)
      GI[i,j,8]=ni8(sm,smv,i,j)
      GI[i,j,9]=ni9(sm,smv,i,j)
      GI[i,j,10]=ni10(sm,smv,i,j)
      GI[i,j,11]=ni11(sm,smv,i,j)
      GI[i,j,12]=ni12(sm,smv,i,j)
      GI[i,j,13]=ni13(sm,smv,i,j)
      GI[i,j,14]=ni14(sm,smv,i,j)
      GI[i,j,15]=ni15(sm,smv,i,j)
      GI[i,j,16]=ni16(sm,smv,i,j)
      GI[i,j,17]=ni17(sm,smv,i,j)
      GI[i,j,18]=ni18(sm,smv,i,j)
      GI[i,j,19]=ni19(sm,smv,i,j)
      GI[i,j,20]=ni20(sm,smv,i,j)
      GI[i,j,21]=ni21(sm,smv,i,j)
      GI[i,j,22]=ni22(sm,smv,i,j)
      GI[i,j,23]=ni23(sm,smv,i,j)
      GI[i,j,24]=ni24(sm,smv,i,j)
      GI[i,j,25]=ni25(sm,smv,i,j)
      GI[i,j,26]=ni26(sm,smv,i,j)
      GI[i,j,27]=ni27(sm,smv,i,j)
sumi=GI.sum(axis=1)

GID=np.zeros([r2,r2,28])
sumid=np.zeros([r2,28])

for i in range(r2):
  for j in range(r2):
    if j!=i:
      GID[i,j,0]=ni0(sd,i,j)
      GID[i,j,1]=ni1(sd,i,j)
      GID[i,j,2]=ni2(sd,sdv,i,j)
      GID[i,j,3]=ni3(sd,sdv,i,j)
      GID[i,j,4]=ni4(sd,sdv,i,j)
      GID[i,j,5]=ni5(sd,sdv,i,j)
      GID[i,j,6]=ni6(sd,sdv,i,j)
      GID[i,j,7]=ni7(sd,sdv,i,j)
      GID[i,j,8]=ni8(sd,sdv,i,j)
      GID[i,j,9]=ni9(sd,sdv,i,j)
      GID[i,j,10]=ni10(sd,sdv,i,j)
      GID[i,j,11]=ni11(sd,sdv,i,j)
      GID[i,j,12]=ni12(sd,sdv,i,j)
      GID[i,j,13]=ni13(sd,sdv,i,j)
      GID[i,j,14]=ni14(sd,sdv,i,j)
      GID[i,j,15]=ni15(sd,sdv,i,j)
      GID[i,j,16]=ni16(sd,sdv,i,j)
      GID[i,j,17]=ni17(sd,sdv,i,j)
      GID[i,j,18]=ni18(sd,sdv,i,j)
      GID[i,j,19]=ni19(sd,sdv,i,j)
      GID[i,j,20]=ni20(sd,sdv,i,j)
      GID[i,j,21]=ni21(sd,sdv,i,j)
      GID[i,j,22]=ni22(sd,sdv,i,j)
      GID[i,j,23]=ni23(sd,sdv,i,j)
      GID[i,j,24]=ni24(sd,sdv,i,j)
      GID[i,j,25]=ni25(sd,sdv,i,j)
      GID[i,j,26]=ni26(sd,sdv,i,j)
      GID[i,j,27]=ni27(sd,sdv,i,j)
sumid=GID.sum(axis=1)
                        
n=adj.sum(axis=1)                     
nd=adj.sum(axis=0)                     

xt=np.zeros([1,28])
for i in range(r2):
  if int(n[i])>1:                         
    normmj=np.zeros([int(n[i]),28])       
    for mj in range(r1):               
      if adj[i,mj]==1:
        nmk=0                        
        for mk in range(r1):          
          if ((adj[i,mk]==1) and (mk!=mj)):
            for k in range(28):
              if sumi[mk,k]!=0:
                normmj[nmk,k]=(GI[mk,mj,k])/(sumi[mk,k])   
            nmk+=1
        xmj=normmj.sum(axis=0)        
        if (xmj.sum(0))!=0:           
          xt=np.vstack((xt,xmj))   

np.delete(xt,0,0)
x=xt.T
v=np.zeros([28,1])
s=np.ones([x.shape[1],1])

xx=np.dot(x,xt)
xx=np.mat(xx)
xxi=np.linalg.pinv(xx)
xxix=np.dot(xxi,x)
v=np.dot(xxix,s)                      

xtd=np.zeros([1,28])
for j in range(r1):
  if int(nd[j])>1:                         
    normdi=np.zeros([int(nd[j]),28])       
    for di in range(r2):               
      if adj[di,j]==1:
        ndk=0                         
        for dk in range(r2):          
          if ((adj[dk,j]==1) and (dk!=di)):
            for k in range(28):
              if sumid[dk,k]!=0:
                normdi[ndk,k]=(GID[dk,di,k])/(sumid[dk,k])   
            ndk+=1
        xdi=normdi.sum(axis=0)       
        if (xdi.sum(0))!=0:          
          xtd=np.vstack((xtd,xdi))   

np.delete(xtd,0,0)
xd=xtd.T
vd=np.zeros([28,1])
sd=np.ones([xd.shape[1],1])

xxd=np.dot(xd,xtd)
xxd=np.mat(xxd)
xxid=np.linalg.pinv(xxd)
xxixd=np.dot(xxid,xd)
vd=np.dot(xxixd,sd)                      


st=np.zeros([r2,r1])
stm=np.zeros([r2,r1])
std=np.zeros([r2,r1])
lists=[]
for i in range(r2):
  xmj=np.array([28])
  normmj=np.zeros([int(n[i]),28])
  for mj in range(r1):
    normdi=np.zeros([int(nd[mj]),28])
    if adj[i,mj]==0:                   
      if int(n[i])!=0:                      
        nmk=0
        for mk in range(r1):
          if adj[i,mk]==1:
            for k in range(28):
              if sumi[mk,k]!=0:
                normmj[nmk,k]=(GI[mk,mj,k])/(sumi[mk,k])
            nmk+=1
        xmj=normmj.sum(axis=0)    
        stm[i,mj]=np.dot(xmj,v)    #score1

      if int(nd[mj])!=0:
        ndk=0
        for dk in range(r2):
          if adj[dk,mj]==1:
            for k in range(28):
              if sumid[dk,k]!=0:
                normdi[ndk,k]=(GID[dk,i,k])/(sumid[dk,k])
            ndk+=1
        xdi=normdi.sum(axis=0)    
        std[i,mj]=np.dot(xdi,vd)    #score2

      st[i,mj]=(stm[i,mj]+std[i,mj])/2

      listst=[st[i,mj],i,mj]
      lists.append(listst)          



lists.sort(reverse=True)            
table1=fnsm.sheet_by_name('sm') 
table2=fnmi.sheet_by_name('mirna')  
for i in range(len(lists)):         
  fpredict.writelines([str(table1.cell(lists[i][1],0).value),'\t',str(table2.cell(lists[i][2],0).value),'\t',str(lists[i][0]),'\n'])        
fpredict.close()

        
        
   
      
          
          
                  
                  
                  
      
