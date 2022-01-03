import sys
import string
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
#from __future__ import print_function

global periodfs
global maximum
global offset

def press(event):
    global wavenumber
    global periodfs
    global maximum
    global offset
    #print('press', event.key)
    sys.stdout.flush()
    if event.key != '':
        plt.cla()
        #ax.set_title('Key pressed')
        ax.plot(times,angles)
        ax.errorbar(times,angles,xerr=timestd,yerr=yerr)
        ax.set_xlabel('Delay time [fs]')
        ax.set_ylabel('Torsion angle [$^\circ$]')
        newperiodfs=periodfs
        newmaximum=maximum
        newoffset=offset
        if event.key== u'ctrl+right':
            newperiodfs=periodfs+1.0
        if event.key== u'ctrl+left':
            newperiodfs=periodfs-1.0
        if event.key== u'left':
            newmaximum=maximum+5.0
        if event.key== u'right':
            newmaximum=maximum-5.0
        if event.key== u'up':
            newoffset=offset+0.1
        if event.key== u'down':
            newoffset=offset-0.1
        if event.key== 's':
            px=np.linspace(times[0],times[len(times)-1],10000)
            py=newoffset+(amplitude*np.sin( (2*np.pi*px/newperiodfs) + (newmaximum*0.5*np.pi/newperiodfs) ) )
            ax.plot(px,py,'r--')
            period=newperiodfs*1.0e-15
            frequency=1.0/period
            wavenumber=frequency/29979245800
            ax.set_title(title+'; '+str(int(wavenumber[0]))+r' cm$^{-1}$')

            ax.set_title(title+'; '+str(int(wavenumber[0]))+r' cm$^{-1}$')
            filename='C:/Users/Thomas/Desktop/PYMOL/'+title+'.eps'
            plt.savefig(filename)
            
        px=np.linspace(times[0],times[len(times)-1],10000)
        py=newoffset+(amplitude*np.sin( (2*np.pi*px/newperiodfs) + (newmaximum*0.5*np.pi/newperiodfs) ) )
        ax.plot(px,py,'r--')
        period=newperiodfs*1.0e-15
        frequency=1.0/period
        wavenumber=frequency/29979245800
        ax.set_title(title+'; '+str(int(wavenumber[0]))+r' cm$^{-1}$')

        fig.canvas.draw()
        
        periodfs=newperiodfs
        maximum=newmaximum
        offset=newoffset
        
def gettorsion(infile, residue, atoms):
    pdb=open(infile,'r')
    C=np.zeros((4,3))
    for line in pdb:
        #print line
        if not 'ATOM' and not 'HETATM 'in line:
            continue
        if 'LINK' in line:
            continue
        if line[17:26]==residue and (line[13:16] == atoms[0]):
            C[0,0]=float(line[28:38])
            C[0,1]=float(line[38:46])
            C[0,2]=float(line[46:54])
        if line[17:26]==residue and (line[13:16] == atoms[1]):
            C[1,0]=float(line[28:38])
            C[1,1]=float(line[38:46])
            C[1,2]=float(line[46:54])
        if line[17:26]==residue and (line[13:16] == atoms[2]):
            C[2,0]=float(line[28:38])
            C[2,1]=float(line[38:46])
            C[2,2]=float(line[46:54])
        if line[17:26]==residue and (line[13:16] == atoms[3]):
            C[3,0]=float(line[28:38])
            C[3,1]=float(line[38:46])
            C[3,2]=float(line[46:54])


            
    coords=np.transpose(C)
    at1=coords[:,0]
    at2=coords[:,1]
    at3=coords[:,2]
    at4=coords[:,3]

    b1=at2-at1
    b2=at3-at2
    b3=at4-at3

    a= np.cross(b1,b2)

    n1=np.cross(b1,b2)/np.linalg.norm(np.cross(b1,b2))
    n2=np.cross(b2/np.linalg.norm(b2),b3/np.linalg.norm(b3))
    m1=np.cross(n1,b2/np.linalg.norm(b2))

    x=np.dot(n1,n2)
    y=np.dot(m1,n2)

    torsion=-1*np.arctan2(y,x)*180.0/3.14159265359

    return torsion

def sinewave(x,amplitude,period,phase,offset):
    s=offset+amplitude*np.sin(2*3.14159*(x-phase)/period)
    return s

dirname=r'Z:/P17807/EXTRAPOLATE_AUTO/RESULTS_ALL0.30/'
pdbs=[
    r'MbCO_May_Reprocessed_Dark.pdb',
    r'extrapolated_150fs_refined1.pdb_fitted.pdb',
    r'extrapolated_225fs_refined1.pdb_fitted.pdb',
    r'extrapolated_300fs_refined1.pdb_fitted.pdb',
    r'extrapolated_375fs_refined1.pdb_fitted.pdb',
    r'extrapolated_450fs_refined1.pdb_fitted.pdb',
    r'extrapolated_525fs_refined1.pdb_fitted.pdb',
    r'extrapolated_600fs_refined1.pdb_fitted.pdb',
    r'extrapolated_750fs_refined1.pdb_fitted.pdb',
    r'extrapolated_900fs_refined1.pdb_fitted.pdb',
    r'extrapolated_1300fs_refined1.pdb_fitted.pdb'
        ]




#times=[ 0 , 150 , 225 , 300 , 375 , 450 , 525 , 600 , 750 , 900 , 1300]
times=[ 0 , 254, 327, 402, 471, 604, 627, 702, 847, 1001, 1401]
#timestd=[60,60,60,60,60,60,60,60,60,60,60]
timestd=[ 0, 33, 32, 33, 31,35,33, 33, 34, 33, 34]

#times=[240,330,390,430,460,490,530,560,590,630,740]
#timestd=[125,76,48,41,38,37,37,39,42,46,87]    
#residue='R10 A 300'
#residue='TRP A 182'
residue='HIS A  64'
#residue='HIS A  93'
#residue='LEU A  29'
#residue='HEM A 201'
#residue='LYS A 216'
#residue='MET A 118'
#residue='TRP A  86'
#residue='THR A  90'
#residue='ARG A  82'
#residue='ALA A 215'
#residue='ASP A 212'
#residue='TYR A 185'
#residue='VAL A  68'
#residue='PHE A  43'
#atoms=['CBD','CAD','C3D','O1D']
#atoms=['N  ','CA ','CB ','CG '] # VAL CHI1
atoms=['CA ','CB ','CG ','CD2'] # HIS CHI2
#atoms=['C2B','C3B','CAB','CBB'] # HEM VINYL GROUP
#atoms=['CA ','CB ','CG ','CD2'] # PHE CHI2
#atoms=['C12','C13','C14','C15'] # BOND, ISOMERIZING BOND
#atoms=['N  ','CA ','C  ','O  '] #
#atoms=['C18','C5 ','C13','C20'] # RET OMEGA 2
#atoms=['C18','C5 ','C9 ','C19'] # RET OMEGA 1
#atoms=['C20','C13','C9 ','C19'] # RET OMEGA 3
#atoms=['C8 ','C9 ','C10','C11'] # RET OMEGA 4
#atoms=['C7 ','C8 ','C9 ','C19'] # RET OMEGA 5
#atoms=['C7 ','C6 ','C1 ','C16']
#atoms=['C7 ','C6 ','C5 ','C18']
#atoms=['CB ','CG ','CD ','CE ']
#atoms=['CG ','CD ','CE ','NZ ']
#atoms=['CB ','CG ','SD ','CE '] # MET CHI3
#atoms=['CA ','CB ','CG ','SD '] # MET CHI2
#atoms=['N  ','CA ','CB ','CG '] # MET CHI1
#atoms=['CA ','CB ','CG ','CD2'] # PHE CHI2
#atoms=['CA ','CB ','CG ','CD2'] # TRP CHI2
#atoms=['N  ','CA ','CB ','CG '] # TRP CHI1
#atoms=['CA ','CB ','CG ','CD2'] # TYR CHI2
#atoms=['N  ','CA ','CB ','CG '] # TYR CHI1
#atoms=['N  ','CA ','CB ','CG2'] # THR CHI1
#atoms=['C18','C5 ','C6 ','C7 '] # RET head group
#atoms=['CA ','CB ','CG ','CD2'] # TYR CHI2
#atoms=['CA ','CB ','CG ','OD1'] # ASP CHI2
#atoms=['N  ','CA ','CB ','CG '] # ASP CHI1
#atoms=['C20','C13','C14','C15'] # RET business end
#atoms=['CG ','CD ','NE ','CZ ']
title=residue+'-'+atoms[0]+'-'+atoms[1]+'-'+atoms[2]+'-'+atoms[3]
title=title.replace(' ','')
print(title)

angles=[]
for pdb in pdbs:
    angle=gettorsion(dirname+pdb,residue,atoms)
    #angle=angle % 360.0
    #print pdb, angle
    angles.append(angle)
#angles=scipy.signal.detrend(angles)
for value in angles:
    print(value)


yerr=[0,0,0,0,0,0,0,0,0,0,0]

fig, ax = plt.subplots()  
fig.canvas.mpl_connect('key_press_event', press)



ax.plot(times,angles)
ax.errorbar(times,angles,xerr=timestd,yerr=yerr)
ax.set_xlabel('Delay time [fs]')
ax.set_ylabel('Torsion angle [$^\circ$]')

#################GET PERIOD AND WAVENUMBER################
freqs=np.linspace(0.000001,0.1,200000)
xp=np.asarray(times,dtype=np.float64)
yp=np.asarray(angles,dtype=np.float64)
ls=scipy.signal.lombscargle(xp,scipy.signal.detrend(yp),freqs)

freq=freqs[ls==np.max(ls)]
#print 'The maximum is at',freq
period=2*np.pi/freq
print('Period [fs]')
print(period[0])
periodfs=period
period=period*1.0e-15
frequency=1.0/period
wavenumber=frequency/29979245800
print('Wavenumber [cm-1]')
print(wavenumber[0])
print('--------------------------')
##########################################################

offset=np.average(angles)
maximum=np.max(scipy.signal.detrend(yp))
minimum=np.min(scipy.signal.detrend(yp))
amplitude=(maximum-minimum)/2
maxima=np.where(yp==np.max(yp))
m0=int(maxima[0])
print(m0)
maximum= times[m0]

px=np.linspace(times[0],times[len(times)-1],10000)
py=offset+(amplitude*np.sin( (2*np.pi*px/periodfs) + (maximum*0.5*np.pi/periodfs) ) )
ax.plot(px,py,'r--')

print('Amplitude:',amplitude)
print('Offset:',offset)
print('Maximum position:',maximum)


ax.set_title(r'1 $\mu$J: '+title+'; '+str(int(wavenumber[0]))+r' cm$^{-1}$')
print("Left      : shift left")
print("Right     : shift right")
print("Ctrl-Left : decrease period")
print("Ctrl-Right: increase period")
print("Up        : shift up")
print("Down      : shift down")
print("s         : save eps figure")
plt.show()


plt.plot(freqs,ls)
plt.title('Periodogram')
plt.xlabel('Frequency')
plt.ylabel('Lomb-Scargle value')
plt.show()
