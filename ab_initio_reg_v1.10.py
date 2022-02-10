
#import de libraries
import tkinter as tk
from tkinter import font as tkFont
from regquantv5 import *
from regPlot import *

#utilisation des nombres complexes qui un type dans python
#a=1+1j
#a=1+b*1j
#a.real
#a.imag


xScore=310;yLine=10

xcol=-1;ycol=-1
pi=3.14159265

fen1=tk.Tk()
fen1.title('Programme 3D                     ')
fen1.geometry('800x400+50+100')
#On definit la police de caractere:
arial14b = tkFont.Font(family='Arial', size=14,weight='bold')

arialHUGEb = tkFont.Font(family='Arial', size=50, weight='bold')

can1 = tk.Canvas(fen1, width =780, height =380,bg='white')

can1.grid(row=0,column=0,rowspan=2,padx=5,pady=2)

#Dessine la ligne rouge
can1.create_line(0,yLine,780,yLine,width=6,fill='red')


eps=0.1 #Precision dans la racine carree

def HarmSph(t,phi,l,m):
    if l==0 and m ==0: return 0.2820947918 #scale*1/Rsqrt(4*pi,eps)
    if l==1 and m ==0: return 0.4886025119*Rcos(t) #scale*Rsqrt(3/(4*pi),eps)*Rcos(t)
    if l==1 and m ==1: return 0.3454941495*Rsin(t)*(Rcos(phi)+Rsin(phi)*1j) #scale*Rsqrt(3/(8*pi),eps)*(Rsin(phi)*Rsin(t))
    if l==1 and m ==-1:return 0.3454941495*Rsin(t)*(Rcos(phi)-Rsin(phi)*1j) #scale*Rsqrt(3/(8*pi),eps)*(Rcos(phi)*Rsin(t))
    if l==2 and m ==0: return 0.3153915653*(1-3*Rcos(t)**2) #scale*Rsqrt(5/(16*pi),eps)*(1-3*Rcos(t)**2)
    if l==2 and m ==1: return 0.772548404*Rsin(t)*Rcos(t)*(Rcos(phi)+Rsin(phi)*1j) #scale*Rsqrt(15/(8*pi),eps)*Rsin(t)*Rcos(t)*(Rcos(phi))
    if l==2 and m ==-1: return 0.772548404*Rsin(t)*Rcos(t)*(Rcos(phi)-Rsin(phi)*1j)
    if l==2 and m ==2: return 0.386274202*Rsin(t)**2 * (Rcos(2*phi)+Rsin(2*phi)*1j)
    if l==2 and m ==-2: return 0.386274202*Rsin(t)**2 * (Rcos(2*phi)-Rsin(2*phi)*1j)
    
def Radial(r,n,l):
    eps=0.1
    if l==0:
        if n==1 :  return 2*Rexp(-r)
        if n==2 :  return 0.7071067812*(1-r/2)*Rexp(-r/2) #1/Rsqrt(2,eps)*(1-r/2)*Rexp(-r/2)
        if n==3 :  return 0.3849001795*(1-2*r/3+2*r*r/27)*Rexp(-r/3) #2/(3*Rsqrt(3,eps))*(1-2*r/3+2*r*r/27)*Rexp(-r/3)
        if n==4 :  return 0.25*(1-3*r/4+r*r/8-r*r*r/192)*Rexp(-r/4)
    if l==1:
        if n==2:   return 0.2041241452*r*Rexp(-r/2) #1/(2*Rsqrt(6,eps))*r*Rexp(-r/2)
        if n==3:   return 0.1209624564*r*(1-r/6)*Rexp(-r/3) #8/27/Rsqrt(6,eps)*r*(1-r/6)*Rexp(-r/3)
        if n==4:   return 0.08068715305*r*(1-r/4+r*r/80)*Rexp(-r/4) #Rsqrt(5/3,eps)/16*r*(1-r/4+r*r/80)*Rexp(-r/4)
    if l==2:
        if n==3:   return 9.016009177E-03*r*r*Rexp(-r/3) #4/81/Rsqrt(30,eps)*r*r*Rexp(-r/3)
        if n==4:   return 6.98771243E-03*r*r*(1-r/12)*Rexp(-r/4) #1/64/Rsqrt(5,eps)*r*r*(1-r/12)*Rexp(-r/4)
    if l==3:
        if n==4:   return 1/768/Rsqrt(35,eps)*r*r*r*Rexp(-r/4)

def OrbAtomic(r,theta,phi,n,l,m):
    return HarmSph(theta,phi,l,m)*Radial(r,n,l) #Attention c'est un nombre complexe

def orbitale_Dz(t,phi):
    return OrbAtomic(6,t,phi,3,2,0)
def orbitale_Dx(t,phi):
    return OrbAtomic(6,t,phi,3,2,1)


# ---------------------- Integrales ------------
# Objectif: vérifier la normalisation des orbitales

def Integrale_r(n,l):
    N=100
    r=0
    rmax=50
    dr=rmax/N
    sum=0
    #Methode des trapezes
    sum += r*r**0.5*(Radial(0,n,l)**2+Radial(rmax,n,l)**2 )
    while r < rmax :
        f=Radial(r,n,l)
        sum +=r*r*f*f
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr #Pour l'integration sur phi
    print("\nNormalisation: L'integrale de Radial(",n,",",l,") vaut ",sum)

def calculate():
    a=0
    while a==0:
        str=input("(q=quit) Radial : n,l = ?   Répondre par n,l : ")
        if str=='q': a=1
        if a==0:
            n,l=int(str.split(',')[0]),int(str.split(',')[1])
            Integrale_r(n,l)
#calculate()
#exit()

def CalculateIntegrale(n,l,m):
    #integration par la methode des trapezes    

    def OrbAtomToInt(r,theta,phi):
        return OrbAtomic(r,theta,phi,n,l,m) # Pour n,l,m=2,1,0
        #return HarmSph(theta,phi,2,0)
    
    rmax=25
    N=100
    dr=rmax/N
    dtheta=pi/N*2
    dphi=2*pi/N*2
    r=0
    sum=0
    #Methode des trapezes
    while r < rmax :
        theta=0
        sum += r*r*Rsin(theta)*0.5*(OrbAtomToInt(r,0,0)**2+OrbAtomToInt(r,pi,0)**2)
        while theta+dtheta < pi:
            phi=0
            phisum=0
            while phi < 2*pi:
                f=OrbAtomToInt(r,theta,phi)
                phisum += f*f.conjugate() #on calcule la norme
                phi +=dphi
            sum +=r*r*Rsin(theta)*phisum
            theta += dtheta
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
        #if a % 10 == 0 and a != 0: print(int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr*dtheta*dphi #Pour l'integration sur phi
    print("L'integrale vaut ",sum)




def calculate2():
    a=0
    while a==0:
        str=input("(q=quit)Int(Rnl(r)*HarmSph_lm(theta,phi) : n,l,m = ?   Répondre par n,l,m : ")
        if str=='q': a=1
        if a==0:
            n,l,m=int(str.split(',')[0]),int(str.split(',')[1]),int(str.split(',')[2])
            CalculateIntegrale(n,l,m)

#calculate2() Pour lancer le calcul complet de la normalisation d'une OA


# ---------------------Intégration <psi(n,l,m) | H1e | psi(n,l,m)> ---------------



def DeriveRadial1(r,dr,n,l):
    dr1=dr/2
    return (Radial(r+dr1,n,l)-Radial(r-dr1,n,l))/dr
def DeriveRadial2(r,dr,n,l):
    dr1=dr/2
    return (Radial(r+dr1,n,l)+Radial(r-dr1,n,l) - 2*Radial(r,n,l))/dr1**2
def LaplacienRadial(r,dr,n,l):
    #1/r²dr(r²dr(Rnl))= d2r(Rnl)+2/r.dr(Rnl)
    return DeriveRadial2(r,dr,n,l)+2*DeriveRadial1(r,dr,n,l)/r

def DeriveTheta1(theta,dtheta,phi,l,m):
    dt1=dtheta/2
    return Rsin(theta)*(HarmSph(theta+dt1,phi,l,m)-HarmSph(theta-dt1,phi,l,m))/dtheta
def DeriveTheta2(theta,dtheta,phi,l,m):
    dt1=dtheta/2
    return (DeriveTheta1(theta+dt1,dtheta,phi,l,m) - DeriveTheta1(theta-dt1,dtheta,phi,l,m))/dtheta

def DerivePhi2(theta,phi,dphi,l,m):
    dp1=dphi/2
    return (HarmSph(theta,phi+dp1,l,m) + HarmSph(theta,phi-dp1,l,m) - 2*HarmSph(theta,phi,l,m))/dp1**2


def CalculateBraKetDeH(n,l,m):
    #integration par la methode des trapezes    

    def OrbAtomToInt(r,theta,phi):
        return OrbAtomic(r,theta,phi,n,l,m) # Pour n,l,m=2,1,0
        #return HarmSph(theta,phi,2,0)
    
    rmax=75
    N=125
    dr=rmax/N
    dtheta=pi/N*2
    dphi=2*pi/N*2
    r=0+dr/1000
    sum=0
    #Methode des rectangles
    while r < rmax :
        theta=dtheta/1000
        sum += 0 #r*r*Rsin(theta)*0.5*(OrbAtomToInt(r,0,0)**2+OrbAtomToInt(r,pi,0)**2)
        while theta < pi:
            phi=0
            phisum=0
            while phi < 2*pi:
                f=OrbAtomToInt(r,theta,phi)
                phisum += f.conjugate()*(  (-Rsin(theta))*(0.5*r**2*LaplacienRadial(r,dr/1000,n,l)+r*Radial(r,n,l))*HarmSph(theta,phi,l,m) \
                                           - Radial(r,n,l)*(DeriveTheta2(theta,dtheta/1000,phi,l,m) \
                                           + DerivePhi2(theta,phi,dphi/1000,l,m)/Rsin(theta)) )
                phi +=dphi
            sum +=phisum
            theta += dtheta
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
        #if a % 10 == 0 and a != 0: print(int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr*dtheta*dphi #Pour l'integration sur phi
    print("L'integrale vaut ",sum)


def Calculate3():
    a=0
    while a==0:
        str=input("(q=quit)<psi(n,l,m)|H1e|psi(n,l,m)> : n,l,m = ?   Répondre par n,l,m : ")
        if str=='q': a=1
        if a==0:
            n,l,m=int(str.split(',')[0]),int(str.split(',')[1]),int(str.split(',')[2])
            CalculateBraKetDeH(n,l,m)

#Calculate3()

RMAX=8
NPTS=20
#------------------ Integral de recouvrement pour H2+ ----------
def CalculateRecover(Rz):
    #R est un vecteur R=(0,0,Rz) l'ion H2+ est oriente selon (Oz)
    #integration par la methode des rectangles integrale de recouvrement <PhiA(r)|PhiB(r)
    #Les atomes H sont espaces de R
    
    rmax=75
    N=125
    dr=rmax/N
    dtheta=pi/N*2
    dphi=2*pi/N*2
    r=0+dr/10000
    sum=0
    eps=0.01 #pour le calcul de la racine carree
    #Methode des rectangles
    while r < rmax :
        theta=dtheta/1000
        sum += 0
        while theta < pi:
            phi=0
            phisum=0
            while phi < 2*pi:
                term1=r*Rcos(theta)*Rz
                term2=Rz*Rz/4+r*r
                f= 0.3183098861837907*Rexp(-Rsqrt(term1+term2,eps)-Rsqrt(term2-term1,eps) )
                phisum += f
                phi +=dphi
            sum +=r**2*Rsin(theta)*phisum
            theta += dtheta
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
        #if a % 10 == 0 and a != 0: print(int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr*dtheta*dphi #Pour l'integration sur phi
    print("L'integrale vaut ",sum)
    return sum

def RecoverInt():
    RecInt=[]
    Rmax=RMAX # en Bohr radius
    n=NPTS
    for i in range (0,n):
            RecInt.append(CalculateRecover(i*Rmax/n))
    print("Valeur de l'integrale de recouvrement RecInt =",RecInt)

#RecoverInt()
#(10->4)RecInt = [1.0000218719340805, 0.976579878752984, 0.9109595749407907, 0.8144565567304402, 0.7022784427608751, 0.5907676453530467, 0.4866057752575087, 0.392919248092388, 0.3126366918945206, 0.24565812339091386]
RecInt = [1.0000218719340805, 0.976579878752984, 0.9109595749407907, 0.8144565567304402, 0.7022784427608751, 0.5907676453530467, 0.4866057752575087, 0.392919248092388, 0.3126366918945206, 0.24565812339091386, 0.19067972015076715, 0.14660363666307574, 0.11179188504434401, 0.08454891265140123, 0.06352383658977005, 0.04745143880826213, 0.03523831144198642, 0.02604231113705095, 0.019163751680545973, 0.014040333260325314]

#------------------ Integral J pour H2+ ----------
def CalculateJ(Rz):
    #R est un vecteur R=(0,0,Rz) l'ion H2+ est oriente selon (Oz)
    #integration par la methode des rectangles integrale de recouvrement <PhiA(r)|PhiB(r)
    #Les atomes H sont espaces de R
    
    rmax=75
    N=125
    dr=rmax/N
    dtheta=pi/N*2
    dphi=2*pi/N*2
    r=0+dr/10000
    sum=0
    eps=0.01 #pour le calcul de la racine carree
    #Methode des rectangles
    while r < rmax :
        theta=dtheta/1000
        sum += 0
        while theta < pi:
            phi=0
            phisum=0
            while phi < 2*pi:
                term1=r*Rcos(theta)*Rz
                term2=Rz*Rz/4+r*r
                f= 0.3183098861837907*Rexp(-2*Rsqrt(term1+term2,eps) )/Rsqrt(term2-term1,eps)
                phisum += f
                phi +=dphi
            sum +=r**2*Rsin(theta)*phisum
            theta += dtheta
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
        #if a % 10 == 0 and a != 0: print(int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr*dtheta*dphi #Pour l'integration sur phi
    print("L'integrale vaut ",sum)
    return sum

def JInt():
    JInt=[]
    Rmax=RMAX # en Bohr radius
    n=NPTS
    for i in range (0,n):
            JInt.append(CalculateJ(i*Rmax/n))
    print("Valeur de l'integrale de recouvrement JInt =",JInt)

#JInt()

#JInt = [0.8952139234586354, 0.8647352475901295, 0.7869546639511883, 0.6839799146534887, 0.5606205145631958, 0.47473346967676777, 0.410018419745145, 0.35421609076458077, 0.3123157737843805, 0.2797949691689589]
JInt = [0.8952139234586354, 0.8647352475901295, 0.7869546639511883, 0.6839799146534887, 0.5606205145631958, 0.47473346967676777, 0.410018419745145, 0.35421609076458077, 0.3123157737843805, 0.2797949691689589, 0.2511431387636676, 0.22858563229940063, 0.21043033174803644, 0.19347448343689766, 0.17971657457979395, 0.1683764751396771, 0.15721296528912915, 0.14799423662180272, 0.14029538696519575, 0.13237856749366927]
#------------------ Integral de K pour H2+ ----------
def CalculateK(Rz):
    #R est un vecteur R=(0,0,Rz) l'ion H2+ est oriente selon (Oz)
    #integration par la methode des rectangles integrale de recouvrement <PhiA(r)|PhiB(r)
    #Les atomes H sont espaces de R
    
    rmax=75
    N=125
    dr=rmax/N
    dtheta=pi/N*2
    dphi=2*pi/N*2
    r=0+dr/10000
    sum=0
    eps=0.01 #pour le calcul de la racine carree
    #Methode des rectangles
    while r < rmax :
        theta=dtheta/1000
        sum += 0
        while theta < pi:
            phi=0
            phisum=0
            while phi < 2*pi:
                term1=r*Rcos(theta)*Rz
                term2=Rz*Rz/4+r*r
                f= 0.3183098861837907*Rexp(-Rsqrt(term1+term2,eps)-Rsqrt(term2-term1,eps) )/Rsqrt(term2+term1,eps)
                phisum += f
                phi +=dphi
            sum +=r**2*Rsin(theta)*phisum
            theta += dtheta
        r += dr
        a = int(100*r/rmax)
        if a % 10 == 0 and a != 0: print("Progression: "+int(a/10)*"="+(10-int(a/10))*" "+"|",int(100*r/rmax),r"%",end="\r",flush=True)
        #if a % 10 == 0 and a != 0: print(int(100*r/rmax),r"%",end="\r",flush=True)
    sum *= dr*dtheta*dphi #Pour l'integration sur phi
    print("L'integrale vaut ",sum)
    return sum

def KInt():
    KInt=[]
    Rmax=RMAX # en Bohr radius
    n=NPTS
    for i in range (0,n):
            KInt.append(CalculateK(i*Rmax/n))
    print("Valeur de l'integrale de recouvrement KInt =",KInt)

#KInt()
#KInt = [0.8952139234586354, 0.8683217498020176, 0.798333800596906, 0.703362514677367, 0.5200771208325459, 0.4034270080455887, 0.3204308929629806, 0.2305940665028692, 0.17069083505589583, 0.12927396336430122]
KInt = [0.8952139234586354, 0.8683217498020176, 0.798333800596906, 0.703362514677367, 0.5200771208325459, 0.4034270080455887, 0.3204308929629806, 0.2305940665028692, 0.17069083505589583, 0.12927396336430122, 0.09164541548109945, 0.0662388190883088, 0.04880966789263915, 0.034279352034604145, 0.024417991906849283, 0.0176753639531791, 0.012334398179272076, 0.008700113860284268, 0.006219547508768227, 0.004320084439158528]

# --- Tracer du profil de l'energie en fonction de la distance internucleaire H2+ -----
def TracerDisInterNucl():
    npts=NPTS
    rmin=0.001
    rmax=RMAX
    fct=[] ; fct2=[] ; fct3=[]; fct4=[]
    t0=0 #clock()
    ymin=0
    ymax=0
    for i in range(npts):
        r=rmin+rmax*i/npts
        y1=RecInt[i]
        y2=KInt[i]
        y3=JInt[i]
        y4=(-0.5+1/r) - ( JInt[i] + KInt[i] )/(1+RecInt[i])
        #y2=1/(1-RecInt[i])*((-0.5+1/r)*(1-RecInt[i]) - ( JInt[i] - KInt[i] ))
        if y1<ymin:  ymin=y1
        if y1>ymax : ymax=y1
        if y2<ymin:  ymin=y2
        if y2>ymax : ymax=y2
        if y3<ymin:  ymin=y3
        if y3>ymax : ymax=y3
        if y4*2<ymin:  ymin=y4*2
        
        fct.append( (r,y1) )
        fct2.append( (r,y2) )
        fct3.append( (r,y3) )
        fct4.append( (r,y4*2) )
        
    g=Graphe(can1,500,200,260,180)
    g.box(rmin,rmax,ymin,1.2)
    g.step(1,ymax/10)
    g.trace(fct,'blue',2)
    txt="S(bleue), J(rouge) et K(noire) en fonction de R "
    g.caption(txt,'black','Arial',10)
    g.trace(fct2,'black',2)
    g.trace(fct3,'red',2)
    g.trace(fct4,'green',2)
    
TracerDisInterNucl()



# ---------------------------------------------
scale=1000 #170 #taille en pixel de representation
def implicite(r,theta,phi):
    return scale*OrbAtomic(r/10,theta,phi,3,1,0).real #r/30 pour n,l,m=2,1,0
#n,l,m=3,1,0 est interessant car on a deux changements de signe.

def AfficherorbitaleMoleculaire():
    global pi,c,mode
    defcan.can=can1 #definir le canevas de tkinter pour l'affichage 3D
    #d=Obj3D() #On doit mettre comme parametres trois listes vides lors de l'initialisation.
    #e=Obj3D()
    #f=Obj3D()
    g=Obj3D()

    #g.courbe3D((230,200,0),orbitale_Dz,30,20)
    #e.courbe3D((300,200,0),orbitale_Dx,20,40)
    #f.courbe3D((370,200,100),orbitale_Px,20,20)

    #zeroCube [r,theta,phi  r2,theta2,phi2
    zeroCube=[((1,0,0,0),(2900,pi,2*pi,0))]
    isotrigger=10 #70000 pour 2,1,0   ; 6 ou 7 pour 3,2,0 ;2800 pour 3,1,0  
    g.rmin=6
    g.PhiAnglemax=pi/12
    g.npassmin=4 #Le nombre de pass de réduction du cube de r (augmenter si les cubes sont vides au demarrage
    g.isosurface((300,200,0),implicite,zeroCube,isotrigger)


    c=g#d+e#+f #merge les objets
    #Attention, l'objet cree est modifie dans c, car l'objet cree envoie des listes de self. au lieu d'en creer de nouvelles
    c.orig((300,200,0)) #Fixer le point de rotation
    pi=3.14159265
    #mode="wired"
    #mode='plain'
    mode='above'  #Affichage vu de dessus
    c.draw(mode)

#AfficherorbitaleMoleculaire()

def DrawFonctions():
    # ---------------- Tracer fonction 2D ----------------------
    npts=1000
    xmin=0
    xmax=20
    fct=[];fct2=[]
    t0=0 #clock()
    ymin=0
    ymax=0
    n=3; l=0
    for i in range(npts):
        x=xmin+i*(xmax-xmin)/npts
        y=2000*Radial(x,n,l)
        if y<ymin: ymin=y
        if y>ymax : ymax=y
        fct.append((x,y))
        
    for i in range(npts):
        x=xmin+i*(xmax-xmin)/npts
        y=2000*Radial(x,3,1)
        if y<ymin: ymin=y
        if y>ymax : ymax=y
        fct2.append((x,y))


    g=Graphe(can1,600,200,160,180)
    d=Rexp(xmax-1)+2
    g.box(xmin,xmax,ymin,ymax)
    g.step(2,ymax/10)
    g.trace(fct,'blue',1)
    txt='Tracé de Radial( n='+str(n)+', l='+str(l)+' )'
    g.caption(txt,'black','Arial',10)
    g.trace(fct2,'black',1)

# ---------------------------------------------------------------



    





tagx=can1.create_text(50,30,anchor='w',text='Rx',font=arial14b,activefill='red',fill='blue')
tagy=can1.create_text(80,30,anchor='w',text='Ry',font=arial14b,activefill='red',fill='blue')
tagz=can1.create_text(110,30,anchor='w',text='Rz',font=arial14b,activefill='red',fill='blue')
tagxm=can1.create_text(45,50,anchor='w',text='-Rx',font=arial14b,activefill='red',fill='blue')
tagym=can1.create_text(75,50,anchor='w',text='-Ry',font=arial14b,activefill='red',fill='blue')
tagzm=can1.create_text(105,50,anchor='w',text='-Rz',font=arial14b,activefill='red',fill='blue')
tagAnime=can1.create_text(150,30,anchor='w',text='Anime',font=arial14b,activefill='red',fill='blue')
tagChange=can1.create_text(150,50,anchor='w',text='Change',font=arial14b,activefill='red',fill='blue')

release=0
n=20
def BPRelease(event):
	global Release,n
	Release=1

def BPx(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((pi/n,0,0))
		c.draw(mode)
def BPy(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((0,pi/n,0))
		c.draw(mode)
def BPz(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((0,0,pi/n))
		c.draw(mode)
def BPxm(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((-pi/n,0,0))
		c.draw(mode)
def BPym(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((0,-pi/n,0))
		c.draw(mode)
def BPzm(event):
	global Release,n
	Release=0
	while Release == 0:
		c.clear()
		c.rot((0,0,-pi/n))
		c.draw(mode)

def BPAnime(event):
	p=400
	for i in range(p):
		c.clear()
	#c.rot((pi/n,0,0)) #Rotation	
		if 0 < i < p/3: c.rot((pi/n,2*pi/n,0.5*pi/n)) #Rotation
		elif p/3 <= i < 2*p/3: c.rot((0.5*pi/n,3*pi/n,pi/n)) #Rotation
		else:  c.rot((0.25*pi/n,2*pi/n,0)) #Rotation
		c.draw(mode)
toggleObj=1
def BPChange(event):
	global toggleObj,c
	if toggleObj==-1:
		c.cube((200,200,0),80) #Centre du cube et largeur des côtés
		toggleObj=1
		return
	else:
		c.sphere((200,200,0),100,20,20)
		toggleObj=-1

can1.tag_bind(tagx, '<ButtonPress-1>', BPx)
can1.tag_bind(tagy, '<ButtonPress-1>', BPy)
can1.tag_bind(tagz, '<ButtonPress-1>', BPz)
can1.tag_bind(tagxm, '<ButtonPress-1>', BPxm)
can1.tag_bind(tagym, '<ButtonPress-1>', BPym)
can1.tag_bind(tagzm, '<ButtonPress-1>', BPzm)
can1.tag_bind(tagAnime, '<ButtonPress-1>', BPAnime)
can1.tag_bind(tagChange, '<ButtonPress-1>', BPChange)
can1.bind('<ButtonRelease>',BPRelease) #Pour la selection d'exercice à droite







#boucle qui lance l'affichage de la fenetre avec son contenu
fen1.mainloop()





