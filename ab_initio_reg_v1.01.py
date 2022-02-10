
#import de libraries
import tkinter as tk
from tkinter import font as tkFont
from regquantv2 import *

xScore=310;yScore=10

carreau=0
damier=0
xcol=-1;ycol=-1
MatDim=0

fen1=tk.Tk()
fen1.title('Programme 3D                     ')
fen1.geometry('600x400+50+100')
#On definit la police de caractere:
arial14b = tkFont.Font(family='Arial', size=14,weight='bold')

arialHUGEb = tkFont.Font(family='Arial', size=50, weight='bold')

can1 = tk.Canvas(fen1, width =580, height =380,bg='white')

can1.grid(row=0,column=0,rowspan=2,padx=5,pady=2)

#Dessine la ligne rouge
can1.create_line(0,yScore+10,580,yScore+10,width=6,fill='red')


eps=0.1 #Precision dans la racine carree
scale=170 #taille en pixel de representation
def orbitale_S(t,phi):
    return scale*1/Rsqrt(4*pi,eps)
def orbitale_Pz(t,phi):
    return scale*Rsqrt(3/(4*pi),eps)*Rcos(t)
def orbitale_Py(t,phi):
    return scale*Rsqrt(3/(4*pi),eps)*(Rsin(phi)*Rsin(t))
def orbitale_Px(t,phi):
    return scale*Rsqrt(3/(4*pi),eps)*(Rcos(phi)*Rsin(t))
def orbitale_Dz(t,phi):
    return scale*Rsqrt(3/(8*pi),eps)*(1-3*Rcos(t)**2)
def orbitale_Dy(t,phi):
    return scale*Rsqrt(3/(4*pi),eps)*(Rsin(2*phi)*Rsin(t))
def orbitale_Dx(t,phi):
    return scale*Rsqrt(3/(4*pi),eps)*(Rcos(2*phi)*Rsin(t))

def Radial(r,n,l):
    if l==0:
        if n==1 :  return 2*Rexp(-r)
        if n==2 :  return 1/Rsqrt(2,0.1)*(1-r/2)*Rexp(-r/2)
        if n==3 :  return 2/3/Rsqrt(2)*(1-2*r/3+2*r*r/27)*Rexp(-r/3)
        if n==4 :  return 1/4*(1-3*r/4+r*r/8-r*r*r/192)*Rexp(-r/4)
    if l==1:
        if n==2:   return 2/Rsqrt(6)*r*Rexp(-r/2)
        if n==3:   return 8/27/Rsqrt(6)*r*(1-r/6)*Rexp(-r/3)
        if n==4:   return Rsqrt(5/3)/16*r*(1-r/4+r*r/80)*Rexp(-r/4)
    if l==2:
        if n==3:   return 4/Rfact(8)/Rsqrt(30)*r*r*Rexp(-r/3)
        if n==4:   return 1/64/Rsqrt(5)*r*r*(1-r/12)*Rexp(-r/4)
    if l==3:
        if n==4:   return 1/768/Rsqrt(35)*r*r*r*Rexp(-r/4)

defcan.can=can1 #definir le canevas de tkinter pour l'affichage 3D
d=Obj3D() #On doit mettre comme parametres trois listes vides lors de l'initialisation.
e=Obj3D()
#f=Obj3D()
g=Obj3D()

d.courbe3D((230,200,0),orbitale_Dz,30,20)
e.courbe3D((300,200,0),orbitale_Dx,20,40)
#f.courbe3D((370,200,100),orbitale_Px,20,20)

zeroCube=[((0,0,0),(600,600,600))]
print("ZeroCube=",zeroCube)
g.isosurface((370,200,100),orbitale_Px,zeroCube,10)

c=d+e+f #merge les objets
#Attention, l'objet cree est modifie dans c, car l'objet cree envoie des listes de self. au lieu d'en creer de nouvelles
c.orig((300,200,0)) #Fixer le point de rotation
pi=3.14159265
#mode='plain'
mode='above'  #Affichage vu de dessus

c.draw(mode)


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





