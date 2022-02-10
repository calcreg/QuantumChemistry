
#from random import randrange
import time
#time.sleep(0.1) #dort 0.1 seconde


pi=3.14159265

def Rsqrt(a,eps): #esp is the precision required
	#returns sqrt(a); We use Newton's algorithm
	#x0=u0-f(u0)/f'(u0)
	#we go until dx<=eps
	if a <0:
		print('error trying sqrt of negative number')
		return -1
	if a==0: return 0
	x=a #We start large for x>1 x>sqrtx (but careful for x<1  sqrt x >x )
	if a < 1:
		x=1/a
	x0=x 
	dx=1+eps
	#print('dx=',dx)
	while dx > eps:
		x1=0.5*(x0+x/x0)
		dx=x0-x1 #the suite of the xi is decreasing
		#print('a,x1,dx=',a,x1,dx)
		x0=x1
	if a >= 1: return x0
	else: return 1/x0
def Rcos(theta):
	#print(theta)
	#Only angles in [0,2*pi]
	if theta <0: theta=-theta #Cosine parity
	if theta>2*pi: theta=theta%(2*pi)
	if theta>2*pi: print('Problem! Rcos: theta > 2*pi')
	#a=pi/20
	#Lcos=[(0,1),(a,0.9876),(2*a,0.951),(3*a,0.891),(4*a,0.809),(5*a,0.7071),(6*a,0.5877),(7*a,0.4539),(8*a,0.309),(9*a,0.1564),(10*a,0)]
	a=pi/80
	Lcos=[(0*a,1),(1*a,0.999229036240723),(2*a,0.996917333733128),(3*a,0.993068456954926),
		(4*a,0.987688340595138),(5*a,0.98078528040323),(6*a,0.972369920397676),(7*a,0.962455236453647),
		(8*a,0.951056516295153),(9*a,0.938191335922484),(10*a,0.923879532511287),(11*a,0.908143173825081),
		(12*a,0.891006524188368),(13*a,0.872496007072797),(14*a,0.852640164354092),(15*a,0.831469612302545),
		(16*a,0.809016994374947),(17*a,0.785316930880745),(18*a,0.760405965600031),(19*a,0.734322509435686),
		(20*a,0.707106781186548),(21*a,0.678800745532942),(22*a,0.649448048330184),(23*a,0.619093949309834),
		(24*a,0.587785252292473),(25*a,0.555570233019602),(26*a,0.522498564715949),(27*a,0.488621241496955),
		(28*a,0.453990499739547),(29*a,0.418659737537428),(30*a,0.38268343236509),(31*a,0.346117057077493),
		(32*a,0.309016994374947),(33*a,0.271440449865074),(34*a,0.233445363855905),(35*a,0.195090322016128),
		(36*a,0.156434465040231),(37*a,0.117537397457838),(38*a,0.078459095727845),(39*a,0.0392598157590687),
		(40*a+0.0000000000000008,0)] #ajout d'un peu à cause des valeurs bizarres de Python
	#theta is in radian
	#We define the cosine roughtly
	def minicos(theta):
		for i in range(len(Lcos)-1):
			if Lcos[i][0] <= theta <= Lcos[i+1][0]:
				return Lcos[i][1]+(theta-Lcos[i][0])*(Lcos[i+1][1]-Lcos[i][1])/(Lcos[i+1][0]-Lcos[i][0])
		print('problème dans minicos pour theta=',theta)
		return 0

	if 0 <= theta <= pi/2 :
		return minicos(theta)
	if pi/2 <= theta <= pi :
		return -minicos(pi-theta)
	if pi <= theta <= 3*pi/2 :
		return -minicos(theta-pi)
	if 3*pi/2 <= theta <= 2*pi :
		return minicos(2*pi-theta)

def Rsin(theta):
	return Rcos(-pi/2+theta)


# -----------------------------------------------------------------------
class defcan:
        can=0  #tkinter canvas used

#Il faut introduire la couleur selon un flag posé que l'on ajuste avec le signe, pouvoir merger des objets par +
class Obj3D:
    def __init__(self):
    #ListDePts est la liste contenant les points bien ordonnees lors de la creation d'un objet particulier
    #ListDeFaces contient les parallélogrammes à tracer, bien ordonnes.
            #Nature est 'sphere' 'cube' 'cylinder' etc...
            self.nature='unknown object'
            self.handle=[]
    #Version précédente: mydisp=Display(can1,(200,200,500),(200,-400,-1000)) #Canevas d'affichage et point de fuite
            self.ptf=(200,200,500)
            self.can=defcan.can
            self.light=(200,-400,-1000)
            self.faces=[]
            self.pts=[]
            self.ListColor=[]
            self.color1=1
            self.color2=2

    def __add__(self,other):
            obj=Obj3D()
            obj.pts=self.pts+other.pts
            obj.faces=self.faces+other.faces
            obj.ListColor=self.ListColor+other.ListColor
            obj.nature="merged object"
            return obj
    
    def orig(self,Origine):
            self.p0=Origine #Utile pour la rotation autour de self.p0 = (x0,y0,z0)
    
    def display(self,canevas,PointDeFuite,PointLight): #Si on veut changer les parametres d'affichage.
            #PointDeFuite est donné par ses coordonnées (x,y,z)
            #PointLight est le 3-uplet de coordonnées de la source lumineuse
            self.ptf=PointDeFuite
            self.can=canevas
            self.light=PointLight

    def color(self,col1,col2):
            self.color1=col1
            self.color2=col2

    def courbe3D(self,Origine,fonct,Nx,Ny):
            self.nature='object created by a curve'
            self.p0=Origine
            self.ntheta=Nx #theta
            self.nphi=Ny #phi
            pts=[]
            pi=3.14159265
            (x0,y0,z0)=self.p0
            #Si on n rectangles en latitude alors il y a n+1 points
            #par contre il y a n points en longitude
            for i in range(self.nphi):
                    anglePhi=2*i*pi/self.nphi
                    cosphi=Rcos(anglePhi)
                    sinphi=Rsin(anglePhi)
                    for j in range(self.ntheta+1):
                            angleTheta=(j+0.3)*pi/(self.ntheta)#On rajouote 0.3 pour fermer le dessous
                            costheta=Rcos(angleTheta)
                            sintheta=Rsin(angleTheta)
                            r=fonct(angleTheta,anglePhi)
                            y=y0+abs(r)*costheta
                            x=x0+abs(r)*sintheta*cosphi
                            z=z0-abs(r)*sintheta*sinphi
                            pts.append((x,y,z))
                            if j < self.ntheta : #Attention il ne faut pas aller jusque ntheta+1
                                            if r < 0 : self.ListColor.append(self.color1)
                                            else: self.ListColor.append(self.color2) 
            self.pts=pts #Attention il y a (Ntheta+1)*Nphi points pour Ntheta*Nphi rectangles
            facelist=[]
            for i in range(self.nphi):
                    for j in range(self.ntheta):
                            facelist.append(self.pts[i*(self.ntheta+1)+j])
                            facelist.append(self.pts[i*(self.ntheta+1)+j+1])
                            facelist.append(self.pts[((i+1)%self.nphi)*(self.ntheta+1)+j+1])
                            facelist.append(self.pts[((i+1)%self.nphi)*(self.ntheta+1)+j])
            self.faces=facelist
            #print('ListColor 0 =',self.ListColor)

    def sphere(self,Origine,Radius,Nx,Ny):
            self.nature='sphere'
            self.p0=Origine
            self.r=Radius
            self.ntheta=Nx #theta
            self.nphi=Ny #phi
            self.ListColor.clear()
            pts=[]
            pi=3.14159265
            (x0,y0,z0)=self.p0
            #Si on n rectangles en latitude alors il y a n+1 points
            #par contre il y a n points en longitude
            for i in range(self.nphi):
                    cosphi=Rcos(2*i*pi/self.nphi)
                    sinphi=Rsin(2*i*pi/self.nphi)
                    for j in range(self.ntheta+1):
                            costheta=Rcos((j+0.3)*pi/(self.ntheta))
                            sintheta=Rsin((j+0.3)*pi/(self.ntheta))
                            y=y0+self.r*costheta
                            x=x0+self.r*sintheta*cosphi
                            z=z0-self.r*sintheta*sinphi
                            pts.append((x,y,z))
                            self.ListColor.append(self.color1)
            self.pts=pts #Attention il y a (Ntheta+1)*Nphi points pour Ntheta*Nphi rectangles
            facelist=[]
            for i in range(self.nphi):
                    for j in range(self.ntheta):
                            facelist.append(self.pts[i*(self.ntheta+1)+j])
                            facelist.append(self.pts[i*(self.ntheta+1)+j+1])
                            facelist.append(self.pts[((i+1)%self.nphi)*(self.ntheta+1)+j+1])
                            facelist.append(self.pts[((i+1)%self.nphi)*(self.ntheta+1)+j])
            self.faces=facelist
            #print('Nombre de rectangles=',len(self.faces))

    def cube(self,Origine,width):
            self.nature='cube'
            self.p0=Origine
            self.w=width
            #Coordonnées du cube
            a=self.w
            (x0,y0,z0)=self.p0
            self.pts=[(-a+x0,-a+y0,a+z0),(a+x0,-a+y0,a+z0),(a+x0,a+y0,a+z0),(-a+x0,a+y0,a+z0),
                    (-a+x0,-a+y0,-a+z0),(a+x0,-a+y0,-a+z0),(a+x0,a+y0,-a+z0),(-a+x0,a+y0,-a+z0)]
            #Il faut respecter l'orientation des rectangles, normale vers l'exterieur
            faceList=[self.pts[0],self.pts[1],self.pts[2],self.pts[3],
                            self.pts[7],self.pts[6],self.pts[5],self.pts[4],
                            self.pts[4],self.pts[5],self.pts[1],self.pts[0],
                            self.pts[3],self.pts[2],self.pts[6],self.pts[7],
                            self.pts[0],self.pts[3],self.pts[7],self.pts[4],
                            self.pts[2],self.pts[1],self.pts[5],self.pts[6]]
            self.faces=faceList

    def rot(self,angle):
            #attention angle est un 3-uplet
            (tx,ty,tz)=(angle[0],angle[1],angle[2])
            costx=Rcos(tx)
            sintx=Rsin(tx)
            costy=Rcos(ty)
            sinty=Rsin(ty)
            costz=Rcos(tz)
            sintz=Rsin(tz)
            (x0,y0,z0)=self.p0 #Centre de rotation
            
            if tx != 0 :
                newfaces=[]
                for i in range(len ( self.faces ) ):
                    (x,y,z)=self.faces[i]
                    dy=y-y0
                    dz=z-z0
                    x1 =  x
                    y1 =  y0 + dy*costx  + dz*sintx
                    z1 =  z0 - dy*sintx  + dz*costx
                    newfaces.append((x1,y1,z1))
                self.faces.clear()
                self.faces=newfaces

            if ty != 0:
                newfaces2=[]
                for i in range(len ( self.faces ) ):
                    (x,y,z)=self.faces[i]
                    dx=x-x0
                    dz=z-z0
                    x1 =  x0 + dx*costy  + dz*sinty
                    y1 =  y
                    z1 =  z0 - dx*sinty  + dz*costy
                    newfaces2.append((x1,y1,z1))
                self.faces.clear()
                self.faces=newfaces2

            if tz != 0:
                newfaces3=[]                    
                for i in range(len ( self.faces ) ):
                    (x,y,z)=self.faces[i]
                    dx=x-x0
                    dy=y-y0
                    x1 =  x0 + dx*costz  + dy*sintz
                    y1 =  y0 - dx*sintz  + dy*costz
                    z1 =  z
                    newfaces3.append((x1,y1,z1))
                self.faces.clear()
                self.faces=newfaces3


    def draw(self,mode):
            #différents mode: 'wired' ; 'plain'
            canevas=self.can
            if self.nature == 'unknown object':
                    print("Object should be named")
                    return
            def Conv3Dto2D(u):
                    (x,y,z)=u[0],u[1],u[2]
                    (xf,yf,zf)=self.ptf #On récupère le point de fuite
                    return (x+(z/zf)*(xf-x),y+(z/zf)*(yf-y))
            def orderlist(liste):
                    #la list est conçu ainsi:
                    #[(x0,y0,z0),(x1,y1,z1), ... ,(xn,yn,zn)]
                    #On ordonne selon le z
                    #Plus tard on fera selon le produit scalaire avec le vecteur lié au point de fuite le plus grand
                    zlist=[]; ordlist=[] ; colorTmpList=[];
                    for i in range(int(len(liste)/4)):
                        (x0,y0,z0)=liste[4*i]
                        (x1,y1,z1)=liste[4*i+1]
                        (x2,y2,z2)=liste[4*i+2]
                        (x3,y3,z3)=liste[4*i+3]
                        #xm=int((x0+x1+x2+x3)/4)
                        #ym=int((y0+y1+y2+y3)/4)
                        zm=int((z0+z1+z2+z3)/4)
                        zlist.append((zm,i,self.ListColor[i])) #On va ordonner par rapport à zm
                    #print(self.ListColor)
                    slist=sorted(zlist)
                    #print(slist)
                    for i in range(int(len(liste)/4)):
                        for j in range(4):
                            ordlist.append(  (int(liste[4*slist[i][1]+j][0]),
                                              int(liste[4*slist[i][1]+j][1]),
                                              int(liste[4*slist[i][1]+j][2]),
                                              self.ListColor[4*slist[i][2]+j] ) ) #On ajoute la couleur dans le p-uplet
                    #print("ordlist (size = ", len(ordlist)," )= ",ordlist)
                    slist.clear()
                    zlist.clear()
                    ordlist.reverse()
                    return ordlist
            #if self.nature == 'cube' and mode == 'wired':
            if mode == 'wired':
                    for i in range(int(len(self.faces)/4)):
                            x0,y0,z0 = self.faces[4*i][0],self.faces[4*i][1],self.faces[4*i][2]
                            x1,y1,z1 = self.faces[4*i+1][0],self.faces[4*i+1][1],self.faces[4*i+1][2]
                            x2,y2,z2 = self.faces[4*i+2][0],self.faces[4*i+2][1],self.faces[4*i+2][2]
                            ux,uy,uz=(x1-x0,y1-y0,z1-z0)
                            vx,vy,vz=(x2-x0,y2-y0,z2-z0)
                            #produit vectoriel
                            wx=uy*vz-uz*vy
                            wy=uz*vx-ux*vz
                            wz=ux*vy-uy*vx
                            #scalaire le vecteur lié au point de fuite
                            (xf,yf,zf)=self.ptf #On récupère les coordonnées du point de fuite
                            Fx,Fy,Fz= (xf-x0 , yf-y0, -zf-z0) #! c'est -zf car on récupère le vecteur qui regarde par dessus
                            prod_scal=wx*Fx+wy*Fy+wz*Fz
                            #Fonctionne bien si zf est de l'ordre de plusieurs milliers
                            for j in range(4):
                                    #(x0,y0,z0)=self.faces[4*i][j]
                                    #(x1,y1,z1)=self.faces[4*i][(j+1)%4] #Le dernier boucle sur le premier
                                    #print (self.faces[4*i+j])
                                    (x0,y0)=Conv3Dto2D( self.faces[4*i+j] )
                                    (x1,y1)=Conv3Dto2D( self.faces[4*i+(j+1)%4] )
                                    #print("4*i+j=",4*i+j," et 4*i+(j+1)%4=",4*i+(j+1)%4,"line",(x0,y0),"-",(x1,y1))
                                    if prod_scal > 0: self.handle.append(canevas.create_line(x0,y0,x1,y1,width=2,fill='black'))
                                    else: self.handle.append(canevas.create_line(x0,y0,x1,y1,width=2,fill='#aaaaaa'))
            #if self.nature == 'cube' and mode == 'plain':
            if mode == 'plain':
                    #print('in order')
                    ordlist=orderlist(self.faces)
                    #print('out order')
                    color=['black','#111111','#222222','#333333','#444444','#555555','#666666','#777777','#888888','#999999','#aaaaaa','#bbbbbb','#cccccc','#dddddd','#eeeeee','white']
                    for i in range(int(len(ordlist)/4)):
                            (x0,y0,z0)=ordlist[4*i]
                            (x1,y1,z1)=ordlist[4*i+1]
                            (x2,y2,z2)=ordlist[4*i+2]
                            (x3,y3,z3)=ordlist[4*i+3]
                            xm=int((x0+x1+x2+x3)/4)
                            ym=int((y0+y1+y2+y3)/4)
                            zm=int((z0+z1+z2+z3)/4)
                            ux,uy,uz=x1-x0,y1-y0,z1-z0
                            vx,vy,vz=x3-x0,y3-y0,z3-z0
                            #Redéfinition de x0,y0,...
                            (X0,Y0)=Conv3Dto2D( ordlist[4*i] )
                            (x1,y1)=Conv3Dto2D( ordlist[4*i+1] )
                            (x2,y2)=Conv3Dto2D( ordlist[4*i+2] )
                            (x3,y3)=Conv3Dto2D( ordlist[4*i+3] )
                            pts=[(X0,Y0),(x1,y1),(x2,y2),(x3,y3),(X0,Y0)]
                            #produit vectoriel
                            wx=uy*vz-uz*vy
                            wy=uz*vx-ux*vz
                            wz=ux*vy-uy*vx
                            NormW=Rsqrt(wx*wx+wy*wy+wz*wz,0.01)
                            #vecteur normal au rectangle
                            if NormW==0:
                                    #print('NormW=0. Set to 1.')
                                    NormW=1
                            wx=wx/NormW
                            wy=wy/NormW
                            wz=wz/NormW
                            #scalaire le vecteur lié au point de fuite
                            (xl,yl,zl)=self.light #On récupère les coordonnées du point de fuite	
                            Lx,Ly,Lz= (xl-xm , yl-ym, zl-zm)
                            NormL=Rsqrt(Lx*Lx+Ly*Ly+Lz*Lz,0.01)
                            #NormL=1
                            if NormL==0:
                                    #print('NormL=0. Set to 1.')
                                    NormL=1
                            Lx=Lx/NormL
                            Ly=Ly/NormL
                            Lz=Lz/NormL
                            prod_scal=wz
                            #prod_scal=-(Lz-abs(Lx*wx+Ly*wy+Lz*wz)*wz) #Produit scalaire avec -ez de L+2w
                            if prod_scal < 0: prod_scal=0 
                            #print('prod_scal=',prod_scal*20)
                            col=int(15*prod_scal)
                            if col >= 15: col=15
                            Outcol=col-1
                            if Outcol < 0: Outcol=0
                            #print(col)
                            self.handle.append(canevas.create_polygon(pts,width=2,outline=color[Outcol],fill=color[col]))
                                    #else: self.handle.append(canevas.create_line(x0,y0,x1,y1,width=2,fill='#aaaaaa'))
            if mode == 'allwhite':
                    ordlist=orderlist(self.faces)
                    #print("ordlist (size = ", len(ordlist)," )= ",ordlist)
                    
                    color1=['black','#111111','#222222','#333333','#444444','#555555','#666666','#777777','#888888','#999999','#aaaaaa','#bbbbbb','#cccccc','#dddddd','#eeeeee','white']
                    color2=['black','#111111','#222222','#222233','#333344','#333355','#444466','#444477','#444488','#505099','#5555aa','#6060bb','#6565cc','#7070dd','#7777ee','#8888FF']
                    color3=['black','#111111','#222222','#223322','#334433','#335533','#446644','#447744','#448844','#559955','#55aa55','#55bb55','#66cc66','#66dd66','#66ee66','#66FF66']

                    for i in range(int(len(ordlist)/4)):
                            (x0,y0,z0,cl1)=ordlist[4*i]
                            (x1,y1,z1,cl2)=ordlist[4*i+1]
                            (x2,y2,z2,cl3)=ordlist[4*i+2]
                            (x3,y3,z3,cl4)=ordlist[4*i+3]
                            if (cl1+cl2+cl3+cl4) == 4*self.color1: #Sur un parallelogramme ce n'est pas toujours la meme couleur
                              color=color1
                            else:
                              color=color2
                            xm=int((x0+x1+x2+x3)/4)
                            ym=int((y0+y1+y2+y3)/4)
                            zm=int((z0+z1+z2+z3)/4)
                            ux,uy,uz=x1-x0,y1-y0,z1-z0
                            vx,vy,vz=x3-x0,y3-y0,z3-z0
                            #Redéfinition de x0,y0,...
                            (X0,Y0)=Conv3Dto2D( ordlist[4*i] )
                            (x1,y1)=Conv3Dto2D( ordlist[4*i+1] )
                            (x2,y2)=Conv3Dto2D( ordlist[4*i+2] )
                            (x3,y3)=Conv3Dto2D( ordlist[4*i+3] )
                            pts=[(X0,Y0),(x1,y1),(x2,y2),(x3,y3),(X0,Y0)]
                            #produit vectoriel
                            wx=uy*vz-uz*vy
                            wy=uz*vx-ux*vz
                            wz=ux*vy-uy*vx
                            NormW=Rsqrt(wx*wx+wy*wy+wz*wz,0.01)
                            #vecteur normal au rectangle
                            if NormW==0:
                                    #print('NormW=0. Set to 1.')
                                    NormW=1
                            wx=wx/NormW
                            wy=wy/NormW
                            wz=wz/NormW
                            #scalaire le vecteur lié au point de fuite
                            (xl,yl,zl)=self.light #On récupère les coordonnées du point de fuite	
                            Lx,Ly,Lz= (xl-xm , yl-ym, zl-zm)
                            NormL=Rsqrt(Lx*Lx+Ly*Ly+Lz*Lz,0.01)
                            #NormL=1
                            if NormL==0:
                                    #print('NormL=0. Set to 1.')
                                    NormL=1
                            Lx=Lx/NormL
                            Ly=Ly/NormL
                            Lz=Lz/NormL
                            prod_scal=wz
                            #prod_scal=-(Lz-abs(Lx*wx+Ly*wy+Lz*wz)*wz) #Produit scalaire avec -ez de L+2w
                            if prod_scal < 0: prod_scal=-prod_scal 
                            #print('prod_scal=',prod_scal*20)
                            col=int(15*prod_scal)
                            if col >= 15: col=15
                            Outcol=col-1
                            if Outcol < 0: Outcol=0
                            #print(col)
                            self.handle.append(canevas.create_polygon(pts,width=2,outline=color[Outcol],fill=color[col]))
                                    #else: self.handle.append(canevas.create_line(x0,y0,x1,y1,width=2,fill='#aaaaaa'))
            canevas.update()
            time.sleep(0.01)

    def clear(self): #nettoie l'affichage mais ne libere pas la memoire de l'objet
            canevas=self.can
            for hdle in self.handle:
                    if hdle!=0: canevas.delete(hdle)
            self.handle.clear()
            #Il ne faut pas updater le canevas sinon effet de sintillement 
            #avant le réaffichage.
            #canevas.update() 






