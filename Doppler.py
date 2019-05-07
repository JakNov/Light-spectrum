import time
import numpy
import cv2
import imageio
import string
import os
#importované funkce, které jsou popsané v sekci 4.1 
from prfc import V
from prfc import Vv
from prfc import pro
from prfc import inter
from prfc import scale
from prfc import backscale
from prfc import VV

#Soustava souřadnic zdroje:
#Počátek [0,0,0] je v levém horním rohu obrázku (tady začíná
# Python obrázek načítat).
#Osa x je orientovaná směrem k dolnímu okraji obrázku. Osa y je
# orientovaná směrem
#k pravému okraji obrázku. Osa z je orientovaná kolmo na obrázek
# ve směru k pozorovateli.


#Jednotky jsou voleny tak, že c = 1; jenotky v rovině xy lze
# změnit konstantou pp na pixel/pp
Vmax = [0,0,-0.5] 
Vmin = [0,0.05,0]
Start = [0.5,0.5,500] #Počáteční poloha pozorovatele; x,y
# souřadnice jsou násobky výšky a šířky obrázku; souřadnice z je
# pak vzdálenost od obrázku

pp = 1000 

Accel = [0,0,-0.1] #Zrychlení, se kterým se mění rychlost, pro
# kterou je obrázek přepočítán
#snahou není simulovat pohyb obecně zrychlujícího pozorovatele, 
#ale jen napočítat dopplerovské posunutí pro více rychlostí 

Frames = 100 #Počet snímků, které budou vytvořeny

time = (Vmax[2] - Vmin[2])/Accel[2] #Doba, po kterou se 
#pozorovatel bude pohybovat je momentálně určena jen podle složky
# rychlosti v ose z
t=time/Frames #Doba, která uplyne mezi dvěma snímky
jpg = '.jpg' #obrázky budou ukládány ve formátu jpg

images = []
Pic = r'Picture.jpg' #Adresa obrázku, se kterým se bude pracovat
img =cv2.imread(Pic) #Načtení původního obrázku






images = images + [img] #Zde se ukládají jednotlivé snímky;
# z tohoto seznamu lze pak například vytvořit animaci přímo
# v Pythonu

#Tato fáze je popsaná v předchozí sekce o funkci Spectr

num_lines = sum(1 for line in open(r'Lambda.txt'))

file = open(r'Lambda.txt')


Lambda = []
Min = 420
Max = 703
for i in range(num_lines): 
    numbers = file.readline()
    numbers = numbers.split()
    for j in range(3):
        numbers[j] = float(numbers[j])
   
    Lambda.append(numbers)
file.close()

num_lines = sum(1 for line in open(r'Spektrum.txt'))

file = open(r'Spektrum.txt')

Spektrum = []
for i in range(num_lines): 
    numbers = file.readline()
    numbers = numbers.split()
    for j in range(3):
        numbers[j] = float(numbers[j])
   
    Spektrum.append(numbers)
file.close()

Min2 = 420
Max2 = 700
Krok = 14

LS = len(Spektrum)
L = len(Lambda)
col = []

height, width, channels = img.shape
stred = [height/2,width/2]

col=[[[[0 for x in range(4)] for ii in range(LS+2)] for j in range(width)] for i in range(height)] 
    
for i in range(height):
    for j in range(width):
        print(i,j) #Výpis pixelu, který je načítán    
                
        B = images[0][i,j,0]
        G = images[0][i,j,1]
        R = images[0][i,j,2]
        
        souc = int(R)+int(B)+int(G)
        if  souc > 0:
            R = scale(R)
            G = scale(G)
            B = scale(B)
            
            col[i][j][LS + 1][0] = 1 #Pokud bude pixel černý,
# hodnota bude změněna na 0
            X = (0.4124564*R + 0.3575761*G + 0.1804375*B)
            Y = (0.2126729*R + 0.7151522*G + 0.0721750*B)
            Z = (0.0193339*R + 0.1191920*G + 0.9503041*B)
            
            x = X /(X + Y + Z)
            y = Y /(X + Y + Z)
            lum = Y/LS

            for ii in range(LS):
                A = [Spektrum[ii][1],Spektrum[ii][2]]
                XX = [x,y]
                    
                for jj in range(L - 2):
                    index = (ii*Krok + jj + 1)%L
                    C = [Lambda[index][1],Lambda[index][2]]
                    D = [Lambda[(index + 1)%L][1],Lambda[(index + 1)%L][2]]
                
                    P = inter(A,XX,C,D)

                    if V(C,P) + V(D,P) - V(C,D) < 0.0000000001:
                                           
                        KK = A[1]*(P[0]-XX[0])/(P[1]*(XX[0] - A[0]))
                         
                        col[i][j][ii][0] = lum * KK / (KK + 1)
                        col[i][j][ii][1] = index
                        KK = C[1]*(D[0] - P[0])/(D[1]*(P[0] - C[0]))

                        col[i][j][ii][2] = (lum - col[i][j][ii][0])* KK / (KK + 1)
                        col[i][j][ii][3] = lum-col[i][j][ii][0]- col[i][j][ii][2]

                        jj = L - 2
                
        else:
            #Pixel byl černý
            col[i][j][LS + 1][0] = 0

            
#Zadání tří bodů, kterými je dán prostor sRGB; viz sekci 1.3
sR = [0.64 , 0.33] 
sG = [0.30 , 0.60]
sB = [0.15 , 0.06]
#Vdálenost bodů trojúhelníku, který definuje sRGB prostor
RG = V(sR,sG)
GB = V(sG,sB)
RB = V(sR,sB)
#Plocha trojůhelníku, který vymezuje sRGB prostor; bude
# zapotřebí, až budeme provádět projekce barev do tohoto prostoru
S_tr0 = RG * Vv(sR,sG,sB)/2


for w in range(Frames):
    #Výpočet složek rychlosti pro daný snímek
    Vx[0] = Vmin[0] + Accel[0]*w*t
    Vx[1] = Vmin[1] + Accel[1]*w*t
    Vx[2] = Vmin[2] + Accel[2]*w*t

    #Nastavení jména daného snímku
    name = r'Doppler_effect'
    name = name + str(w) + jpg

    for i in range(height):
        for j in range(width):
            
            Y = 0  #Intenzita výsledné barvy v prostoru xyY                                                         
            if col[i][j][LS + 1][0] == 1: #Pokud pixel nebyl 
#černý
            
                L_citX = 0
                L_citY = 0
                L_jmen = 0                                                                                              
                #Výpočet poměru vlnových délek signálu 
#reprezentovaných v soustavě pozorovatele P a zdroje Z
                Z =[i,j,0] #Poloha zdroje; jedná se o daný pixel
# na obrázku. Ten  je dán jako rovina z = 0. 

                #Nyní je dopočítána poloha pozorovatele
# v soustavě zdroje.
                P = [Start[0]*height +t*Vx[0]+Accel[0]*t**2/ 2,Start[1]*width + t*Vx[1] + Accel[1]*t**2 / 2,Start[2] + t*Vx[2] + Accel[2]*t**2 / 2]
                #Výpočet vlnového vektoru v soustavě zdroje;
# je dán jako rozdíl poloh pozorovatele a zdroje 
                kv = [P[0] - Z[0],P[1] - Z[1],P[2] - Z[2]]
                #Pokud je vektor rychlosti nebo vlnový vektor
# nulový, kosinus úhlu ve vztahu 3.6 je 0
                if VV(Vx) == 0 or VV(kv) ==0:
                    cost = 0
                else:
                #Kosinus [hlu theta; vztah 3.7
                    cost = (Vx[0]*kv[0]+Vx[1]*kv[1]+Vx[2]*kv[2]) / (VV(Vx)*VV(kv))

                #Signál v našem případě není ve frekvenční 
#reprezentaci, ale je vyjádřen pomocí vlnové délky. Jedná se
# o výraz 3.6 
                K = (numpy.sqrt(1 - VV(Vx)**2))/(1 - VV(Vx)*cost)

                for ii in range(LS):
                    #Dopplerovské posunutí vlnových délek
                    Ind2 = int((col[i][j][ii][1] + Min)*K) - Min
#Vlnová délka, jejíž hodnota je uložena na pozici [i][j][ii][1]
# pole col je dopplerovksy posunuta
                    Ind3 = int(((col[i][j][ii][1] + 1)%L + Min)*K) - Min #Posunutí vedlejší vlnové délky
                    Ind = col[i][j][ii][1]

                    #Pokud se posunuté vlnové délky nacházejí
# ve viditelném spektru, je proveden výpočet 1.1 smíchání barev 
                    if K * Spektrum[ii][0] >= Min2 and K * Spektrum[ii][0] <= Max2:
                        L_citX =L_citX+(Lambda[int((ii*Krok+ Min2)*K) - Min2][1]/Lambda[int((ii*Krok +Min2)*K)- Min2][2])*col[i][j][ii][0] #Jsou započteny intenzity
# původních barev, ty však mají nyní nové souřadnice odpovídající
# posunutým barvám
                        L_citY = L_citY + col[i][j][ii][0]
                        L_jmen = L_jmen + col[i][j][ii][0]/Lambda[int((ii*Krok + Min2)*K) - Min2][2]
                        Y = Y + col[i][j][ii][0]
                       
                                                                                                             
                    if K * Lambda[Ind][0] >= Min and K * Lambda[Ind][0] <= Max:
                        L_citX = L_citX +  (Lambda[Ind2][1]/Lambda[Ind2][2])*col[i][j][ii][2]
                        L_citY = L_citY + col[i][j][ii][2]
                        L_jmen = L_jmen + col[i][j][ii][2]/Lambda[Ind2][2]
                        Y = Y + col[i][j][ii][2]
                      
                    if K * Lambda[(Ind + 1)%L][0] >= Min and K * Lambda[(Ind + 1)%L][0] <= Max:
                        L_citX = L_citX + (Lambda[Ind3][1]/Lambda[Ind3][2])*col[i][j][ii][3]
                        L_citY = L_citY + col[i][j][ii][3]
                        L_jmen = L_jmen +  col[i][j][ii][3]/Lambda[Ind3][2]
                        Y = Y + col[i][j][ii][3]

                    if L_jmen > 0:
                        x = L_citX/L_jmen
                        y = L_citY/L_jmen                    
                 #Nyní máme smíchané celé spektrum a můžeme barvu
# převést do prostoru xyY                           
                if L_jmen > 0:
                     #Jinak se jedná o černý pixel
                    x = L_citX/L_jmen
                    y = L_citY/L_jmen
                    #Nyní známe polohu bodu v rovině xy; ještě
# musíme ověřit, že je i v prostoru sRGB
                    XX = [x,y]
                
                #Nyní se musíme podívat na polohu naší barvy vůči
# trojúhelníku, který určuje prostor sRGB

                    #Vzdálenost bodu XX od přímek, na kterých
# leží strany trojúhelníku sR sG sB
                    XXrg = Vv(sR,sG,XX)
                    XXgb = Vv(sG,sB,XX)
                    XXrb = Vv(sR,sB,XX)
                    #Vzdálenosti bodu XX od vrcholů trojúhelníku
                    RXX = V(sR,XX)
                    GXX = V(sG,XX)
                    BXX = V(sB,XX)
                    
                    sXrg = pro(sR,sG,XX) #Kolmá projekce bodu XX
# na přímku danou body sR, sG
                    sXgb = pro(sG,sB,XX) #Kolmá projekce bodu XX
# na přímku danou body sG, sB
                    sXrb = pro(sR,sB,XX) #Kolmá projekce bodu XX
# na přímku danou body sR, sB

                    #Vzdálenosti bodu XX a jeho projekcí na dané
# přímky
                    sXXrg = V(XX,sXrg)
                    sXXgb = V(XX,sXgb)
                    sXXrb = V(XX,sXrb)

                    
                    if abs((V(sR,sXrg) + V(sG,sXrg)) - RG)> 0.0000000001: #Pokud je argument roven nule, projekce bodu XX 
#leží na úsečce sRsG, tedy na straně trojúhelníku; jinak jeho 
#projekce neleží na této straně trojůhelníku a nemá smysl 
#uvožovat o této projekci

                        sXXrg = sXXrg + sXXgb + sXXrb+RXX+GXX+BXX
                    if abs((V(sG,sXgb) + V(sB,sXgb)) - GB)> 0.0000000001:#Pokud je argument roven nule, projekce bodu 
#XX leží na úsečce sBsG
                        
                        sXXgb = sXXrg + sXXgb + sXXrb+RXX+GXX+BXX
                    if abs((V(sR,sXrb) + V(sB,sXrb)) - RB) > 0.0000000001:#Pokud je argument roven nule, projekce bodu XX
# leží na úsečce sRsB
                        
                        sXXrb = sXXrg + sXXgb + sXXrb+RXX+GXX+BXX
 
                    S_tr  = 0.5 * (RG * XXrg + GB*XXgb + RB*XXrb)
#Součet ploch tří trojúhelníků sRsGXX, sGsBXX a sRsBXX;
                    
                    #Pokud je součet jejich obsahů roven obsahu
# trojůhelníku sRsGsB, pak bod XX v něm leží,
                    #a tedy i v prostoru sRGB a není třeba jej
# projektovat
                    
                    if S_tr0 == S_tr or abs(S_tr0 - S_tr) < 0.0000000001  :
                        NotOK = 0*K #Obsahy se shodují a bod XX
# leží v prostoru sRGB
                            
                    else:
                    #Bod XX leží mimo náš trojúhelník;
# na trojúhelníku se pokusíme najít bod, který je mu nejblíž
                    
                        if sXXrg == min(sXXrg, sXXgb, sXXrb, RXX, GXX, BXX):
                            x = sXrg[0] 
                            y = sXrg[1]

                            
                        elif sXXgb == min(sXXgb, sXXrb, RXX,GXX,BXX):
                            x = sXgb[0]
                            y = sXgb[1]
                            
                        elif sXXrb == min(sXXrb,RXX,GXX,BXX):
                            x = sXrb[0]
                            y = sXrb[1]
                            
                        elif RXX == min(RXX,GXX,BXX):
                            x = sR[0]
                            y = sR[1]
                            
                        elif BXX == min(GXX,BXX):
                            x = sB[0]
                            y = sB[1]
                            
                        else:
                            x = sG[0]
                            y = sG[1]

                    #Teď už můžeme převést barvu z xyY prostoru
# do XYZ a z něj do RGB        
                    XX = [x,y]
                    
                    X = Y*x / y
                    Z = Y*(1 - x - y)/y
                    
                    R = 3.2404542*X - 1.5371385*Y - 0.4985314*Z
                    G =-0.9692660*X + 1.8760108*Y + 0.0415560*Z
                    B = 0.0556434*X - 0.2040259*Y + 1.0572252*Z

                    #Normalizace
                    R = backscale(R)
                    G = backscale(G)
                    B = backscale(B)
                
                else:
                    R = 0
                    G = 0
                    B = 0
            #V průběhu výpočtu se mohly hodnoty v rámci strojové
# přesnosti dostat z intervalu [0,255], tak provedeme kontrolu
            else:
                R = 0
                G = 0
                B = 0
            if R < 0:
                R = 0
            if G < 0:
                G = 0
            if B < 0:
                B = 0

            if R > 255:
                G = int(G*255/R)
                B = int(B*255/R)
                R = 255
            if G > 255:
                R = int(R*255/G)
                B = int(B*255/G)
                G = 255
            if B > 255:
                G = int(G*255/B)
                R = int(R*255/B)
                B = 255
            
                
            images[w][i,j] = [B,G,R] #Uložení nových hodnot
# RGB pixelu

            #Do dostatečně velkých obrázků zapíšeme do spodního
# rohu rychlost, pro kterou byly dopočítány
            if height > 100 and width > 100:
                cv2.putText(images[w],str(Vx),(int(width/20),height -int(height/20)),cv2.FONT_HERSHEY_SIMPLEX, 0.5,(0,200,200),2)
                cv2.putText(images[w],str(Vx),(int(width/20),height -int(height/20)), cv2.FONT_HERSHEY_SIMPLEX, 0.5,(0,0,0),1)

    #Vytvoření pracně vytvořeného obrázku s posunutými barvami
    cv2.imwrite(name,images[w])
    #Načtení nového obrázku, který bude přepočten pro novou
# rychlost
    img = cv2.imread(Pic)
    images = images + [img]
           
print(r'konec') #To už je vše
