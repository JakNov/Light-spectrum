
import time
import numpy
import cv2
import imageio
import string

from prfc import Vv
from prfc import pro
from prfc import inter

def Spektr(R,G,B):
    #vstupem této funkce jsou RGB údaje o barvě pixelu
    #načtení hodnot pro vlnové délky z textového souboru,
    #ve které jsou uloženy informace o gamutu formát je 
    #toto:
    # vlnová délka 	x	y
    # 420	0.159581463	0.015892612

    #délka souboru v řádcích
    num_lines = sum(1 for line in open(r'Lambda.txt'))

    file = open(r'Lambda.txt')

	 #pole, ve kterém jsou údaje o okraji gamutu uloženy
    Lambda = []
    Min = 420 #maximální vlnová délka spektra v nm
    Max = 703 #minimální vlnová délka spektra v nm

    for i in range(num_lines): 
        numbers = file.readline()
        numbers = numbers.split()
        for j in range(3):
            numbers[j] = float(numbers[j])
       
        Lambda.append(numbers)
    file.close()


    #načtení druhé sady spektra, pro které budeme hledat
    #partnerky; již menší množství vlnových délek


    num_lines = sum(1 for line in open(r'Spektrum.txt'))

    file = open(r'Spektrum.txt')
	 #pole, ve kterém jsou uloženy vlnové délky, pro které
    #bude hledat partnerky
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
    Krok = 14       #vzdálenost vlnových délek spektrálních barev,
    #které procházíme s maximální a minimální vlnovou délkou

    #definování dL
    #LS počet všech spektrálních barev kterým budou vybírány
    #partnerky
    LS = len(Spektrum)

    #počet všech spektrálních barev, které známe
    #a ze kterých budou vybírány partnerky
    L = len(Lambda)

    #tady budou uloženy informace o barvách
    #[barva A][0 = intenzita barvy A z pole Spektrum,
    # 1 = pořadí barvy C v poli Lambda, 2 = intenzita barvy C,
    #3 = intenzita barvy D]; barvy C a D jsou sousední
    #spektrální barvy, jejichž smícháním dostaneme
    #partnerku k barvě A
    col = []
    col = [[0 for i in range(4)] for j in range(LS+2)]
        souc = int(R)+int(B)+int(G)
    if  souc > 0:
        R = scale(R)
        G = scale(G)
        B = scale(B)
        col[LS + 1][0] = 1
        #převod na prostor XYZ
        X = (0.4124564*R + 0.3575761*G + 0.1804375*B)
        Y = (0.2126729*R + 0.7151522*G + 0.0721750*B)
        Z = (0.0193339*R + 0.1191920*G + 0.9503041*B)

        #převod na prostor xyY
        x = X /(X + Y + Z)
        y = Y /(X + Y + Z)
        
        #definování dL;
        #proběhne LS hledání partnerek k LS spektrálních barev
        lum = Y/LS
        
        
        for i in range(LS):
            #načtení souřadnic spektrální barvy,
				#jijíž partnerku budeme hledat
            A = [Spektrum[i][1],Spektrum[i][2]]
            #souřadnice barvy, kterou chceme namíchat
				#v prostoru xyY
            XX = [x,y]
                
            for j in range(L - 2):
#najde barvu spektrální barvu A a začne od ní zkoušet,
#do kterého intervalu se trefí
                index = (i*Krok + j + 1)%L 
                C = [Lambda[index][1],Lambda[index][2]]
#dvě sousední barvy a my vyzkoušíme, jestli jejich mix
#bude partnerkou pro A
                D = [Lambda[(index +1)%L][1],Lambda[(index...

						...  +  1)%L][2]]
            
                P = inter(A,XX,C,D)#průsečík přímky AXX a CD
                
                if V(C,P) + V(D,P) - V(C,D) < 0.0000000001:
#pokud je průsečík na přímce CD, pak lze partnerku namíchat
#z barev CD a tedy je započteme
                    
                    K = A[1]*(P[0] - XX[0])/(P[1]*(XX[0] - A[0]))
#vztah 2.5
                     
                    col[i][0] = lum * K / (K + 1)
#vztah 2.5; intenzita první barvy A
                    col[i][1] = index
#pozice barvy C; pozice D je automaticky index +1
                    K = C[1]*(D[0] - P[0])/(D[1]*(P[0] - C[0]))
#vztah 2.5 pro mix C a D

                    col[i][2] = (lum - col[i][0])* K / (K + 1)
#intenzita barvy C
                    col[i][3] = lum - col[i][0] - col[i][2]
#intenzita barvy D

                    
                    j = L - 2 #už není potřeba hledat
   
    Spc = []
    Spc = [[0 for i in range(2)] for j in range(L)]
    for i in range(L):
        Spc[i][0] = Lambda[i][0]
   
    for i in range(LS):
        Spc[int(i*Krok)][1] = Spc[int(i*Krok)][1] + col[i][0]
        Spc[col[i][1]][1] = Spc[col[i][1]][1] + col[i][2]
        Spc[(col[i][1]+1)%L][1] = Spc[(col[i][1]+1)%L][1]...
... + col[i][3]

    return(Spc) #vrátí pole vlnových délek s intenzitou,
#kterou jsou ve spektru dané barvy se souřadnicemi RGB zastoupeny
