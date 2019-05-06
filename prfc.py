import time
import numpy
import cv2
import imageio
import string

def V(A,B):
    v = numpy.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2 )
    return v

def VV(A):
    v = numpy.sqrt((A[0])**2 + (A[1])**2 + (A[2])**2)
    return v

def Vv(A,B,C):
    p = [A[0] - B[0],A[1] - B[1]]
    v = abs(-p[1]*C[0] + p[0]*C[1] - (-p[1]*A[0] + p[0]*A[1]))/numpy.sqrt(p[0]**2 + p[1]**2)
    return v

def pro(A,B,C):
    p = [-A[1] + B[1],A[0] - B[0]]
    c = -(p[0]*A[0] + p[1]*A[1])
    x1 = (p[1]**2 * C[0] - p[0]*p[1]*C[1] - p[0]*c)/(p[0]**2 + p[1]**2)
    x2 = (p[0]**2 * C[1] - p[0]*p[1]*C[0] - p[1]*c)/(p[0]**2 + p[1]**2)
    X = [x1,x2]
    return X
    
def inter(A,B,C,D):
    p = [-A[1] + B[1],A[0] - B[0]]
    c = -(p[0]*A[0] + p[1]*A[1])

    q = [-C[1] + D[1],C[0] - D[0]]
    cc = -(q[0]*C[0] + q[1]*C[1])
    
    x1 = (cc * p[1] - c*q[1])/(q[1] * p[0] - p[1]*q[0])
    x2 = (c * q[0] - cc * p[0])/(q[1] * p[0] - p[1]*q[0])
    X = [x1,x2]
    return X

def scale(V):
    if (V/255) <= 0.04045:
        return ((V/255)/12.92)
    else:
        return (((V/255) + 0.055)/1.055)**2.4

def backscale(V):
    if (V) <= 0.0031308:
        return ((V*255)*12.92)
    else:
        return ((V**(1/2.4))*1.055 - 0.055)*255   
        
