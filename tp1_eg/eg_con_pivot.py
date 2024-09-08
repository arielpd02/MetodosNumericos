import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from numpy import *
from numpy.linalg import *


eps=sys.float_info.epsilon

#--------------EG con pivot---------------

def pivot(a:np.array,i)->int:
    max=abs(a[i,i])
    res=0
    for k in range(i+1,len(a)):
        if(abs(a[k,i])>max):
            max=a[k,i]
            res=k
    return res

def producto_punto(v:np.array,w:np.array):
    res=0
    for i in range(0,v.size):
        res+=v[i]*w[i]
    return res

#Para poder limitar el error numérico y que el algoritmo sea correcto (i.e. devuelva S.I. cuando no hay solucion o alguna si existe) , consideramos que un x<=e-6 es efectivamente 0.Es decir , nos basamos en una precision de 6 digitos en los resultados.

def eliminacion_gaussiana_pivot(m:np.array,b:np.array)->np.array:
    if(len(m)!=len(m[0])):
        print("Error. Matriz no cuadrada")

    cota_decimal=eps*pow(10,6)

    #Construimos A' , extension de A con b
    n=len(m)
    m_ext=np.hstack((m,np.zeros((n,1))))
    m_ext[:,n]=b

    for i in range(0,n):

        fila_a_permutar=pivot(m_ext,i)
        if(fila_a_permutar>i):      #Para el caso que devuelve 0
            copy=np.array(m_ext[i,:])
            m_ext[i,:]=m_ext[fila_a_permutar,:]
            m_ext[fila_a_permutar,:]=copy

        #Si la columna ya es nula,paso.
        if(np.all(m_ext[i+1:,i] == 0)):
            continue
        
        for j in range(i+1,n):
            if(np.abs(m_ext[i,i])<cota_decimal):
                print("Posible Error Numérico.Cociente con divisor cercano a 0.\n")
            multiplier=m_ext[j,i]/m_ext[i,i]
            m_ext[j,:]=m_ext[j,:]-m_ext[i,:]*multiplier     #M_j=M_j - M_i * multiplier
            
    #print(m_ext)

    #Excepcion a devolver
    for i in range(0,n):
        if(np.abs(m_ext[i,i])>cota_decimal): #M_ii>0        
            continue
        elif(np.abs(m_ext[i,n])>cota_decimal):#    b_ii>0                    
            print("Sistema incompatible.No existen soluciones.")
            return 
        else:               #TODO:que hago si es SCI
            print("Hay infinitas soluciones.")            
            return


    #Ya A' escalonada -> despejamos x
    x=np.zeros(n)
    for i in range(n-1,-1,-1):
        b_i=m_ext[i,n]
        if(i==n-1):
            x[i]=b_i/m_ext[i,i]
            continue
        prod=producto_punto(m_ext[i,i+1:n],x[i+1:n])
        x[i]=(b_i-prod)/m_ext[i,i]

    #Devolvemos X=[x1,...,xn]
    return x


