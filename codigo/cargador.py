import numpy as np
import random


def test_bc_3_1():
    file = open("test/input_wc_10.txt","w")
    cantElemMax = 28
    cantElem = 2
    i = cantElem
    while i <= cantElemMax:
        valorObjetivo = np.random.randint(i, i**2)
        file.write(str(i) + " " + str(valorObjetivo) + "\n")
        #file.write(str(valorObjetivo) + "\n")
        for j in range(0,i):
            if(i==cantElemMax and (j==i-1)):
                file.write(str(np.random.randint(0,i)))
            else:
                file.write(str(np.random.randint(0,i)) + "\n")
        i+=5
    file.close()


def test_wc_2():
    file=open("test/input_wc_2_2.txt","w")
    cantElem = 26
    valorObjetivo = cantElem
    file.write(str(cantElem) + " " + str(valorObjetivo) + "\n")
    for j in range(0,2):
        for k in range(0,cantElem):
            if(j==1 and k==cantElem-1):
                file.write(str(valorObjetivo))
            else:
                if(j==0 and 0==k):
                    file.write(str(0) + "\n")
                else:
                    file.write(str(np.random.randint(0,cantElem)) + "\n")

def test_wc_3():
    file = open("test/input_wc_14.txt","w")
    n = 20
    while n<=20:
        v=n
        while v<=n**3:
            file.write(str(n) + " " + str(v) + "\n")
            valor = np.random.normal(loc=v,scale=v/n**3,size=n)
            for i in range(0,n):
                file.write(str(int(abs(valor[i])))+ "\n")
            v+=100
        n+=5

def test_wc_2_2():
    file = open("test/input_wc_2_2.txt","w")
    n = 50
    vo = n
    while vo <= n**2:
        file.write(str(n) + " " + str(vo) + "\n")
        for j in range(0, n):
            file.write(str(np.random.randint(n,n**2)) + "\n")
        vo+=n

    file.close()

def main():
    test_wc_3()
main()
