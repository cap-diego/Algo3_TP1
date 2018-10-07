import numpy as np
import random


def test_peorcaso_var_n():
    file = open("test/test_peor_caso_variando_n_v1000.txt","w")
    cantElemMax = 140
    cantElem = 0
    i = cantElem
    while i <= cantElemMax:
        valorObjetivo = 1000
        #valorObjetivo = np.random.randint(0,1+i**2)
        file.write(str(i) + " " + str(valorObjetivo) + "\n")
        #file.write(str(valorObjetivo) + "\n")
        for j in range(0,i):
            val = np.random.randint(0,valorObjetivo)
            if((j==cantElemMax-1)):
                #file.write(str(np.random.randint(0,i)))
                file.write(str(val))
            else:
                file.write(str(val)+"\n")
                #file.write(str(np.random.randint(0,i)) + "\n")
        i+=2
    file.close()

def test_peorcaso_var_ambos():
    file = open("test/test_peor_caso_variando_ambos.txt","w")
    cantElemMax = 60
    cantElem = 5
    i = cantElem
    while i <= cantElemMax:
        vo = 0
        while vo<=50000:
            file.write(str(i) + " " + str(vo) + "\n")
            #file.write(str(valorObjetivo) + "\n")
            for j in range(0,i):
                val = np.random.randint(0,i)
                if(vo==50000 and (j==i-1) and i==cantElemMax):
                    file.write(str(val))
                    #file.write(str(np.random.randint(0,i**2)))
                else:
                    file.write(str(val) + "\n")
                    #file.write(str(np.random.randint(0,i**2)) + "\n")
            vo+=1000
        i+=3

    file.close()


def test_peorcaso_var_v():
        file = open("test/test_peor_caso_variando_v_n15.txt","w")
        cantElem = 25
        vo = 100
        while vo <= 50000:
            file.write(str(cantElem) + " " + str(vo) +"\n")
            val=np.random.randint(0,cantElem)
            for j in range(0,cantElem):
                if(vo==50000 and (j==cantElem-1)):
                    file.write(str(val))
                    #file.write(str(np.random.randint(0,cantElem)))
                else:
                    file.write(str(val) + "\n")
                    #file.write(str(np.random.randint(0,cantElem)) + "\n")
            vo+=1000
        file.close()

def main():
    print("Generado correctamente")
    test_peorcaso_var_ambos()
main()
