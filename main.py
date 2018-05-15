import matriz as mtz
import time

def gauss_seidel(m, b, x0=None, eps=1e-4, maxiter=100):

    a = mtz.altura(m)
    # l = mtz.largura()

    x0 = [0]*a if x0 is None else x0
    x1 = x0[:]
    if not mtz.issparse(m):
        print("É esparsa")
        csr = mtz.csr(m)
        aa = csr[0]
        ia = csr[1]
        ja = csr[2]
        print(aa)
        print(ia)
        print(ja)
        for _ in range(maxiter):
            for indi in range(len(ja)):
                if indi == len(ja) - 1:
                    indf = len(aa)
                    aux = aa[ja[indi]: indf]
                    aux2 = ia[ja[indi]: indf]
                else:
                    indf = ja[indi+1]
                    aux = aa[ja[indi]: indf]
                    aux2 = ia[ja[indi]: indf]
                soma = 0
                for i in range(len(aux)):
                    soma -= aux[i]*x1[aux2[i]] if indi != aux2[i] else 0
                x1[indi] = (soma + float(b[indi][0])) / aux[0]
            if all(abs(x1[i]-x0[i]) < eps for i in range(a)):
                res = open("resposta.txt", "w")
                for i in range(len(x1)):
                    res.write(str(x1[i]))
                    res.write("\n")
                res.close()
                return x1
            x0 = x1[:]
        raise ValueError("A solucao nao converge")
    else:
        print ("Não é esparsa")
        for _ in range(maxiter):
            for i in range(a):
                soma = sum(-m[i][j]*x1[j] for j in range(a) if i != j)
                x1[i] = (b[i][0] + soma) / m[i][i]
            if all(abs(x1[i]-x0[i]) < eps for i in range(a)):
                return x1
            x0 = x1[:]


if __name__ == '__main__':
    coeficientes = input("Insira o nome do primeiro arquivo (coeficientes da Matriz):")
    valores = input("Insira o nome do segundo arquivo (valores b):")
    mat = mtz.criarmatriz(coeficientes)
    val = mtz.criarmatriz(valores)
    start_time = time.time()
    print(gauss_seidel(mat, val))
    print("--- %s segundos de execução ---" % (time.time() - start_time))
