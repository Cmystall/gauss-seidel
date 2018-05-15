from typing import List, Any


def altura(m):
    return len(m)


def largura(m):
    return len(m[0])


def issparse(m):
    a = altura(m)
    l = largura(m)
    maxnumelem = a * l
    contador = 0
    for i in (range(a)):
        for j in range(l):
            if i == j and m[i][j] == 0.0:
                return False
            if m[i][j] != 0.0:
                contador += 1
    if contador > int(maxnumelem / 2):
        return False
    else:
        return True


def csr(m):
    csr = [[], [], []]  # type: List[List[float]]
    a = altura(m)
    l = largura(m)
    if a == l:
        for i in range(a):
            if m[i][i] != 0.0:
                csr[0].append(m[i][i])
                csr[1].append(i)
                csr[2].append(len(csr[0])-1)
            for j in range(l):
                if j != i and m[i][j] != 0.0:
                    csr[0].append(m[i][j])
                    csr[1].append(j)

    else:
        for i in range(a):
            for j in range(l):
                if m[i][j] != 0.0:
                    csr[0].append(m[i][j])
                    csr[1].append(i)
                    csr[2].append(len(csr[0]))
    return csr

def criarmatriz(m):
    entrada = open(m, 'r')
    novo = []
    for linha in entrada.readlines():
        linha = [float(x) for x in linha.split(" ")]
        novo.append(linha)
    entrada.close()
    print(novo)
    return novo
