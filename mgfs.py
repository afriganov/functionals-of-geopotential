from numpy import *
elipsoid = 1
if elipsoid == 1: #MARSOV
    aell =  3396000
    fel =  1/196.877360
    w = 7.088218066303858E-05
    el = math.sqrt(fel*(2-fel))
    bell = aell*math.sqrt(1-math.pow(el,2))
    GMell = 0.4282837285418775E+14
    m = (math.pow(w,2)*math.pow(aell,2)*bell)/GMell
    elc = math.sqrt(math.pow(aell,2)-math.pow(bell,2))/bell
    q0 = 0.5*(((1+(3/math.pow(elc,2)))*math.atan(elc))-(3/elc))
    J2 = (math.pow(el,2)/3)*(1-((2/15)*(m*elc/q0)))
    Cp20 = -J2/math.sqrt(5)
def Pnm():
    if n == 0 and m == 0:
        return 1
    elif n == 1 and m == 1:
        return math.sqrt(3)*u
    elif n == m and m > 1:
        return u*math.sqrt((2*m+1)/(2*m))*tablepnm[m-1,m-1]
    else:
        anm = math.sqrt(((2*n-1)*(2*n+1))/((n-m)*(n+m)))
        bnm = math.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((n-m)*(n+m)*(2*n-3)))
        return anm*t*tablepnm[n-1,m]-bnm*tablepnm[n-2,m]

#INPUT
nmax = int(input("Do kojeg reda zelite razviti ?"))
fie1 = int(input("Unesite pocetnu sirinu:"))
lae1 = int(input("Unesite pocetnu duzinu:"))
fie2 = int(input("Unesite zavrsnu sirinu:"))
lae2 = int(input("Unesite zavrsnu duzinu:"))
laem = lae1
fiem = fie1
k = 0
l = 0
R = 3396000
GM = GMell
tablerje = np.zeros((abs(fie1)+fie2+1,lae2+1))
Cell = np.zeros((nmax + 1, nmax + 1))
Sell = np.zeros((nmax + 1, nmax + 1))
dC = np.zeros((nmax + 1, nmax + 1))
dS = np.zeros((nmax + 1, nmax + 1))
tablepnm = np.zeros((nmax + 1, nmax + 1))
tablerkr = np.zeros(nmax + 1)
Cm = np.zeros((nmax + 1, nmax + 1))
Sm = np.zeros((nmax + 1, nmax + 1))
with open("ggm3.txt") as f:
    C = np.zeros((nmax + 1, nmax + 1))
    S = np.zeros((nmax + 1, nmax + 1))
    for lines in f:
        lista = lines.split()
        n = int(lista[0])
        m = int(lista[1])
        C[n, m] = lista[2]
        S[n, m] = lista[3]
        if n == nmax and n == m:
            break

for i in range(0, 11):
    Cell[2 * i, 0] = math.pow(-1, i) * ((3 * math.pow(el, (2 * i))) / ((2 * i + 1) * (2 * i + 3) * math.sqrt(4 * i + 1))) * (
                                     1 - i - (math.pow(5, (3 / 2)) * i * (Cp20 / math.pow(el, 2))))
for n in range(0, nmax + 1):
    for m in range(0, n + 1):
        dC[n, m] = C[n, m] - Cell[n, m] * (GMell / GM) * math.pow((aell / R), n)
        dS[n, m] = S[n, m]
        m += 1
    n += 1
for fie1 in range(fiem,fie2+1):
    fie = math.radians(fie1)
    for lae1 in range(laem,lae2+1):
        lae = math.radians(lae1)
        # Varijable i izracun
        V = 0
        h = 0
        la = lae
        N = aell / (math.sqrt(1 - math.pow(el, 2) * math.pow(math.sin(fie), 2)))
        x = (N + h) * math.cos(fie) * math.cos(lae)
        y = (N + h) * math.cos(fie) * math.sin(lae)
        z = ((1 - math.pow(el, 2)) * N + h) * sin(fie)
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))
        fi = math.atan(z / (math.sqrt(math.pow(x, 2) + math.pow(y, 2))))
        t = math.sin(fi)
        u = math.cos(fi)




        for n in range(0, nmax + 1):
            for m in range(0, n + 1):
                tablepnm[n, m] = Pnm()
                m += 1
            n += 1
        for n in range(0, nmax + 1):
            rkr = pow(R / r, n)
            tablerkr[n] = rkr

        for n in range(0, nmax + 1):
            for m in range(0, n + 1):
                V = V + tablerkr[n] * (n - 1) * (
                        dC[n, m] * math.cos(math.radians(m * math.degrees(la))) + dS[n, m] * math.sin(
                    math.radians(m * math.degrees(la)))) * tablepnm[n, m]
                m += 1
            n += 1
        dg = (GM / math.pow(r, 2)) * V * 10E+4
        tablerje[k,l] = dg
        l = l + 1
        lae1 = lae1 +1
    l = 0
    fie1 = fie1 + 1
    k = k + 1
print(tablerje)
with open("finalni120.txt", "w") as f:
    for i in range(0,abs(fiem)+fie2+1):
        j = 0
        for j in range(0,laem+lae2+1):
            f.write(str(tablerje[i][j]))
            f.write(" ")
        f.write("\n")
