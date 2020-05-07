def generatepq(P,Q):
    """
    (Algorithm 3.12: Lenstra, Verheul: An overview of the XTR public key system)

    Generates the parameters p and q of XTR.

    Parameters:

    P (int): bits of p;

    Q (int): bits of q

    Returns:

    int: p,q, primes such that q|p^2-p+1

    """

    q=prime_range(2**Q,2**(Q+1))
    for i in q:
        if i.mod(12)==7:
            q=i
            break
    r=(1-sqrt(GF(q)(-3)))*2.inverse_mod(q)
    r=r.mod(q)

    p=prime_range(2**(P),2**(P+1))
    for i in p:
        if i.mod(3)==2 and (i-r).mod(q)==0:
            p=i
            break
    ##Si no funciona se puede intentar implementar con el otro valor de r.
    return(p,q)

def irreducibility(c):
    """
    (Algorithm 3.33: Lenstra, Verheul: An overview of the XTR public key system)

    Checks if F(c,X)=X^3-cX^2+c^pX-1 is irreducible over GF(p^2). If c chosen randomly,
    then P(F(c,X) irreducible)=1/3.

    Parameters:

    c (GF(p^2))

    Returns:
    True if F(c,X) is irreducible
    False otherwise

    """
    if c not in F2:
        return(irreducibility(F2(c)))
    f0=(-27+9*c**(p+1)-2*c**3)/27
    f1=c**p-c**2/3
    disc=f0**2+4*(f1/3)**3
    if kronecker(ZZ(disc),p)==-1:
        return(false)
    r1=(-f0+sqrt(F(disc)))/2
    y=r1**((p+1)/3)
    if y!=y**p:
        return(true)
    else:
        return(false)

def computetrace(p,q):
    """
    (Algorithm 3.34: Lenstra, Verheul: An overview of the XTR public key system)

    Computes Trace(g), where g is a generator of the XTR subgroup.

    Parameters:

    (int) p,q: with the conditions of XTR
    
    Returns:
    (GF(p^2)) Trace(g)

    """
    c=F2.random_element()
    while irreducibility(c)==false or c in F:
        c=F2.random_element()
    n=(p**2-p+1)/q
    d=computeS(n,c)[1]
    if d==3:
        computetrace()
    return(d)
def computeS(n,c):
    """
    (Algorithm 2.35: Lenstra, Verheul: An overview of the XTR public key system)

    Computes S_n(c) as defined in Lenstra, Verheul.

    Parameters:

    (int) n>0;

    (GF(p^2)) c

    Returns:
    (tup) S_n(c)=(c_{n-1},c_n,c_{n+1})
    """
    if type(c) != type(F2(1)):
        return(computeS(n,F2(c)))
    if n==0:
        return(c**p,F2(3),c)
    elif n<0:
        (a1,a2,a3)=computeS(-n,c)
        return(a3**p,a2**p,a1**p)
    elif n==1:
        return(F2(3),c,c**2-2*c**p)
    elif n==2:
        (c0,c1,c2)=computeS(1,c)
        return(c1,c2,c*c2-c**p*c1+c0)
    elif n==3:
        (c1,c2,c3)=computeS(2,c)
        return(c2,c3,c*c3-c**p*c2+c1)
    else:
        if n%2==1:
            m=n
        else:
            m=n-1
        k=1
        Sk=computeS(3,c) #Sk=(c_{2k},c_{2k+1},c_{2k+2})

        binario=bin((m-1)/2)[2::]
        r=len(binario)
        contador=0
        if r!= 1:
            for j in range(r):
                if j==0:
                    continue
                if eval(binario[j])==0:
                    Sk=(Sk[0]**2-2*Sk[0]**p,
                         Sk[0]*Sk[1]-c**p*Sk[1]**p+Sk[2]**p,
                         Sk[1]**2-2*Sk[1]**p) #S2k=(c_{4k},c_{4k+1},c_{4k+2})
                    k=2*k
                else:
                    Sk=(Sk[1]**2-2*Sk[1]**p,
                         Sk[2]*Sk[1]-c*Sk[1]**p+Sk[0]**p,
                         Sk[2]**2-2*Sk[2]**p) #S2k+1=(c_{4k+2},c_{4k+3},c_{4k+4})
                    k=2*k+1
        if m==n:
            return(Sk)
        else:
            return(Sk[1],Sk[2],c*Sk[2]-c**p*Sk[1]+Sk[0])
