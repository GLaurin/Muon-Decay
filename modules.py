#%%

import numpy as np

#%%

def OptimiseFunc(Func, x, data, d_err, n, conv, N=10, pm=0.1, pc=0.8, sigm=0.2, cc_max=500, cr_co=1, prt=True):
    """
    Optimisation d'une fonction sur des donnees par minimisation du Chi carre
    Func  : fonction a optimiser;  Func(x,*args)
    x     : domaine discret de data
    data  : valeurs sur lesquelles ajuster la fonction
    d_err : erreurs sur les valeurs; len(d_err)==1 ou d_err.size==data.size
    n     : nombre de parametre a optimiser; n==len(args)
    conv  : facteurs de conversation des solutions normalisees 
    N     : nombre de vecteurs tests
    pm    : probabilite de mutation
    pc    : probabilite de croisement
    sigm  : largeur gaussienne des mutations (normalise)
    cc_max: nombre maximal d'iteration
    cr_co : critere de convergence
    """

    def Roulette(rang,N):
        """
        Selection aleatoire mais favoritiste d'un indice parent
        """
        S = N*(N+1)/2                   # Valeur maximale de la boucle
        r = np.random.uniform()*S       # Critere aleatoire de selection
        s_cum = 0                       # Somme cumulative des rangs
        for k in range(N):              # Boucle sur le nombre de solutions
            s_cum += (1+rang[k])        # MaJ de la somme cumulative
            if s_cum >= r and k >= 2:   # Test de selection
                return k

    def Reproduction(u1,u2,pm,pc,sigm=0.2):
        """
        Generation de deux nouveaux vecteurs solutions a partir de
        deux vecteurs-parents de parametres normalises. 
        """
        
        ## Initialisation des nouveaux vecteurs 
        nu1 = np.zeros_like(u1)
        nu2 = np.zeros_like(u1)
        
        ## Croisement
        rc = np.random.uniform(size=u1.size)
        ic = (rc > pc)                          # Test probabiliste
        gc = np.random.uniform(size=u1.size)    # Facteur de croisement
        gc[ic == False] = 0
        
        ## Mutation
        rm = np.random.uniform(size=u1.size)
        im = (rm > pm)                          # Test probabiliste
        g1 = np.random.normal(0,sigm,u1.size)   # Facteurs de mutation
        g2 = np.random.normal(0,sigm,u1.size)
        g1[im] = 0
        g2[im] = 0
        
        ## Détermination des nouveaux vecteurs
        nu1 = abs((   gc *u1 + (1-gc)*u2) + g1)
        nu2 = abs(((1-gc)*u1 +    gc *u2) + g2)
        
        ## Application des bornes
        nu1[nu1>=1] = (nu1-nu1//1)[nu1>=1]
        nu2[nu2>=1] = (nu2-nu2//1)[nu2>=1]
        nu1[nu1<=0] = (nu1-nu1//1)[nu1<=0]
        nu2[nu2<=0] = (nu2-nu2//1)[nu2<=0]
        
        return nu1, nu2

    def Evaluation(u, conv, data, d_err, x, N, Func):
        """
        Evalutation du Chi carre pour chaque solution
        """
        ev = np.zeros(N)
        for i in range(N):
            args  = u[:,i]*conv
            y     = Func(x, *args)
            ev[i] = 1/np.sum((y-data)**2/d_err**2)
        return ev

    ### MAIN
    u    = np.random.uniform(0,1,size=(n,N))              # Solutions aleatoires
    f    = Evaluation(u, conv, data, d_err, x, N, Func)   # Evaluation des solutions
    rang = np.argsort(f)                                  # Classement initial

    ## Resolution
    cc = 0                                          # Compteur d'iterations
    while (f.max()<cr_co) and (cc<cc_max):
  
        ## Selection des solutions parents
        k1 = Roulette(rang,N)                       # Premier parent
        k2 = Roulette(rang,N)                       # Deuxieme parent
        while k2==k1:                               # Verification de parents differents
            k2 = Roulette(rang,N)                   # MaJ si necessaire

        ## Reproduction
        nu = np.zeros((n,N))                        # Initialisation des nouvelles solutions
        for i in range(N//2):                       # Calculs des nouvelles solutions
            nu[:,i*2], nu[:,i*2+1] = Reproduction(u[:,k1], u[:,k2], pm, pc, sigm)

        ## Remplacement des solutions
        nf             = Evaluation(nu, conv, data, d_err, x, N, Func)  # Nouvelles evaluations
        nrang          = np.argsort(nf)             # Nouveau rang
        u[:,rang[:-1]] = nu[:,nrang[1:]]            # MaJ elitiste des solutions
        f[rang[:-1]]   = nf[nrang[1:]]              # MaJ elitiste des evaluations
        rang           = np.argsort(f)              # MaJ elitiste des rangs

        if cc%(cc_max//10)==0 and prt==True:        # Affichage de la progression
            print(f"{cc/cc_max*100:.1f}% - {f.max():.4f}")
        cc +=1                                      # MaJ du compteur
        
    return u[:,f.argmax()]*conv, f.max()
