# -*- coding: iso-8859-15 -*-


def TracePR(mP,mR) :
    plt.ioff()
    Pmax_level=mP.data.max()
    Rmax_level=mR.data.max()
    # ==========================================
    plt.figure()
    plt.contourf(X,Y,mP.todense(),np.linspace(0.0001,Pmax_level,101))
    V1=plt.axis()
    plt.plot(vec_teta,Th*(vec_teta-1.)/vec_teta,'r-')
    plt.plot(vec_teta,Tc*(vec_teta-1.),'r-')
    plt.axis(V1)
    plt.ylabel(r'$\rm{\Delta T\ \left[K\right]}$')
    plt.xlabel(r'$\theta= {T_M}/{T_m}\left[\mathrm{-}\right]$')
    plt.colorbar(ticks=np.arange(0,Pmax_level+1,5))
    plt.title(r'$\rm{Puissance\ \dot{W}\ \left[kW\right]}$')
    # ==========================================
    plt.figure()
    plt.contourf(X,Y,mR.todense(),np.linspace(0.0001,Rmax_level,101))
    V=plt.axis()
    plt.plot(vec_teta,Th*(vec_teta-1.)/vec_teta,'r-')
    plt.plot(vec_teta,Tc*(vec_teta-1.),'r-')
    plt.axis(V)
    plt.ylabel(r'$\rm{\Delta T\ \left[K\right]}$')
    plt.xlabel(r'$\theta= {T_M}/{T_m}\left[\mathrm{-}\right]$')
    plt.colorbar(ticks=np.arange(0,Rmax_level+0.001,0.01))
    plt.title(r'$\rm{Rendement\ \mu}$')
    # ==========================================
    plt.figure()
    col=(mR.data**2+mP.data**2)**0.5
    plt.scatter(mR.data,mP.data,s=20,c=col,marker='o',edgecolors="none")
    plt.axis([0,Rmax_level+0.01,0,Pmax_level+2])
    plt.colorbar()
    plt.xlabel(r'$\rm{Rendement\ \left(\mu\right)\   \left[-\right] }$')
    plt.ylabel(r'$\rm{Puissance\ \left(\dot{W}\right)\ \left[kW\right]}$')
    plt.title(r'$\rm{Courbes\ Puissance\ -\ Rendement}$')
    Formule=r'$\rm{Couleur=}\ \sqrt{P^2+\mu^2}$'
    plt.text(0.45,0.1,Formule,transform = plt.gca().transAxes)
    plt.ion()
    plt.show()


def Trace_Optim(fignum,M) :
    plt.figure(fignum)
    x=M[:,colx]
    y=M[:,coly]
    del M
    plt.subplot(3,2,1)
    plt.scatter(x,y,c=M[:,0],s=10,edgecolors="none")
    plt.colorbar()
    plt.text(0.05,120,r'$\alpha$',fontsize=24)

    plt.subplot(3,2,2)
    plt.scatter(x,y,c=M[:,1],s=10,edgecolors="none")
    plt.colorbar()
    plt.text(0.05,120,r'$K_h$',fontsize=24)

    plt.subplot(3,2,3)
    plt.scatter(x,y,c=M[:,2],s=10,edgecolors="none")
    plt.colorbar()
    plt.text(0.05,120,r'$K_c$',fontsize=24)

    plt.subplot(3,2,4)
    plt.scatter(x,y,c=M[:,3],s=10,edgecolors="none")
    plt.colorbar()
    plt.text(0.05,120,r'$K_l$',fontsize=24)

    plt.subplot(3,2,5)
    plt.scatter(x,y,c=M[:,4],s=10,edgecolors="none")
    plt.colorbar()
    plt.text(0.05,120,r'$K_{reg}$',fontsize=24)
