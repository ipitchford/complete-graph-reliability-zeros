"""Plot the 78 nontrivial zeros of K14; check one disjoint disk per zero."""
from pathlib import Path
from math import comb
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt

from verify import reliability_polynomials
from certified_intervals import iv


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)


def main():
    n = 14
    r = n-1
    N = n*r//2
    d = N-r
    rel = reliability_polynomials(n)[n]
    h = [sum(rel[j]*comb(r+k-j-1, r-1) for j in range(k+1)) for k in range(d+1)]
    require(all(x>0 for x in h))
    upper = sum(k*(k-1)*h[k] for k in range(2,d+1))
    mp.mp.dps = 60
    eps = iv.mpf('1e-22')
    centres, visible = [], []
    for initial in np.roots(h[::-1]):
        root = mp.mpc(initial)
        for _ in range(35):
            value, derivative = mp.polyval(h[::-1],root,derivative=True)
            delta = value/derivative
            root -= delta
            if abs(delta)<mp.mpf('1e-50'):
                break
        z = iv.mpc(mp.nstr(root.real,53),mp.nstr(root.imag,53))
        require(abs(z)+eps < 1)
        f, fp = iv.mpc(0), iv.mpc(0)
        for k in range(d,-1,-1):
            fp = fp*z+f
            f = f*z+h[k]
        require(abs(f)+upper*eps**2/2 < abs(fp)*eps)
        require(all(abs(z-c)>2*eps for c in centres))
        centres.append(z)
        visible.append(complex(root))
    require(len(centres)==d)
    print('K14:',d,'disjoint certified root disks, each of radius 1e-22; 13 further roots at 1.')
    zeros = np.array(visible)
    fig, (first,second) = plt.subplots(1,2,figsize=(8.6,3.05),gridspec_kw={'width_ratios':[1,1.36]})
    ink,teal,berry,muted = '#272532','#00796b','#99375c','#77717c'
    theta = np.linspace(0,2*np.pi,400)
    first.plot(np.cos(theta),np.sin(theta),c=muted,lw=.9)
    first.scatter(zeros.real,zeros.imag,s=8,color=teal)
    first.scatter([1],[0],s=20,c=berry)
    first.annotate('multiplicity 13',(1,0),(.01,.14),fontsize=7,color=berry,arrowprops={'arrowstyle':'-','lw':.7,'color':berry})
    first.set(xlabel='Real part of q',ylabel='Imaginary part of q',xlim=(-1.14,1.14),ylim=(-1.14,1.14),aspect='equal',xticks=[-1,0,1],yticks=[-1,0,1])
    upperroots = zeros[zeros.imag>0]
    second.scatter(np.angle(upperroots),abs(upperroots),s=13,color=teal)
    second.axvline(2.393242274134895,c=berry,ls='--',lw=1)
    second.text(2.35,.839,'limiting direction',rotation=90,va='top',ha='right',fontsize=8,color=berry)
    second.set(xlim=(0,np.pi),ylim=(.715,.855),xlabel='Argument of q in the upper half-plane (radians)',ylabel='Modulus of q',xticks=[0,1,2,3])
    for ax in (first,second):
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=8,colors=ink)
        ax.xaxis.label.set_size(9); ax.yaxis.label.set_size(9)
    fig.tight_layout(pad=.55,w_pad=1.3)
    out = Path(__file__).resolve().parent.parent / 'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out / 'clique-roots.pdf', bbox_inches='tight')
    fig.savefig(out / 'clique-roots.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':
    main()

