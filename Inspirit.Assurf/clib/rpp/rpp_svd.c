#include <math.h>

#include "inc/rpp_types.h"
#include "inc/rpp_vecmat.h"

typedef double SVD_FLOAT;


static SVD_FLOAT at,bt,ct;
static SVD_FLOAT maxarg1,maxarg2;

#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
    (ct=bt/at,at*_sqrt((SVD_FLOAT)(1.0)+ct*ct)) : (bt ? (ct=at/bt,bt*_sqrt((SVD_FLOAT)(1.0)+ct*ct)): (SVD_FLOAT)(0.0)))

#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= (SVD_FLOAT)(0.0) ? fabs(a) : -fabs(a))

int svdcmp( double *a, int m,int n, double *w,double *v)
{
    int flag,i,its,j,jj,k,ii=0,nm=0;
    SVD_FLOAT c,f,h,s,x,y,z;
    SVD_FLOAT anorm=0.0,g=0.0,scale=0.0;

    if (m < n) return -1;	// must augment A with extra zero rows

	//assert(n==3);
	SVD_FLOAT rv1[3];		// was: rv1=G_alloc_vector(n);

    n--;
    m--;

    for (i=0;i<=n;i++) {
        ii=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += (SVD_FLOAT)fabs(a[k*3+i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k*3+i] /= (SVD_FLOAT)scale;
                    s += a[k*3+i]*a[k*3+i];
                }
                f=a[i*3+i];
                g = -SIGN(_sqrt(s),f);
                h=f*g-s;
                a[i*3+i]=(f-g);
                if (i != n) {
                    for (j=ii;j<=n;j++) {
                        for (s=0.0,k=i;k<=m;k++) s += a[k*3+i]*a[k*3+j];
                        f=s/h;
                        for (k=i;k<=m;k++) a[k*3+j] += (SVD_FLOAT)(f)*a[k*3+i];
                    }
                }
                for (k=i;k<=m;k++) a[k*3+i] *= (SVD_FLOAT)(scale);
            }
        }
        w[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m && i != n) {
            for (k=ii;k<=n;k++) scale += fabs(a[i*3+k]);
            if (scale) {
                for (k=ii;k<=n;k++) {
                    a[i*3+k] /= scale;
                    s += a[i*3+k]*a[i*3+k];
                }
                f=a[i*3+ii];
                g = -SIGN(_sqrt(s),f);
                h=f*g-s;
                a[i*3+ii]=f-g;
                for (k=ii;k<=n;k++) rv1[k]=a[i*3+k]/h;
                if (i != m) {
                    for (j=ii;j<=m;j++) {
                        for (s=0.0,k=ii;k<=n;k++) s += a[j*3+k]*a[i*3+k];
                        for (k=ii;k<=n;k++) a[j*3+k] += s*rv1[k];
                    }
                }
                for (k=ii;k<=n;k++) a[i*3+k] *= scale;
            }
        }
        anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=0;i--) {
        if (i < n) {
            if (g) {
                for (j=ii;j<=n;j++)
                    v[j*3+i]=(a[i*3+j]/a[i*3+ii])/g;
                for (j=ii;j<=n;j++) {
                    for (s=0.0,k=ii;k<=n;k++) s += a[i*3+k]*v[k*3+j];
                    for (k=ii;k<=n;k++) v[k*3+j] += s*v[k*3+i];
                }
            }
            for (j=ii;j<=n;j++) v[i*3+j]=v[j*3+i]=0.0;
        }
        v[i*3+i]=1.0;
        g=rv1[i];
        ii=i;
    }
    for (i=n;i>=0;i--) {
        ii=i+1;
        g=w[i];
        if (i < n)
            for (j=ii;j<=n;j++) a[i*3+j]=0.0;
        if (g) {
            g=(SVD_FLOAT)(1.0)/g;
            if (i != n) {
                for (j=ii;j<=n;j++) {
                    for (s=0.0,k=ii;k<=m;k++) s += a[k*3+i]*a[k*3+j];
                    f=(s/a[i*3+i])*g;
                    for (k=i;k<=m;k++) a[k*3+j] += f*a[k*3+i];
                }
            }
            for (j=i;j<=m;j++) a[j*3+i] *= g;
        } else {
            for (j=i;j<=m;j++) a[j*3+i]=0.0;
        }
        ++a[i*3+i];
    }
    for (k=n;k>=0;k--) {
        for (its=1;its<=30;its++) {
            flag=1;
            for (ii=k;ii>=0;ii--) {
                nm=ii-1;
                if (fabs(rv1[ii])+anorm == anorm) {
                    flag=0;
                    break;
                }
                if (fabs(w[nm])+anorm == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=ii;i<=k;i++) {
                    f=s*rv1[i];
                    if (fabs(f)+anorm != anorm) {
                        g=w[i];
                        h=PYTHAG(f,g);
                        w[i]=h;
                        h=(SVD_FLOAT)(1.0)/h;
                        c=g*h;
                        s=(-f*h);
                        for (j=0;j<=m;j++) {
                            y=a[j*3+nm];
                            z=a[j*3+i];
                            a[j*3+nm]=y*c+z*s;
                            a[j*3+i]=z*c-y*s;
                        }
                    }
                }
            }
            z=w[k];
            if (ii == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<=n;j++) v[j*3+k]=(-v[j*3+k]);
                }
                break;
            }
            if (its == 30) return -2; /*No convergence in 30 SVDCMP iterations*/
            x=w[ii];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/((SVD_FLOAT)(2.0)*h*y);
            g=PYTHAG(f,(SVD_FLOAT)(1.0));
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=(SVD_FLOAT)(1.0);
            for (j=ii;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=PYTHAG(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y=y*c;
                for (jj=0;jj<=n;jj++) {
                    x=v[jj*3+j];
                    z=v[jj*3+i];
                    v[jj*3+j]=x*c+z*s;
                    v[jj*3+i]=z*c-x*s;
                }
                z=PYTHAG(f,h);
                w[j]=z;
                if (z) {
                    z=(SVD_FLOAT)(1.0)/z;
                    c=f*z;
                    s=h*z;
                }
                f=(c*g)+(s*y);
                x=(c*y)-(s*g);
                for (jj=0;jj<=m;jj++) {
                    y=a[jj*3+j];
                    z=a[jj*3+i];
                    a[jj*3+j]=y*c+z*s;
                    a[jj*3+i]=z*c-y*s;
                }
            }
            rv1[ii]=(SVD_FLOAT)(0.0);
            rv1[k]=f;
            w[k]=x;
        }
    }

    return 0;
}

#undef SIGN
#undef MAX
#undef PYTHAG