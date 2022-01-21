from numpy import ones,copy,cos,tan,pi,linspace

def gaus_x_w(N):

    # approximation to roots of Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # roots using Newton's method
    eps = 1e-15
    delta = 1.0
    while delta > eps:
        p0 = ones(N,float)
        p1 = copy(x)
        for i in range(1,N):
            p0,p1 = p1,((2*i+1)*x*p1-i*p0)/(i+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

#function being intgrated
def f(x):
    return x**4 - 2*x + 1

N = 3
a = 0.0
b = 2.0

# map points to domain
x,w = gaus_x_w(N)
xp = 0.5*(b-a)*x + 0.5*(b+a)
wp = 0.5*(b-a)*w

# perform integration
s = 0.0
for i in range(N):
    s += wp[i]*f(xp[i])

print(s)

